#!/usr/bin/env python3
"""
ProSeek - Proteome Seeker (taxonomy-aware)
=========================================

Pipeline:
  1. Resolve each taxon name to a current scientific name + taxid + synonyms.
     UniProt taxonomy is the PRIMARY resolver (same namespace the proteomes are
     filed under, and it resolves old binomials NCBI esearch cannot find); NCBI
     is an OPTIONAL secondary for names UniProt misses. On-disk caching.
  2. Look the resolved name/taxid up in the UniProt **reference** proteome table
     (all domains of life - not limited to bacteria, no keyword filtering).
  3. Download each matched reference proteome FASTA into its own file, then
     build a single combined FASTA.
  4. Produce a QC report (HTML + CSVs) covering name replacements, matching,
     and per-organism download counts.

Matching tries taxid -> current name -> synonyms, and compares the genus+species
of the *leading* binomial only (ignoring parentheticals), so "Citrobacter
freundii" no longer matches "Citrobacter rodentium (Citrobacter freundii
biotype)".

Usage:
    python proseek.py taxa.csv --output_dir outputs          # UniProt + NCBI
    python proseek.py taxa.csv --no_ncbi                      # UniProt only
"""

import argparse
import asyncio
import datetime
import gzip
import math
import os
import re
import time
from io import StringIO

import aiohttp
import pandas as pd
import requests
from tqdm import tqdm

try:
    from colorama import Fore, Style, init as colorama_init
    colorama_init()
except Exception:  # colorama is optional
    class _NoColor:
        def __getattr__(self, _):
            return ""
    Fore = Style = _NoColor()


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
# No cap on TOTAL elapsed time: the UniProt reference-proteome catalogue is a
# single tens-of-thousands-row gzipped stream, and some per-organism FASTAs are
# large, so on slow or proxied networks a 120 s total cap can expire mid-read
# (TimeoutError in response.read()). Instead we bound only the connect and the
# per-read STALL: as long as bytes keep arriving the transfer may run as long as
# it needs, but a genuinely dead socket (no data for sock_read s) still bails and
# the _fetch_text retry/backoff kicks in.
TIMEOUT = aiohttp.ClientTimeout(total=None, connect=60, sock_connect=60, sock_read=300)

# Reference proteomes only. UniProt's old 'proteome_type:1' filter now returns
# nothing; the current query term for reference proteomes is 'reference:true'.
# taxonomy_id is added at runtime only if the user restricts to a clade.
PROTEOME_REFERENCE_FILTER = "(reference:true)"
PROTEOME_FIELDS = "upid,organism,organism_id,protein_count,lineage"
# UniProtKB query that selects every sequence belonging to a proteome.
PROTEOME_FASTA_QUERY = "proteome:{upid}"
PROTEOME_BASE = "https://rest.uniprot.org/proteomes/stream"

# UniProt taxonomy search - fallback name resolver for old binomials that NCBI
# esearch cannot find (NCBI only indexes some synonyms with their authority
# string attached, so a bare "prevotella copri" returns nothing).
UNIPROT_TAX_SEARCH = "https://rest.uniprot.org/taxonomy/search"
UNIPROT_TAX_ENTRY = "https://rest.uniprot.org/taxonomy"   # /{taxid}.json
UNIPROT_SPARQL = "https://sparql.uniprot.org/sparql"
UNIPROT_HEADERS = {"User-Agent": "ProSeek/1.0"}

# GTDB / QIIME-style rank prefixes for the metaumbra lineage file. NCBI/UniProt
# put 'Bacteria' at rank superkingdom (or domain); metaumbra's example uses k__,
# so superkingdom/domain/kingdom all map to k__.
RANK_PREFIX = {
    "superkingdom": "k", "domain": "k", "kingdom": "k",
    "phylum": "p", "class": "c", "order": "o",
    "family": "f", "genus": "g", "species": "s",
}
RANK_ORDER = ["k", "p", "c", "o", "f", "g", "s"]

# NCBI lets anonymous clients make ~3 requests/second.
NCBI_DELAY_NO_KEY = 0.34
NCBI_DELAY_WITH_KEY = 0.11


# --------------------------------------------------------------------------- #
# Name-normalisation helpers (the core of the matching fix)
# --------------------------------------------------------------------------- #
def normalize_name(name: str) -> str:
    """Lower-case, strip, collapse whitespace, and drop provisional-genus
    square brackets so '[Clostridium] innocuum' matches 'Clostridium innocuum'."""
    if not isinstance(name, str):
        return ""
    name = name.replace("[", "").replace("]", "")
    return re.sub(r"\s+", " ", name.strip().lower())


def leading_binomial(organism: str) -> str:
    """
    Return the organism name with any parenthetical content removed.

    "Citrobacter rodentium (Citrobacter freundii biotype)" -> "citrobacter rodentium"
    "Escherichia coli (strain K12)"                         -> "escherichia coli"
    """
    main = re.split(r"\s*\(", organism, maxsplit=1)[0]
    return normalize_name(main)


def genus_species(name: str) -> str:
    """First two tokens of a (already de-parenthesised) name."""
    tokens = normalize_name(name).split()
    return " ".join(tokens[:2])


def strip_parenthetical(name: str) -> str:
    """
    Drop a trailing parenthetical annotation, preserving case.

    "Escherichia coli O1:K1 / APEC (APEC O1)" -> "Escherichia coli O1:K1 / APEC"

    Input lists sometimes carry UniProt *proteome* organism names, whose strain
    annotation in parentheses is NOT part of the *taxonomy* scientific name -
    so resolving the bare name is what actually finds the taxid.
    """
    return re.split(r"\s*\(", name, maxsplit=1)[0].strip()


def clean_id(value) -> str:
    """
    Normalise a numeric id (taxid, organism id) to a plain integer string.

    pandas reads id columns that contain blanks as float, so a whole-number
    taxid like 405955 comes back as 405955.0; this strips the trailing '.0'
    and returns '' for anything missing.
    """
    if value is None:
        return ""
    s = str(value).strip()
    if s.lower() in ("", "nan", "none", "<na>"):
        return ""
    try:
        return str(int(float(s)))
    except (ValueError, OverflowError):
        return s


# --------------------------------------------------------------------------- #
# Step 1 - taxonomy resolution (UniProt primary, NCBI optional)
# --------------------------------------------------------------------------- #
def _extract_synonyms(record, scientific_name):
    """Pull plain-name synonyms out of an efetch taxonomy record."""
    other = record.get("OtherNames", {}) or {}
    names = []
    for key in ("Synonym", "GenbankSynonym", "EquivalentName",
                "Includes", "Anamorph", "Teleomorph"):
        val = other.get(key, [])
        if isinstance(val, str):
            val = [val]
        names.extend(val)
    # The typed 'Name' list carries class-tagged names; keep nomenclatural ones.
    for nm in other.get("Name", []):
        try:
            cls = nm.get("ClassCDE", "")
            disp = nm.get("DispName", "")
        except AttributeError:
            continue
        if cls in ("synonym", "equivalent name", "genbank synonym", "includes") and disp:
            names.append(disp)
    # Drop blanks, the scientific name itself, and de-duplicate (case-insensitive).
    seen, out = set(), []
    sci_norm = normalize_name(scientific_name)
    for n in names:
        nn = normalize_name(n)
        if nn and nn != sci_norm and nn not in seen:
            seen.add(nn)
            out.append(n)
    return out


def _sanitize_for_search(name):
    """Make a name safe to send to UniProt's Lucene query parser.

    Strain names contain query operators (':' '/' '(' ')' etc.) that trigger a
    400 Bad Request, so replace them with spaces. Used for BOTH the query sent
    and the token comparison, so the two stay consistent.
    'Escherichia coli O6:H1 (strain CFT073 / ATCC 700928 / UPEC)'
        -> 'escherichia coli o6 h1 strain cft073 atcc 700928 upec'
    """
    s = normalize_name(name)  # lower, strip brackets, collapse whitespace
    s = re.sub(r'[:/()\[\]{}"~*?\\^+=&|!\-]', " ", s)
    return re.sub(r"\s+", " ", s).strip()


def _token_prefix_match(q_tokens, name):
    """True if the query tokens are a leading token-run of `name`.

    Token-level so 'prevotella copri' matches 'Prevotella copri Hayashi et al.
    2007' (authority suffix) but NOT 'prevotella copristata', and 'o1' never
    matches 'o157'. Both sides are sanitized so strain punctuation lines up.
    """
    nt = _sanitize_for_search(name).split()
    return len(nt) >= len(q_tokens) and nt[:len(q_tokens)] == q_tokens


def _names_of_result(r):
    ns = [r.get("scientificName", ""), r.get("commonName", "")]
    ns += r.get("synonyms", []) or []
    on = r.get("otherNames", []) or []
    ns += on if isinstance(on, list) else [on]
    return [x for x in ns if x]


def _uniprot_tax_entry_names(taxid):
    """Full name list for a taxid from the taxonomy entry endpoint (has the
    complete synonym set the search response sometimes abbreviates)."""
    url = f"{UNIPROT_TAX_ENTRY}/{taxid}.json"
    resp = requests.get(url, headers=UNIPROT_HEADERS, timeout=30)
    resp.raise_for_status()
    return _names_of_result(resp.json())


def resolve_via_uniprot_taxonomy(name, retries=3, delay=1.0):
    """
    Primary resolver: UniProt taxonomy search.

    UniProt taxonomy is the same namespace UniProt files proteomes under, so the
    name it returns is - by construction - the name the proteome catalogue uses,
    and it resolves old binomials (e.g. 'prevotella copri') that a bare NCBI
    esearch returns Count=0 for.

    A node is accepted only if the queried name is a leading token-run of one of
    its names - this catches renames whose synonym carries an authority string
    (e.g. 'Prevotella copri Hayashi et al. 2007') while still rejecting unrelated
    fuzzy hits. When the search response abbreviates a node's synonyms, the full
    taxonomy entry is fetched and re-checked before giving up.

    Returns (scientific_name, taxid, synonyms) for the best (species-level) node,
    or None.
    """
    q_tokens = _sanitize_for_search(name).split()
    if not q_tokens:
        return None
    params = {"format": "json", "size": "25", "query": " ".join(q_tokens)}

    def rank_key(r):
        exact = any(_sanitize_for_search(nm).split() == q_tokens
                    for nm in _names_of_result(r))
        return (not exact,
                r.get("rank", "") != "species",
                len(normalize_name(r.get("scientificName", "")).split()))

    for attempt in range(retries):
        try:
            resp = requests.get(UNIPROT_TAX_SEARCH, params=params,
                                headers=UNIPROT_HEADERS, timeout=30)
            resp.raise_for_status()
            results = resp.json().get("results", [])
            if not results:
                return None

            # 1) accept nodes whose returned names token-prefix-match the query
            matches = [r for r in results
                       if any(_token_prefix_match(q_tokens, nm)
                              for nm in _names_of_result(r))]

            # 2) fallback: the search response may abbreviate synonyms - verify
            #    the top species-rank candidates against their full entry.
            if not matches:
                for r in sorted(results, key=rank_key)[:3]:
                    tid = r.get("taxonId")
                    if not tid:
                        continue
                    try:
                        full = _uniprot_tax_entry_names(tid)
                    except Exception:
                        continue
                    if any(_token_prefix_match(q_tokens, nm) for nm in full):
                        matches = [r]
                        break

            if not matches:
                return None

            chosen = sorted(matches, key=rank_key)[0]
            sci = chosen.get("scientificName", "")
            taxid = str(chosen.get("taxonId", "") or "")
            syns = list(chosen.get("synonyms", []) or [])
            return sci, taxid, syns
        except Exception as e:
            print(f"  UniProt taxonomy error for '{name}' (attempt {attempt + 1}): {e}")
            if attempt + 1 == retries:
                return None
            time.sleep(delay * (attempt + 1))


# Biopython is only needed for the optional NCBI secondary; import it lazily so
# a UniProt-only run has no hard dependency on it.
_ENTREZ = None
_ENTREZ_TRIED = False


def _get_entrez(email=None, api_key=None):
    global _ENTREZ, _ENTREZ_TRIED
    if _ENTREZ_TRIED:
        return _ENTREZ
    _ENTREZ_TRIED = True
    try:
        from Bio import Entrez
        if email:
            Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        _ENTREZ = Entrez
    except Exception as e:
        print(f"  Biopython not available; NCBI secondary disabled ({e}).")
        _ENTREZ = None
    return _ENTREZ


def resolve_via_ncbi(name, retries=3, delay=2):
    """
    Optional secondary resolver: NCBI esearch -> efetch.

    Returns (scientific_name, taxid, synonyms) or None. Used only for names the
    UniProt taxonomy search could not resolve (e.g. very recent records UniProt's
    release has not ingested yet).
    """
    Entrez = _get_entrez()
    if Entrez is None:
        return None

    taxon_id = None
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="taxonomy", term=name)
            record = Entrez.read(handle)
            handle.close()
            if record["IdList"]:
                taxon_id = record["IdList"][0]
            break  # empty list is persistent for NCBI - don't retry
        except Exception as err:
            print(f"  NCBI esearch error for '{name}' (attempt {attempt + 1}): {err}")
            if attempt + 1 == retries:
                return None
            time.sleep(delay * (2 ** attempt))

    if taxon_id is None:
        return None

    for attempt in range(retries):
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            rec = records[0]
            resolved = rec["ScientificName"]
            synonyms = _extract_synonyms(rec, resolved)
            return resolved, str(taxon_id), synonyms
        except Exception as err:
            print(f"  NCBI efetch error for '{name}' (attempt {attempt + 1}): {err}")
            if attempt + 1 == retries:
                return None
            time.sleep(delay * (2 ** attempt))


def resolve_taxonomy(taxon_name, use_ncbi=True):
    """
    Resolve a name to (status, resolved_name, taxid, synonyms, resolver).

    status:   'Valid' | 'Renamed' | 'No records found'
    resolver: 'uniprot' | 'ncbi' | ''

    Order: UniProt taxonomy (primary) -> NCBI (optional secondary).
    """
    resolver = ""
    res = None
    matched_variant = taxon_name

    # Try the name as given, then with any trailing parenthetical stripped
    # (strain annotations break exact taxonomy matching).
    variants = [taxon_name]
    stripped = strip_parenthetical(taxon_name)
    if normalize_name(stripped) != normalize_name(taxon_name) and stripped:
        variants.append(stripped)

    for variant in variants:
        res = resolve_via_uniprot_taxonomy(variant)
        if res:
            resolver, matched_variant = "uniprot", variant
            break
    if not res and use_ncbi:
        for variant in variants:
            res = resolve_via_ncbi(variant)
            if res:
                resolver, matched_variant = "ncbi", variant
                break

    if not res:
        return "No records found", None, None, [], ""

    sci, taxid, syns = res
    # Guarantee the queried name is available to the synonym match tier, in case
    # the proteome catalogue still files the organism under that older name.
    for extra in (taxon_name, stripped):
        if extra and normalize_name(extra) not in {normalize_name(s) for s in syns}:
            syns = syns + [extra]
    status = ("Valid"
              if normalize_name(sci) == normalize_name(matched_variant)
              else "Renamed")
    return status, sci, taxid, syns, resolver


def validate_taxa(taxa_list, cache_file, request_delay, use_ncbi=True):
    """
    Validate every taxon name, using an on-disk cache to avoid re-querying.

    Returns a DataFrame with columns:
        original_taxa, status, resolved_name, taxid, synonyms, resolver
    The resolved (current) name, taxid and synonyms are all carried forward so
    the matcher can try taxid first, then the current name, then synonyms.
    'resolver' records whether UniProt or NCBI produced the answer.
    """
    cache = {}
    if cache_file and os.path.exists(cache_file):
        try:
            cdf = pd.read_csv(cache_file).fillna("")
            for _, r in cdf.iterrows():
                cache[r["original_taxa"]] = (
                    r["status"],
                    r["resolved_name"] or None,
                    clean_id(r["taxid"]) or None,
                    [s for s in str(r["synonyms"]).split("|") if s],
                    r.get("resolver", "") or "",
                )
            print(f"Loaded {len(cache)} cached taxonomy lookups.")
        except Exception as e:
            print(f"Could not read cache ({e}); starting fresh.")

    rows = []
    n_uniprot = n_ncbi = 0
    for taxon in tqdm(taxa_list, desc="Validating names", unit="taxon"):
        if taxon in cache:
            status, resolved, taxid, synonyms, resolver = cache[taxon]
        else:
            status, resolved, taxid, synonyms, resolver = resolve_taxonomy(
                taxon, use_ncbi=use_ncbi
            )
            cache[taxon] = (status, resolved, taxid, synonyms, resolver)
            time.sleep(request_delay)

        if resolver == "uniprot":
            n_uniprot += 1
        elif resolver == "ncbi":
            n_ncbi += 1

        rows.append(
            {
                "original_taxa": taxon,
                "status": status,
                "resolved_name": resolved if resolved else "",
                "taxid": clean_id(taxid),
                "synonyms": "|".join(synonyms) if synonyms else "",
                "resolver": resolver,
            }
        )

    tax_df = pd.DataFrame(rows)
    print(f"Resolved by UniProt: {n_uniprot}; by NCBI secondary: {n_ncbi}.")

    if cache_file:
        cache_out = pd.DataFrame(
            [
                {
                    "original_taxa": k,
                    "status": v[0],
                    "resolved_name": v[1] if v[1] else "",
                    "taxid": clean_id(v[2]),
                    "synonyms": "|".join(v[3]) if v[3] else "",
                    "resolver": v[4],
                }
                for k, v in cache.items()
            ]
        )
        cache_out.to_csv(cache_file, index=False)

    return tax_df


# --------------------------------------------------------------------------- #
# Step 2 - fetch UniProt reference proteome metadata
# --------------------------------------------------------------------------- #
async def _fetch_text(session, url, retries=3):
    for attempt in range(retries):
        try:
            async with session.get(url) as response:
                response.raise_for_status()
                compressed = await response.read()
                return gzip.decompress(compressed).decode("utf-8")
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            # many of these (TimeoutError, CancelledError) stringify to "",
            # so show the type too -- otherwise the log line is blank.
            print(f"  Fetch attempt {attempt + 1} failed: {type(e).__name__}: {e}")
            if attempt == retries - 1:
                raise
            await asyncio.sleep(2 ** attempt)


async def _fetch_reference_proteomes_async(taxonomy_id=None):
    query = PROTEOME_REFERENCE_FILTER
    if taxonomy_id:
        query = f"(taxonomy_id:{taxonomy_id}) AND {query}"

    from urllib.parse import urlencode
    params = {
        "compressed": "true",
        "fields": PROTEOME_FIELDS,
        "format": "tsv",
        "query": query,
    }
    url = f"{PROTEOME_BASE}?{urlencode(params)}"

    # trust_env=True so aiohttp honours HTTP(S)_PROXY / NO_PROXY from the
    # environment, the way requests does. On networks that only reach UniProt
    # through a proxy (e.g. HPC compute nodes), aiohttp would otherwise attempt
    # a direct connection that the firewall drops (CancelledError/TimeoutError).
    async with aiohttp.ClientSession(timeout=TIMEOUT, trust_env=True) as session:
        print("Downloading reference proteome catalogue from UniProt...")
        data = await _fetch_text(session, url)

    df = pd.read_csv(StringIO(data), sep="\t").dropna(subset=["Organism"])
    # Normalise the taxid column to clean integer strings (no trailing '.0').
    if "Organism Id" in df.columns:
        df["Organism Id"] = df["Organism Id"].map(clean_id)
    # Protein count is a whole number; keep it integer (nullable) so it never
    # renders as e.g. 2800.0 in the output CSVs.
    if "Protein count" in df.columns:
        df["Protein count"] = pd.to_numeric(
            df["Protein count"], errors="coerce"
        ).astype("Int64")
    print(f"{df.shape[0]} reference proteome entries retrieved.")

    # Guard: the reference filter returning ~nothing means UniProt changed the
    # query syntax again - fail loudly instead of silently matching nothing.
    # Only enforced for the full catalogue; a --taxonomy_id clade can be small.
    threshold = 1000 if taxonomy_id is None else 1
    if df.shape[0] < threshold:
        raise RuntimeError(
            f"Reference proteome catalogue came back with only {df.shape[0]} "
            f"rows using query {query!r}. UniProt likely changed the reference "
            f"filter again (expected tens of thousands). Check "
            f"PROTEOME_REFERENCE_FILTER."
        )
    return df


def fetch_reference_proteomes(taxonomy_id=None):
    return asyncio.run(_fetch_reference_proteomes_async(taxonomy_id))


# --------------------------------------------------------------------------- #
# Step 3 - matching (fixed)
# --------------------------------------------------------------------------- #
def build_match_index(ref_df):
    """Pre-compute normalised genus+species, genus and taxid for every proteome."""
    ref = ref_df.copy()
    ref["_main"] = ref["Organism"].apply(leading_binomial)
    ref["_gs"] = ref["_main"].apply(genus_species)
    ref["_genus"] = ref["_main"].apply(lambda n: n.split(" ")[0] if n else "")
    if "Protein count" in ref.columns:
        ref["Protein count"] = pd.to_numeric(ref["Protein count"], errors="coerce")
    # Normalise the taxid column to a plain integer string for exact comparison.
    if "Organism Id" in ref.columns:
        ref["_taxid"] = ref["Organism Id"].map(clean_id)
    else:
        ref["_taxid"] = ""
    return ref


def _finalise(candidates, matched_form, match_method, match_level, q_norm=None):
    """Pick the best proteome from a candidate set and build the result dict."""
    if candidates.empty:
        return None
    exact = pd.DataFrame()
    if q_norm is not None:
        exact = candidates[candidates["_main"] == q_norm]
    pool = exact if not exact.empty else candidates
    pool = pool.sort_values("Protein count", ascending=False, na_position="last")
    chosen = pool.iloc[0]
    return {
        "Queried Term": matched_form,
        "Match Method": match_method,
        "Match Level": match_level,
        "Matched Organism": chosen["Organism"],
        "Proteome Id": str(chosen["Proteome Id"]).strip(),
        "Organism Id": clean_id(chosen.get("Organism Id", "")),
        "Protein count": chosen["Protein count"],
        "Exact Name Match": bool(not exact.empty),
        "Candidate Count": int(len(candidates)),
        "All Candidate Ids": ";".join(candidates["Proteome Id"].astype(str).tolist()),
    }


def match_by_taxid(taxid, ref_indexed):
    """Match a proteome whose Organism Id equals the resolved NCBI taxid."""
    if not taxid:
        return None
    candidates = ref_indexed[ref_indexed["_taxid"] == str(taxid)]
    return _finalise(candidates, str(taxid), "taxid", "taxid")


def match_by_name(name, ref_indexed):
    """
    Match a proteome by name, at the resolution the query implies.

    - genus-only / species-level (<=2 tokens): match on genus+species and pick
      the most complete proteome (highest protein count); ambiguity is reported
      via Candidate Count.
    - strain / serotype-level (>2 tokens, e.g. 'Escherichia coli O1:K1 / APEC'):
      require the proteome's leading binomial to begin with the SAME token
      sequence. This is token-level (so 'o1' never matches 'o157') and never
      collapses a strain query onto a different strain of the same species. If
      nothing matches at that resolution it returns None - better left unmatched
      and flagged than matched to the wrong strain.

    The query's own trailing parenthetical is stripped first, since input lists
    sometimes carry UniProt proteome names like '... / APEC (APEC O1)'.
    """
    q_main = leading_binomial(name)          # parenthetical stripped + normalised
    q_tokens = q_main.split()
    if not q_tokens:
        return None

    # ---- genus / species level ----
    if len(q_tokens) <= 2:
        if len(q_tokens) == 2:
            candidates = ref_indexed[ref_indexed["_gs"] == q_main]
            level = "genus+species"
        else:
            candidates = ref_indexed[ref_indexed["_genus"] == q_tokens[0]]
            level = "genus"
        return _finalise(candidates, name, None, level, q_norm=q_main)

    # ---- strain / serotype level ----
    pool = ref_indexed[ref_indexed["_gs"] == " ".join(q_tokens[:2])]
    if pool.empty:
        return None

    def token_prefix(main):
        ct = main.split()
        return len(ct) >= len(q_tokens) and ct[:len(q_tokens)] == q_tokens

    cand = pool[pool["_main"].apply(token_prefix)]
    if cand.empty:
        return None

    # Closest first (exact-length match before longer), then most complete.
    cand = cand.assign(_ntok=cand["_main"].apply(lambda m: len(m.split())))
    cand = cand.sort_values(["_ntok", "Protein count"],
                            ascending=[True, False], na_position="last")
    chosen = cand.iloc[0]
    exact = len(chosen["_main"].split()) == len(q_tokens)
    return {
        "Queried Term": name,
        "Match Method": None,
        "Match Level": "strain",
        "Matched Organism": chosen["Organism"],
        "Proteome Id": str(chosen["Proteome Id"]).strip(),
        "Organism Id": clean_id(chosen.get("Organism Id", "")),
        "Protein count": chosen["Protein count"],
        "Exact Name Match": bool(exact),
        "Candidate Count": int(len(cand)),
        "All Candidate Ids": ";".join(cand["Proteome Id"].astype(str).tolist()),
    }


def resolve_match(row, ref_indexed):
    """
    Layered matching for a single taxon, in priority order:
        1. resolved taxid       -> UniProt Organism Id
        2. current valid name   -> genus+species
        3. each synonym         -> genus+species
        4. original input name  -> genus+species (when no resolver found it)

    Returns the match dict (with 'Match Method' set) or None.
    """
    taxid = str(row.get("taxid", "") or "")
    resolved = row.get("resolved_name", "") or ""
    synonyms = [s for s in str(row.get("synonyms", "") or "").split("|") if s]
    original = row["original_taxa"]

    # 1. taxid
    m = match_by_taxid(taxid, ref_indexed)
    if m:
        return m

    # 2. current valid name
    if resolved:
        m = match_by_name(resolved, ref_indexed)
        if m:
            m["Match Method"] = "current-name"
            return m

    # 3. synonyms (covers the case where UniProt still uses the old name)
    for syn in synonyms:
        m = match_by_name(syn, ref_indexed)
        if m:
            m["Match Method"] = "synonym"
            return m

    # 4. fall back to the original input (NCBI gave nothing usable)
    if not resolved:
        m = match_by_name(original, ref_indexed)
        if m:
            m["Match Method"] = "original-name"
            return m

    return None


def find_all_matches(tax_df, ref_df, output_dir):
    ref_indexed = build_match_index(ref_df)

    matched_rows = []
    unmatched = []

    for _, row in tax_df.iterrows():
        result = resolve_match(row, ref_indexed)
        if result is None:
            unmatched.append(
                {
                    "original_taxa": row["original_taxa"],
                    "resolved_name": row.get("resolved_name", ""),
                    "taxid": row.get("taxid", ""),
                    "taxonomy_status": row["status"],
                }
            )
        else:
            result["original_taxa"] = row["original_taxa"]
            result["taxonomy_status"] = row["status"]
            matched_rows.append(result)

    matched_df = pd.DataFrame(matched_rows)

    os.makedirs(output_dir, exist_ok=True)
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    # --- de-duplicate proteomes so each is downloaded exactly once ---
    n_dupes = 0
    if not matched_df.empty:
        matched_df["Proteome Id"] = matched_df["Proteome Id"].astype(str).str.strip()
        dup_mask = matched_df.duplicated(subset=["Proteome Id"], keep="first")
        n_dupes = int(dup_mask.sum())
        if n_dupes:
            # Record which input taxa collapsed onto the same proteome.
            collisions = (
                matched_df.groupby("Proteome Id")
                .filter(lambda g: len(g) > 1)
                .groupby(["Proteome Id", "Matched Organism"])["original_taxa"]
                .apply(lambda s: "; ".join(sorted(set(s))))
                .reset_index()
                .rename(columns={"original_taxa": "input_taxa_sharing_this_proteome"})
            )
            collisions.to_csv(
                os.path.join(reports_dir, "duplicate_matches.csv"), index=False
            )
        matched_df = matched_df[~dup_mask].reset_index(drop=True)

    # --- attach taxonomic lineage from the catalogue (by Proteome Id) ---
    if not matched_df.empty and "Taxonomic lineage" in ref_df.columns:
        lineage_map = (
            ref_df.assign(_pid=ref_df["Proteome Id"].astype(str).str.strip())
            .drop_duplicates(subset=["_pid"])
            .set_index("_pid")["Taxonomic lineage"]
        )
        matched_df["Taxonomic lineage"] = (
            matched_df["Proteome Id"].map(lineage_map).fillna("")
        )

    if not matched_df.empty:
        matched_df.to_csv(os.path.join(reports_dir, "matches.csv"), index=False)
    if unmatched:
        pd.DataFrame(unmatched).to_csv(
            os.path.join(reports_dir, "no_matches.csv"), index=False
        )

    print(f"{len(matched_df)} matched, {len(unmatched)} unmatched.")
    if n_dupes:
        print(f"{Fore.YELLOW}{n_dupes} input taxon/taxa collapsed onto already-"
              f"matched proteomes - de-duplicated; see duplicate_matches.csv."
              f"{Style.RESET_ALL}")
    if not matched_df.empty:
        by_method = matched_df["Match Method"].value_counts().to_dict()
        print("  Match methods: " +
              ", ".join(f"{k}={v}" for k, v in by_method.items()))
        ambiguous = matched_df[matched_df["Candidate Count"] > 1]
        if not ambiguous.empty:
            print(f"{Fore.YELLOW}{len(ambiguous)} match(es) had multiple candidate "
                  f"proteomes - flagged in the QC report.{Style.RESET_ALL}")

    return matched_df, pd.DataFrame(unmatched)


# --------------------------------------------------------------------------- #
# Step 4 - download FASTAs (individual files + combined)
# --------------------------------------------------------------------------- #
def safe_filename(name):
    name = re.sub(r"[^\w\-.]+", "_", str(name)).strip("_")
    return name[:120] if name else "unknown"


def _parse_next_link(link_header):
    """Extract the rel="next" cursor URL from a UniProt Link header."""
    for part in (link_header or "").split(","):
        seg = part.split(";")
        if len(seg) >= 2 and 'rel="next"' in seg[1]:
            return seg[0].strip().strip("<>")
    return None


def _retry_after_seconds(hdr, default=60.0, cap=120.0):
    """Parse a Retry-After header (seconds form) -> float, clamped to cap."""
    try:
        return min(float(str(hdr).strip()), cap) if hdr else default
    except (TypeError, ValueError):
        return default


async def _fetch_proteome_fasta(session, proteome_id, sem, retries=5):
    """
    Fetch all FASTA sequences for a proteome via the paginated search endpoint.

    The uniprotkb/stream endpoint is unreliable (it can return HTTP 200 with the
    body 'Error encountered when streaming data...'), so we use uniprotkb/search
    and follow the Link: rel="next" cursor 500 records at a time.

    HTTP 429 (rate limit) is handled separately: UniProt asks us to wait (its
    Retry-After header, ~60s), so we sleep that long and retry WITHOUT spending a
    normal retry; only genuine errors consume the retry budget. `sem` caps how
    many of these requests run at once across all proteomes.
    """
    from urllib.parse import urlencode
    q = PROTEOME_FASTA_QUERY.format(upid=proteome_id)
    url = ("https://rest.uniprot.org/uniprotkb/search?"
           + urlencode({"format": "fasta", "query": q, "size": "500"}))

    chunks = []
    while url:
        attempt = 0
        throttles = 0
        while True:
            status = text = link = retry_after = None
            err_wait = None
            async with sem:
                try:
                    async with session.get(url) as resp:
                        status = resp.status
                        link = resp.headers.get("Link", "")
                        retry_after = resp.headers.get("Retry-After")
                        text = await resp.text()
                except (aiohttp.ClientError, asyncio.TimeoutError):
                    attempt += 1
                    if attempt >= retries:
                        raise
                    err_wait = 2 ** attempt
            # sleep OUTSIDE the semaphore so a backing-off task frees a slot
            if err_wait is not None:
                await asyncio.sleep(err_wait)
                continue
            if status == 429:                       # rate-limited: wait as asked
                throttles += 1
                if throttles > 20:                  # safety valve
                    raise aiohttp.ClientError("persistent 429 rate limiting")
                await asyncio.sleep(_retry_after_seconds(retry_after))
                continue                            # retry same url, no attempt spent
            # UniProt sometimes answers 200 with an error/notice body -> retry.
            bad = (status != 200
                   or "Error encountered" in text
                   or (text.strip() and not text.lstrip().startswith(">")))
            if bad:
                attempt += 1
                if attempt >= retries:
                    raise aiohttp.ClientError(
                        f"bad response (status {status}): {text[:80].strip()!r}")
                await asyncio.sleep(2 ** attempt)
                continue
            chunks.append(text)
            url = _parse_next_link(link)
            break
    return "".join(chunks)


async def _download_one(session, row, individual_dir, sem, retries=5):
    proteome_id = row["Proteome Id"]
    organism = row["Matched Organism"]
    query = row["Queried Term"]
    fname = f"{safe_filename(organism)}_{safe_filename(proteome_id)}.fasta"
    fpath = os.path.join(individual_dir, fname)

    fail = {
        "Queried Term": query,
        "Matched Organism": organism,
        "Proteome Id": proteome_id,
        "Expected Protein count": row.get("Protein count"),
        "Downloaded Entries": 0,
        "File": "",
        "Success": False,
    }

    try:
        text = await _fetch_proteome_fasta(session, proteome_id, sem, retries=retries)
    except (aiohttp.ClientError, asyncio.TimeoutError) as e:
        print(f"  Download failed for '{query}' ({proteome_id}): {e}")
        return fail

    # Count FASTA records by header lines (lines starting with '>'), NOT raw
    # '>' characters: UniProt descriptions can contain a literal '>'.
    n_entries = sum(1 for ln in text.splitlines() if ln.startswith(">"))
    if n_entries == 0:
        fail["Note"] = "empty response (0 sequences)"
        return fail

    with open(fpath, "w") as fh:
        fh.write(text)
    return {
        "Queried Term": query,
        "Matched Organism": organism,
        "Proteome Id": proteome_id,
        "Expected Protein count": row.get("Protein count"),
        "Downloaded Entries": n_entries,
        "File": fname,
        "Success": True,
    }


async def _download_all_async(matched_df, individual_dir, batch_size,
                              progress=None, max_concurrent=4):
    results = []
    seen = set()  # defensive: never fetch the same proteome twice in a run
    total = len(matched_df)
    done = {"n": 0}
    # cap requests in flight across ALL proteomes so UniProt doesn't 429 us.
    sem = asyncio.Semaphore(max(1, max_concurrent))

    def _tick(bar):
        bar.update()
        done["n"] += 1
        if progress:
            try:
                progress(done["n"], total)
            except Exception:
                pass

    # trust_env=True: route the per-proteome FASTA downloads through the same
    # proxy settings (see _fetch_reference_proteomes_async).
    async with aiohttp.ClientSession(timeout=TIMEOUT, trust_env=True) as session:
        bar = tqdm(total=total, desc="Downloading proteomes", unit="proteome")
        if progress:
            progress(0, total)
        n_batches = math.ceil(total / batch_size) if total else 0
        for b in range(n_batches):
            batch = matched_df.iloc[b * batch_size:(b + 1) * batch_size]
            tasks = []
            for _, row in batch.iterrows():
                pid = str(row["Proteome Id"]).strip()
                if pid in seen:
                    _tick(bar)
                    continue
                seen.add(pid)
                t = asyncio.create_task(
                    _download_one(session, row, individual_dir, sem))
                t.add_done_callback(lambda _: _tick(bar))
                tasks.append(t)
            if tasks:
                results.extend(await asyncio.gather(*tasks))
            await asyncio.sleep(1)  # be polite to UniProt
        bar.close()
    return results


def download_fastas(matched_df, output_dir, batch_size=25, progress=None,
                    max_concurrent=4):
    individual_dir = os.path.join(output_dir, "individual_fastas")
    os.makedirs(individual_dir, exist_ok=True)

    if matched_df.empty:
        print("Nothing matched - skipping download.")
        return pd.DataFrame(), None

    print(f"Downloading {len(matched_df)} reference proteomes "
          f"(<= {max_concurrent} requests in flight)...")
    results = asyncio.run(
        _download_all_async(matched_df, individual_dir, batch_size, progress,
                            max_concurrent=max_concurrent))
    dl_df = pd.DataFrame(results)

    # Build combined FASTA from the individual files that succeeded.
    combined_path = os.path.join(output_dir, "proseek_combined.fasta")
    with open(combined_path, "w") as out:
        for _, r in dl_df[dl_df["Success"]].iterrows():
            fp = os.path.join(individual_dir, r["File"])
            if os.path.exists(fp):
                with open(fp) as src:
                    out.write(src.read())

    return dl_df, combined_path


# --------------------------------------------------------------------------- #
# Step 4b - metaumbra-compatible lineage file
# --------------------------------------------------------------------------- #
def _gtdb_from_rankmap(rankmap):
    """Build 'k__..;p__..;..;s__..' from a {prefix: name} map, present ranks only."""
    return ";".join(f"{r}__{rankmap[r]}" for r in RANK_ORDER if rankmap.get(r))


def _ranked_lineage_uniprot(taxid, retries=3, delay=1.0):
    """Ranked lineage ({rank prefix: name}) for a taxid from UniProt's taxonomy entry."""
    url = f"{UNIPROT_TAX_ENTRY}/{taxid}.json"
    last = None
    for attempt in range(retries):
        try:
            resp = requests.get(url, headers=UNIPROT_HEADERS, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            rankmap = {}
            entries = list(data.get("lineage", []) or [])
            # The organism node itself carries the most specific rank.
            entries.append({"scientificName": data.get("scientificName", ""),
                            "rank": data.get("rank", "")})
            for e in entries:
                rank = (e.get("rank") or "").lower()
                nm = e.get("scientificName", "")
                pref = RANK_PREFIX.get(rank)
                if pref and nm:
                    rankmap.setdefault(pref, nm)
            return rankmap
        except Exception as e:
            last = e
            if attempt + 1 < retries:
                time.sleep(delay * (attempt + 1))
    raise last if last else RuntimeError("lineage fetch failed")


def _ranked_lineage_ncbi(taxid, retries=3, delay=1.0):
    """Ranked lineage ({rank prefix: name}) for a taxid from NCBI efetch (LineageEx)."""
    Entrez = _get_entrez()
    if Entrez is None:
        return {}
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
            rec = Entrez.read(handle)[0]
            handle.close()
            rankmap = {}
            for node in rec.get("LineageEx", []):
                rank = (node.get("Rank") or "").lower()
                nm = node.get("ScientificName", "")
                pref = RANK_PREFIX.get(rank)
                if pref and nm:
                    rankmap.setdefault(pref, nm)
            pref = RANK_PREFIX.get((rec.get("Rank") or "").lower())
            if pref and rec.get("ScientificName"):
                rankmap.setdefault(pref, rec["ScientificName"])
            return rankmap
        except Exception:
            if attempt + 1 < retries:
                time.sleep(delay * (attempt + 1))
    return {}


def _lineage_for_taxid(taxid, cache, use_ncbi):
    """Resolve a GTDB-style lineage string for a taxid (UniProt primary, NCBI fallback)."""
    if taxid in cache:
        return cache[taxid]
    rankmap = {}
    try:
        rankmap = _ranked_lineage_uniprot(taxid)
    except Exception as e:
        print(f"  UniProt lineage error for taxid {taxid}: {e}")
    if use_ncbi and not rankmap:
        try:
            rankmap = _ranked_lineage_ncbi(taxid)
        except Exception as e:
            print(f"  NCBI lineage error for taxid {taxid}: {e}")
    result = _gtdb_from_rankmap(rankmap)
    cache[taxid] = result
    time.sleep(0.1)  # be polite - these fire once per unique taxid
    return result


def write_metaumbra_lineage(matched_df, dl_df, reports_dir, use_ncbi=True):
    """
    Write reports/lineage.tsv in metaumbra format:

        genome_id <TAB> Lineage
        <fasta stem> k__Bacteria;p__...;g__...;s__...

    genome_id is the individual FASTA filename stem, so it matches the digest
    TSV filename stems metaumbra produces from those FASTAs. Lineage is the
    ranked lineage of the proteome's organism taxid.
    """
    if matched_df.empty or dl_df.empty:
        return None
    ok = dl_df[dl_df["Success"]].copy()
    if ok.empty:
        return None

    pid_to_taxid = (
        matched_df.assign(_pid=matched_df["Proteome Id"].astype(str).str.strip())
        .drop_duplicates(subset=["_pid"])
        .set_index("_pid")["Organism Id"]
        .to_dict()
    )

    cache, rows = {}, []
    for _, r in tqdm(ok.iterrows(), total=len(ok),
                     desc="Lineages", unit="genome"):
        stem = os.path.splitext(str(r["File"]))[0]
        taxid = clean_id(pid_to_taxid.get(str(r["Proteome Id"]).strip(), ""))
        lineage = _lineage_for_taxid(taxid, cache, use_ncbi) if taxid else ""
        rows.append({"genome_id": stem, "Lineage": lineage, "taxid": taxid})

    out = pd.DataFrame(rows)

    # metaumbra-strict file: genome_id + Lineage only.
    path = os.path.join(reports_dir, "lineage.tsv")
    out[["genome_id", "Lineage"]].to_csv(path, sep="\t", index=False)

    n_missing = int((out["Lineage"] == "").sum())
    if n_missing:
        print(f"{Fore.YELLOW}{n_missing} genome(s) have an empty lineage "
              f"(taxid lookup failed) - see lineage.tsv.{Style.RESET_ALL}")
    return path


# --------------------------------------------------------------------------- #
# Step 4c - proteome similarity (SPARQL) for multi-candidate taxa
# --------------------------------------------------------------------------- #
def fetch_proteome_similarity(upid, retries=3, delay=1.0):
    """
    Similarity of a proteome to its related (reference) proteomes.

    UniProt's REST API does not expose this, but the SPARQL endpoint does, via
    up:relatedProteome / up:relatedProteomeSimilarityScore. Scores come back as
    fractions (0.909) and are returned here as percentages (90.9).

    Returns a list of (related_upid, similarity_pct) sorted high to low.
    """
    query = f"""PREFIX up: <http://purl.uniprot.org/core/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
SELECT ?proteome ?similarity WHERE {{
    <http://purl.uniprot.org/proteomes/{upid}> up:relatedProteome ?proteome .
    OPTIONAL {{
        ?statement rdf:subject <http://purl.uniprot.org/proteomes/{upid}> ;
            rdf:predicate up:relatedProteome ; rdf:object ?proteome ;
            up:relatedProteomeSimilarityScore ?similarity .
    }}
}} ORDER BY DESC(?similarity)"""
    headers = {"Accept": "application/sparql-results+json", **UNIPROT_HEADERS}
    for attempt in range(retries):
        try:
            resp = requests.get(UNIPROT_SPARQL,
                                params={"query": query, "format": "json"},
                                headers=headers, timeout=60)
            resp.raise_for_status()
            out = []
            for b in resp.json().get("results", {}).get("bindings", []):
                rel = b["proteome"]["value"].rsplit("/", 1)[-1]
                sim = b.get("similarity", {}).get("value")
                pct = round(float(sim) * 100, 1) if sim not in (None, "") else None
                out.append((rel, pct))
            return out
        except Exception as e:
            if attempt + 1 == retries:
                print(f"  Similarity lookup failed for {upid}: {e}")
                return []
            time.sleep(delay * (attempt + 1))


def write_similarity(matched_df, reports_dir):
    """
    For taxa that matched more than one candidate proteome, fetch the chosen
    proteome's similarity to its related reference proteomes and write
    reports/proteome_similarity.csv. Helps judge whether the auto-picked
    proteome is the right representative.
    """
    if matched_df.empty or "Candidate Count" not in matched_df.columns:
        return None
    multi = matched_df[matched_df["Candidate Count"] > 1]
    if multi.empty:
        return None

    rows = []
    for _, r in tqdm(multi.iterrows(), total=len(multi),
                     desc="Similarity", unit="taxon"):
        upid = str(r["Proteome Id"]).strip()
        sims = fetch_proteome_similarity(upid)
        if not sims:
            rows.append({"original_taxa": r.get("original_taxa", ""),
                         "chosen_proteome": upid, "related_proteome": "",
                         "similarity_pct": ""})
        for rel, pct in sims:
            rows.append({"original_taxa": r.get("original_taxa", ""),
                         "chosen_proteome": upid,
                         "related_proteome": rel,
                         "similarity_pct": pct})
        time.sleep(0.2)

    out = pd.DataFrame(rows)
    path = os.path.join(reports_dir, "proteome_similarity.csv")
    out.to_csv(path, index=False)
    print(f"Similarity scores for {len(multi)} multi-candidate taxon(a) "
          f"-> proteome_similarity.csv")
    return path


# --------------------------------------------------------------------------- #
# Step 5 - QC report (HTML + CSVs)
# --------------------------------------------------------------------------- #
def build_qc(tax_df, matched_df, unmatched_df, dl_df, output_dir, meta):
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    # ---- per-organism count validation ----
    if not dl_df.empty:
        dl_df = dl_df.copy()
        dl_df["Expected Protein count"] = pd.to_numeric(
            dl_df["Expected Protein count"], errors="coerce"
        )
        dl_df["Counts Match"] = dl_df.apply(
            lambda r: "pass"
            if r["Success"] and r["Downloaded Entries"] == r["Expected Protein count"]
            else ("fail" if not r["Success"] else "mismatch"),
            axis=1,
        )
        dl_df.to_csv(os.path.join(reports_dir, "download_qc.csv"), index=False)

    # ---- taxonomic lineage of downloaded taxa (separate CSV) ----
    if not matched_df.empty and "Taxonomic lineage" in matched_df.columns:
        lineage_df = matched_df.copy()
        if not dl_df.empty:
            ok = set(dl_df.loc[dl_df["Success"], "Proteome Id"].astype(str))
            lineage_df = lineage_df[lineage_df["Proteome Id"].astype(str).isin(ok)]
        cols = ["original_taxa", "Queried Term", "Matched Organism",
                "Proteome Id", "Organism Id", "Taxonomic lineage"]
        lineage_df = lineage_df[[c for c in cols if c in lineage_df.columns]]
        lineage_df = lineage_df.rename(columns={"Organism Id": "taxid"})
        lineage_df.to_csv(os.path.join(reports_dir, "lineage_detail.csv"), index=False)

    # ---- name-replacement table ----
    tax_df.to_csv(os.path.join(reports_dir, "taxonomy_check.csv"), index=False)
    renamed = tax_df[tax_df["status"] == "Renamed"]

    # ---- metrics ----
    n_total = len(tax_df)
    n_valid = (tax_df["status"] == "Valid").sum()
    n_renamed = (tax_df["status"] == "Renamed").sum()
    n_norec = (tax_df["status"] == "No records found").sum()
    _resolver = tax_df.get("resolver", pd.Series(dtype=str))
    n_by_uniprot = int((_resolver == "uniprot").sum())
    n_by_ncbi = int((_resolver == "ncbi").sum())
    n_matched = len(matched_df)
    n_unmatched = len(unmatched_df)
    n_ambig = int((matched_df["Candidate Count"] > 1).sum()) if not matched_df.empty else 0
    n_dl_ok = int(dl_df["Success"].sum()) if not dl_df.empty else 0
    n_dl_fail = int((~dl_df["Success"]).sum()) if not dl_df.empty else 0
    n_seqs = int(dl_df["Downloaded Entries"].sum()) if not dl_df.empty else 0
    n_count_pass = int((dl_df["Counts Match"] == "pass").sum()) if not dl_df.empty else 0
    if not matched_df.empty:
        mm = matched_df["Match Method"].value_counts()
        n_by_taxid = int(mm.get("taxid", 0))
        n_by_name = int(mm.get("current-name", 0) + mm.get("original-name", 0))
        n_by_syn = int(mm.get("synonym", 0))
    else:
        n_by_taxid = n_by_name = n_by_syn = 0

    metrics = {
        "Input taxa": n_total,
        "Valid names": n_valid,
        "Renamed names": n_renamed,
        "  resolved by UniProt (primary)": n_by_uniprot,
        "  resolved by NCBI (secondary)": n_by_ncbi,
        "No record (any source)": n_norec,
        "Matched to reference proteome": n_matched,
        "  via taxid": n_by_taxid,
        "  via name": n_by_name,
        "  via synonym": n_by_syn,
        "Unmatched": n_unmatched,
        "Ambiguous matches (multi-candidate)": n_ambig,
        "Downloads succeeded": n_dl_ok,
        "Downloads failed": n_dl_fail,
        "Count check passed (entries == expected)": n_count_pass,
        "Total sequences in combined FASTA": n_seqs,
    }

    # Console summary
    print(f"\n{Fore.YELLOW}===== QC SUMMARY ====={Style.RESET_ALL}")
    for k, v in metrics.items():
        print(f"  {k:<45}: {v}")

    # ---- HTML ----
    def table(df, cols=None):
        if df is None or df.empty:
            return "<p class='muted'>None.</p>"
        if cols:
            df = df[[c for c in cols if c in df.columns]]
        return df.to_html(index=False, escape=True, border=0, classes="tbl")

    cards = "".join(
        f"<div class='card'><div class='num'>{v}</div>"
        f"<div class='lbl'>{k}</div></div>"
        for k, v in metrics.items()
    )

    flagged = pd.DataFrame()
    if not matched_df.empty:
        flagged = matched_df[
            (matched_df["Candidate Count"] > 1) | (~matched_df["Exact Name Match"])
        ]

    mismatches = pd.DataFrame()
    if not dl_df.empty:
        mismatches = dl_df[dl_df["Counts Match"] != "pass"]

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>ProSeek QC Report</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
   margin:0;padding:32px;color:#1a1a2e;background:#f6f7fb;}}
 h1{{margin:0 0 4px;font-size:26px;}}
 h2{{margin-top:40px;font-size:18px;border-bottom:2px solid #e0e3ee;padding-bottom:6px;}}
 .sub{{color:#666;margin-bottom:24px;font-size:13px;}}
 .cards{{display:flex;flex-wrap:wrap;gap:12px;}}
 .card{{background:#fff;border:1px solid #e0e3ee;border-radius:10px;
   padding:14px 18px;min-width:140px;flex:1 1 160px;}}
 .num{{font-size:28px;font-weight:700;color:#3b3b98;}}
 .lbl{{font-size:12px;color:#666;margin-top:4px;}}
 table.tbl{{border-collapse:collapse;width:100%;background:#fff;font-size:13px;
   border:1px solid #e0e3ee;border-radius:8px;overflow:hidden;}}
 table.tbl th{{background:#eef0f8;text-align:left;padding:8px 10px;}}
 table.tbl td{{padding:7px 10px;border-top:1px solid #eef0f8;}}
 .muted{{color:#999;}}
 .warn{{color:#b94a00;}}
</style></head><body>
<h1>ProSeek QC Report</h1>
<div class="sub">Generated {meta['timestamp']} &middot; input: {meta['input_file']}
 &middot; reference proteomes only &middot; runtime {meta['runtime']:.2f} min</div>

<h2>Summary</h2>
<div class="cards">{cards}</div>

<h2>Name validation &amp; replacements ({n_renamed} renamed)</h2>
{table(renamed, ['original_taxa','resolved_name','taxid','resolver','status'])}

<h2 class="warn">Flagged matches (multiple candidates or non-exact)</h2>
{table(flagged, ['original_taxa','Queried Term','Match Method','Matched Organism',
 'Proteome Id','Match Level','Exact Name Match','Candidate Count','All Candidate Ids'])}

<h2 class="warn">Count mismatches / failed downloads</h2>
{table(mismatches, ['Queried Term','Matched Organism','Proteome Id',
 'Expected Protein count','Downloaded Entries','Counts Match'])}

<h2>Unmatched taxa</h2>
{table(unmatched_df, ['original_taxa','resolved_name','taxid','taxonomy_status'])}

<h2>Full match table</h2>
{table(matched_df, ['original_taxa','Queried Term','Match Method','Matched Organism',
 'Proteome Id','Organism Id','Protein count','Match Level','Exact Name Match',
 'Candidate Count'])}

<h2>Full download table</h2>
{table(dl_df, ['Queried Term','Matched Organism','Proteome Id',
 'Expected Protein count','Downloaded Entries','Counts Match','File'])}
</body></html>"""

    html_path = os.path.join(reports_dir, "qc_report.html")
    with open(html_path, "w") as fh:
        fh.write(html)
    print(f"\nQC report written to {html_path}")
    return html_path


# --------------------------------------------------------------------------- #
# Input reading
# --------------------------------------------------------------------------- #
def read_taxa(input_file):
    ext = os.path.splitext(input_file)[1].lower()
    delimiter = "," if ext == ".csv" else "\t" if ext == ".tsv" else None
    if delimiter is None:
        raise ValueError("Unsupported file format - provide a .csv or .tsv file.")

    # Files exported from Excel/web are often not UTF-8 (e.g. Windows-1252 with
    # non-breaking spaces, byte 0xa0). Try encodings in order; latin-1 always
    # succeeds as a last resort since it maps every byte.
    df = None
    last_err = None
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            df = pd.read_csv(input_file, sep=delimiter, encoding=enc)
            break
        except UnicodeDecodeError as e:
            last_err = e
    if df is None:
        raise last_err

    if "taxa" not in df.columns:
        raise ValueError("Input file must contain a 'taxa' column.")

    rows = df.shape[0]
    print(f"{input_file} loaded: {rows} rows.")
    taxa = (
        df["taxa"].dropna().astype(str)
        # Replace non-breaking spaces / stray whitespace, then collapse.
        .str.replace("\xa0", " ", regex=False)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip().str.lower()
        .replace("", pd.NA).dropna()
        .drop_duplicates().tolist()
    )
    dupes = rows - len(taxa)
    print(f"{dupes} duplicate(s) removed." if dupes > 0 else "No duplicates found.")
    return taxa


# --------------------------------------------------------------------------- #
# Args + main
# --------------------------------------------------------------------------- #
def parse_arguments():
    p = argparse.ArgumentParser(
        description="Resolve taxon names (UniProt, optional NCBI), match them "
                    "to UniProt reference proteomes, download FASTAs and QC "
                    "the result."
    )
    p.add_argument("input_file", help="CSV/TSV with a 'taxa' column")
    p.add_argument("--output_dir", default="outputs", help="Output directory")
    p.add_argument("--no_ncbi", action="store_true",
                   help="Disable the optional NCBI secondary resolver "
                        "(UniProt-only; no Biopython needed)")
    p.add_argument("--email", default="anonymous@example.com",
                   help="Email for the NCBI secondary resolver (NCBI policy)")
    p.add_argument("--api_key", default=None, help="NCBI Entrez API key (optional)")
    p.add_argument("--taxonomy_id", default=None,
                   help="Restrict UniProt to a clade taxid (default: all domains)")
    p.add_argument("--batch_size", type=int, default=25,
                   help="Proteome download batch size")
    p.add_argument("--max_concurrent", type=int, default=4,
                   help="Max UniProt requests in flight at once (default 4). "
                        "Lower this if UniProt returns HTTP 429 rate-limit "
                        "errors; raise it cautiously for faster downloads.")
    p.add_argument("--no_similarity", action="store_true",
                   help="Skip SPARQL proteome-similarity lookup for "
                        "multi-candidate taxa")
    p.add_argument("--include_all", action="store_true",
                   help="For a taxon matching several proteomes, keep ALL of "
                        "them (species + every subspecies/strain) instead of "
                        "one auto-pick")
    return p.parse_args()


def candidate_details(matched_df, ref_df, with_similarity=True):
    """
    For taxa that matched more than one candidate proteome, return
        {original_taxa: [ {proteome_id, organism, protein_count, taxid,
                           lineage, similarity, chosen}, ... ]}
    so a UI can let the user pick which proteome to use. `similarity` is the
    chosen proteome's SPARQL similarity to that candidate where available.
    """
    if matched_df.empty or "Candidate Count" not in matched_df.columns:
        return {}
    ref = None
    if "Proteome Id" in ref_df.columns:
        ref = ref_df.copy()
        ref["_pid"] = ref["Proteome Id"].astype(str).str.strip()
        ref = ref.drop_duplicates("_pid").set_index("_pid")

    out = {}
    for _, r in matched_df[matched_df["Candidate Count"] > 1].iterrows():
        ids = [x for x in str(r.get("All Candidate Ids", "")).split(";") if x]
        chosen = str(r["Proteome Id"]).strip()
        sims = {}
        if with_similarity:
            sims = {rel: pct for rel, pct in fetch_proteome_similarity(chosen)}
        cands = []
        for pid in ids:
            row = ref.loc[pid] if (ref is not None and pid in ref.index) else None
            pc = None
            if row is not None and pd.notna(row.get("Protein count")):
                pc = int(row["Protein count"])
            cands.append({
                "proteome_id": pid,
                "organism": row["Organism"] if row is not None else pid,
                "protein_count": pc,
                "taxid": clean_id(row["Organism Id"]) if row is not None else "",
                "lineage": (row.get("Taxonomic lineage", "")
                            if row is not None else ""),
                "similarity": sims.get(pid),
                "chosen": pid == chosen,
            })
        cands.sort(key=lambda c: (not c["chosen"], -(c["protein_count"] or 0)))
        out[r["original_taxa"]] = cands
    return out


def apply_selections(matched_df, candidates, selections):
    """
    Apply per-taxon proteome choices.

    selections maps original_taxa -> a proteome id OR a list of ids. When a
    taxon maps to several ids (e.g. 'Bifidobacterium longum' kept alongside its
    'subsp. longum' and 'subsp. infantis' proteomes), that taxon's single row is
    replaced by one row per selected proteome - so a species query can pull all
    of its subspecies/strain proteomes at once.
    """
    df = matched_df.copy()
    selections = selections or {}
    kept = []
    for _, row in df.iterrows():
        taxa = row["original_taxa"]
        if taxa not in selections:
            kept.append(row.to_dict())
            continue
        pids = selections[taxa]
        if isinstance(pids, str):
            pids = [pids]
        cmap = {c["proteome_id"]: c for c in candidates.get(taxa, [])}
        for pid in pids:
            pick = cmap.get(pid)
            if not pick:
                continue
            nr = row.to_dict()
            nr["Proteome Id"] = pick["proteome_id"]
            nr["Matched Organism"] = pick["organism"]
            nr["Organism Id"] = pick["taxid"]
            nr["Protein count"] = pick["protein_count"]
            if "Taxonomic lineage" in nr:
                nr["Taxonomic lineage"] = pick["lineage"]
            kept.append(nr)
    out = pd.DataFrame(kept) if kept else matched_df.iloc[0:0].copy()
    return out.drop_duplicates(subset=["Proteome Id"], keep="first")


def stage_match(input_file, output_dir="outputs", use_ncbi=True,
                taxonomy_id=None, email="anonymous@example.com", api_key=None,
                with_similarity=True):
    """Phase 1 for the UI: read -> resolve -> fetch catalogue -> match.
    Returns everything needed to render the pick screen and then finalise."""
    os.makedirs(output_dir, exist_ok=True)
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    request_delay = 0.15
    if use_ncbi:
        _get_entrez(email=email, api_key=api_key)
        request_delay = NCBI_DELAY_WITH_KEY if api_key else NCBI_DELAY_NO_KEY

    taxa = read_taxa(input_file)
    tax_df = validate_taxa(taxa, os.path.join(reports_dir, "taxonomy_cache.csv"),
                           request_delay, use_ncbi=use_ncbi)
    ref_df = fetch_reference_proteomes(taxonomy_id)
    matched_df, unmatched_df = find_all_matches(tax_df, ref_df, output_dir)
    candidates = candidate_details(matched_df, ref_df, with_similarity)

    return {
        "tax_df": tax_df, "matched_df": matched_df, "unmatched_df": unmatched_df,
        "candidates": candidates,
        "matched": len(matched_df), "unmatched": len(unmatched_df),
        "multi": len(candidates),
        "input_name": os.path.basename(input_file),
    }


def stage_download(state, selections, output_dir="outputs", use_ncbi=True,
                   batch_size=25, similarity=True, progress=None):
    """Phase 2 for the UI: apply the user's picks, download, build reports."""
    start = time.time()
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    matched_df = apply_selections(state["matched_df"], state["candidates"],
                                  selections)
    dl_df, combined_path = download_fastas(matched_df, output_dir,
                                           batch_size=batch_size,
                                           progress=progress)
    meta = {
        "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "input_file": state.get("input_name", ""),
        "runtime": (time.time() - start) / 60,
    }
    build_qc(state["tax_df"], matched_df, state["unmatched_df"], dl_df,
             output_dir, meta)
    write_metaumbra_lineage(matched_df, dl_df, reports_dir, use_ncbi=use_ncbi)
    if similarity:
        write_similarity(matched_df, reports_dir)

    report_path = os.path.join(reports_dir, "qc_report.html")
    return {
        "output_dir": output_dir,
        "combined_fasta": combined_path,
        "reports_dir": reports_dir,
        "report_html": report_path if os.path.exists(report_path) else None,
        "downloaded": int(dl_df["Success"].sum()) if not dl_df.empty else 0,
        "matched": len(matched_df),
        "runtime_min": (time.time() - start) / 60,
    }


def run_proseek(input_file, output_dir="outputs", use_ncbi=True,
                taxonomy_id=None, batch_size=25, email="anonymous@example.com",
                api_key=None, similarity=True, include_all=False,
                max_concurrent=4):
    """
    Run the full ProSeek pipeline. Importable so a UI (e.g. the Shiny app) can
    call it directly. Writes all outputs under output_dir and returns a dict of
    the key result paths.

    include_all: when True, a taxon that matched several candidate proteomes
    keeps ALL of them (species + every subspecies/strain proteome) instead of
    the single auto-pick.
    """
    start = time.time()
    os.makedirs(output_dir, exist_ok=True)
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    request_delay = 0.15
    if use_ncbi:
        _get_entrez(email=email, api_key=api_key)
        request_delay = NCBI_DELAY_WITH_KEY if api_key else NCBI_DELAY_NO_KEY

    print(f"{Fore.YELLOW}** ProSeek: taxonomy-aware proteome seeker **{Style.RESET_ALL}")

    print(f"{Fore.BLUE}Reading taxa...{Style.RESET_ALL}")
    taxa = read_taxa(input_file)

    resolver_note = "UniProt primary + NCBI secondary" if use_ncbi else "UniProt only"
    print(f"{Fore.BLUE}Resolving names ({resolver_note})...{Style.RESET_ALL}")
    cache_file = os.path.join(reports_dir, "taxonomy_cache.csv")
    tax_df = validate_taxa(taxa, cache_file, request_delay, use_ncbi=use_ncbi)

    print(f"{Fore.BLUE}Fetching UniProt reference proteomes...{Style.RESET_ALL}")
    ref_df = fetch_reference_proteomes(taxonomy_id)

    print(f"{Fore.BLUE}Matching names to proteomes...{Style.RESET_ALL}")
    matched_df, unmatched_df = find_all_matches(tax_df, ref_df, output_dir)

    if include_all and not matched_df.empty:
        cands = candidate_details(matched_df, ref_df, with_similarity=False)
        if cands:
            selections = {t: [c["proteome_id"] for c in lst]
                          for t, lst in cands.items()}
            before = len(matched_df)
            matched_df = apply_selections(matched_df, cands, selections)
            print(f"  include_all: expanded {len(cands)} multi-candidate taxa "
                  f"-> {len(matched_df) - before} extra proteome(s).")

    print(f"{Fore.BLUE}Downloading FASTAs...{Style.RESET_ALL}")
    dl_df, combined_path = download_fastas(matched_df, output_dir,
                                           batch_size=batch_size,
                                           max_concurrent=max_concurrent)

    print(f"{Fore.BLUE}Building QC report...{Style.RESET_ALL}")
    meta = {
        "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "input_file": os.path.basename(input_file),
        "runtime": (time.time() - start) / 60,
    }
    build_qc(tax_df, matched_df, unmatched_df, dl_df, output_dir, meta)

    print(f"{Fore.BLUE}Writing metaumbra lineage file...{Style.RESET_ALL}")
    lineage_path = write_metaumbra_lineage(matched_df, dl_df, reports_dir,
                                           use_ncbi=use_ncbi)

    similarity_path = None
    if similarity:
        print(f"{Fore.BLUE}Fetching proteome similarity (multi-candidate)...{Style.RESET_ALL}")
        similarity_path = write_similarity(matched_df, reports_dir)

    runtime = (time.time() - start) / 60
    print(f"\n{Fore.GREEN}Done. Proteomes acquired!{Style.RESET_ALL}")
    print(f"  Combined FASTA : {combined_path}")
    print(f"  Individual FASTAs: {os.path.join(output_dir, 'individual_fastas')}")
    if lineage_path:
        print(f"  Lineage (metaumbra): {lineage_path}")
    print(f"  Reports         : {reports_dir}")
    print(f"{Fore.YELLOW}Total time: {runtime:.2f} min{Style.RESET_ALL}")

    return {
        "output_dir": output_dir,
        "combined_fasta": combined_path,
        "individual_fastas": os.path.join(output_dir, "individual_fastas"),
        "reports_dir": reports_dir,
        "lineage": lineage_path,
        "similarity": similarity_path,
        "matched": len(matched_df),
        "unmatched": len(unmatched_df),
        "downloaded": int(dl_df["Success"].sum()) if not dl_df.empty else 0,
        "runtime_min": runtime,
    }


def main():
    args = parse_arguments()
    run_proseek(
        input_file=args.input_file,
        output_dir=args.output_dir,
        use_ncbi=not args.no_ncbi,
        taxonomy_id=args.taxonomy_id,
        batch_size=args.batch_size,
        email=args.email,
        api_key=args.api_key,
        similarity=not args.no_similarity,
        include_all=args.include_all,
        max_concurrent=args.max_concurrent,
    )


if __name__ == "__main__":
    main()

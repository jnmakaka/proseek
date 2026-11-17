# ProSeek

Taxonomy-aware downloader for UniProt reference proteomes, for metaproteomics workflows.

Give it a list of taxon names. ProSeek resolves current names and taxids (UniProt, with optional NCBI fallback), matches them to UniProt reference proteomes, downloads each proteome FASTA plus a combined FASTA, and writes a QC report and a lineage file.

## Install
```
pip install pandas aiohttp requests tqdm
# optional: colorama, biopython (NCBI fallback resolver)
```

## Use
```
python proseek.py taxa.csv --output_dir outputs    # taxa.csv needs a 'taxa' column
python proseek.py taxa.csv --no_ncbi               # UniProt only
python proseek.py --help                           # all options
```

## Output
- `individual_fastas/` and `proseek_combined.fasta`
- `reports/`: `qc_report.html`, `matches.csv`, `no_matches.csv`, `taxonomy_check.csv`, `download_qc.csv`, `proteome_similarity.csv`, and `lineage.tsv` (GTDB-style, MetaUmbra-compatible)

It can also be imported (`run_proseek`, `stage_match`, `stage_download`) so a UI can call it.

Proteomes are from UniProt (CC BY 4.0); please cite UniProt.

# BLAST_parser

A short pipeline for BLASTing one or more genes against a set of genomes, then parsing the hits, extracting the sequences and making reports and plots. Adapted from PhaseCaller but not specific to any gene.

It is meant for quickly running BLAST searches for whatever gene you want, checking that the primers hit the gene you expect, and seeing a bit of genomic context (flanking genes). `--gene` with several names (or `all`) just runs several configs on the same BLAST database, it doesn't compare genes.

## Requirements

- Python 3.10+ (tested on 3.11 and 3.13)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (`blastn` and `makeblastdb` need to be on your PATH, tested with 2.16.0)
- Python packages (see `requirements.txt`): biopython, pandas, openpyxl, numpy, matplotlib, seaborn

```
conda env create -f environment.yml
conda activate blast_parser_env
```

or, if you already have BLAST+ installed:

```
pip install -r requirements.txt
```

## Usage

```
python run_pipeline.py --gene <gene_name> [<gene_name> ...] --input <genomes> --output <output_folder>
```

| argument | what it does |
|---|---|
| `--gene` | gene name(s) matching `configs/config_<gene>.py`. `all` runs every config except the template/example |
| `--input` | folder of genome FASTAs (`.fasta`, `.fa`, `.fna`, `.fas`) or a `.zip` of them. One file per isolate, the file name is the isolate ID |
| `--output` | output folder |
| `--keep_blast` | keep the combined genome FASTA and BLAST database (`blast_temp/`) after running |

Example:

```
python run_pipeline.py --gene lic1 --input inputs/test_50 --output outputs
```

## Configs

Make a new config from the template:

```
python create_config.py <gene_name>
```

This copies `configs/config_template.py` to `configs/config_<gene_name>.py`. Each config is for one gene and has:

- `gene_name` and `identifier` (the output sub-folder is named after `identifier` if set, otherwise `gene_name`)
- `gene_queries` - one or more reference sequences for the gene (alternative references hitting the same copy are merged)
- `flank_queries` - flanking genes with `position` (`upstream`/`downstream`) and expected `distance` from the gene (can be negative if they overlap), plus `flank_buffer`
- `gene_start_queries` / `gene_end_queries` - optional, only used as supporting evidence (no boundary refinement)
- `primer_queries` - primers with `position` (`forward`/`reverse`)
- `*_blast_params` - BLAST settings and thresholds (`min_identity`, `min_coverage`, optional `high_identity`/`high_coverage`; primers use `max_mismatches` and `min_alignment_pct`)
- `extraction_params` - flank length (`flanking_bp`), contig boundary cutoffs, `extract_flank` and `extract_igr`

Optional groups can be left as `{}` or with an empty `query_sequence` and they are skipped.

## What the pipeline does

1. Combine all genomes into one FASTA and make a BLAST database
2. Write the query FASTAs and run the BLAST searches (gene, flank, gene start/end, primers)
3. Parse the hits: apply thresholds, merge overlapping alternative references into one locus per gene copy, assign flanks to each gene hit
4. Extract the gene for each locus, and optionally the gene plus `flanking_bp` of flank each side (reverse strand handled, flank clipped at contig edges), check contig boundaries, and remove loci with an N
5. Write reports and plots

The extracted gene is the region covered by the BLAST hit, so it can be a partial gene.

## Outputs

Everything goes in `<output>/<identifier or gene_name>/`:

```
blast/       query FASTAs and raw BLAST tables
sequences/   extracted FASTAs
reports/     Excel and text reports
plots/       PNGs
```

`<output>/run_pipeline.log` has the log, and `<output>/blast_temp/` is removed at the end unless `--keep_blast` is used.

### Sequences

- `<gene>_gene_all_sequences.fasta` - the gene for every retained locus
- `<gene>_gene_with_flank_all_sequences.fasta` - same loci with `flanking_bp` of flank each side (if `extract_flank`)
- `..._boundary_critical.fasta` / `..._boundary_warning.fasta` - loci close to a contig edge
- `multicopy/` - one FASTA per isolate with more than one copy, plus a combined file
- `igr/` - upstream/downstream intergenic regions (if `extract_igr`)

A locus is dropped if the gene or the flank contains an N, so all FASTAs have the same loci.

### Reports (`reports/`)

- `<gene>_blast_summary.xlsx` - stage counts, per-query summary, isolate x query matrix, hits by query, totals, isolate summary, copy number, gene hits, gene hits with assigned flanks, flank hits, primer hits
- `<gene>_sequence_extraction.xlsx` - one row per locus sent to extraction, including excluded ones (outcome, boundary status, N in gene vs flank only)
- `sequences/pcr_fragments/<gene>_pcr_fragments.fasta` - (only if primers searched and a pair is found) in silico PCR products: forward + reverse primer hit on the same contig, facing each other, within `max_pcr_fragment_length`. Includes the primer sites, written from the forward primer, header has coords, primer differences and gene overlap
- `<gene>_primer_report.xlsx` - (only if primers searched) summary per primer, differences per hit and per position, primers per locus, in silico PCR fragments, PCR fragment length summary, every primer hit
- `<gene>_gene_hit_stats.txt` / `<gene>_flank_hit_stats.txt` - identity/coverage bins per query

Searches that weren't run say "not searched" rather than 0.

### Plots (`plots/`)

Raw BLAST identity/coverage per query, query comparison, query heatmap, flank assignment, copy number, span vs reference length, sequence extraction QC, contig boundary QC, retained locus quality, primer mismatches, primer positions around the gene and in silico PCR fragment lengths.

## Counts

Isolates, BLAST alignment rows and loci are counted separately. The stage summary sheet gives the unit for each stage: raw rows, rows passing thresholds, merged loci, loci after the per-isolate cap, extracted before N filtering, removed by the N filter and retained.

## Not included

Start/end boundary refinement.

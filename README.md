# Paediatric respiratory CRPA genomic analysis code

This repository contains a lightweight, publication-oriented code release for the manuscript:

**Global clonal backbones frame local virulence signals in paediatric respiratory carbapenem-resistant *Pseudomonas aeruginosa***

The scripts reproduce the main tabular analyses used to support matched sequence-type (ST) comparisons, candidate gene-level signals, temporal summaries, unsupervised PCA, and exploratory local-versus-public machine-learning analyses.

## Scope

This code release is intentionally minimal. It does not include patient-level clinical data, raw sequencing reads, genome assemblies, or local file paths. Input tables should be generated from the user's own genome-analysis pipeline or requested from the corresponding author where permitted by ethics and data-sharing approvals.

## Repository Layout

```text
code_release/
  README.md
  requirements.txt
  .gitignore
  data/
    README.md
    example_integrated_table.tsv
  scripts/
    01_compare_matched_st_modules.py
    02_gene_level_fisher_fdr.py
    03_temporal_summary.py
    04_pca_presence_absence.py
    05_ml_local_public.py
```

## Input Tables

Most scripts use an integrated per-genome table with one row per isolate/genome. Required and optional columns are described in `data/README.md`.

At minimum, the integrated table should contain:

- `genome_id`: genome or isolate identifier
- `source_label`: `local` or `public`
- `ST_label`: sequence type label, such as `ST244`
- `T3SS_type`: `ExoS`, `ExoU`, `ExoS+U`, or `T3SS-negative`
- `collection_year`: year for local isolates, if available
- `total_vf_genes`: total detected virulence-associated genes
- `total_resist_genes`: total detected resistance-associated genes
- binary module columns such as `Phenazine`, `T6SS`, `Efflux_MexEF`
- binary gene columns for gene-level analyses

## Example Commands

Create matched-ST module summaries:

```bash
python scripts/01_compare_matched_st_modules.py \
  --input data/integrated_table.tsv \
  --outdir results/module_comparison
```

Run gene-level Fisher exact tests with FDR correction:

```bash
python scripts/02_gene_level_fisher_fdr.py \
  --input data/gene_presence_long.tsv \
  --outdir results/gene_level_stats
```

Summarise temporal patterns:

```bash
python scripts/03_temporal_summary.py \
  --input data/integrated_table.tsv \
  --outdir results/temporal
```

Run PCA on a binary presence-absence matrix:

```bash
python scripts/04_pca_presence_absence.py \
  --input data/integrated_table.tsv \
  --feature-prefix gene__ \
  --outdir results/pca
```

Run exploratory local-versus-public machine-learning models:

```bash
python scripts/05_ml_local_public.py \
  --input data/integrated_table.tsv \
  --outdir results/ml
```

## Software

Tested with Python 3.11. Install dependencies with:

```bash
pip install -r requirements.txt
```

## Data Availability

The code does not include restricted clinical or genomic datasets. De-identified analysis tables supporting the manuscript may be available from the corresponding author upon reasonable request and subject to institutional approvals.

## Code Availability

Custom scripts used for the manuscript analyses are provided here for transparency and reuse. Paths and data inputs are parameterised so users can run the workflow with their own de-identified tables.


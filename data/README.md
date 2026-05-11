# Data Input Format

This folder contains a small schema example only. Do not upload patient-level clinical data, raw sequencing reads, genome assemblies, or restricted metadata.

## Integrated Per-Genome Table

Recommended file name: `integrated_table.tsv`

Required columns:

- `genome_id`: genome or isolate identifier
- `source_label`: `local` or `public`
- `ST_label`: sequence type label, for example `ST244`
- `T3SS_type`: `ExoS`, `ExoU`, `ExoS+U`, or `T3SS-negative`
- `collection_year`: sampling year for local isolates if available; may be empty for public genomes
- `total_vf_genes`: total number of detected virulence-associated genes
- `total_resist_genes`: total number of detected resistance-associated genes

Recommended module columns:

- `T3SS_effectors`
- `T3SS_apparatus`
- `Alginate`
- `Biofilm_Fap`
- `Flagella`
- `Pili_T4`
- `Pyoverdine`
- `Pyochelin`
- `Protease`
- `T6SS`
- `QS_Las`
- `QS_Rhl`
- `Toxins`
- `Phenazine`
- `Carbapenem`
- `Beta_lactam`
- `Aminoglycoside`
- `Fluoroquinolone`
- `Efflux_MexAB`
- `Efflux_MexCD`
- `Efflux_MexEF`
- `Efflux_MexXY`

Binary gene columns can be included either as direct gene names or with a prefix such as `gene__phzA1`.

## Long Gene-Presence Table

Recommended file name: `gene_presence_long.tsv`

Required columns:

- `genome_id`
- `source_label`
- `ST_label`
- `gene`
- `present`
- `kind`: optional, for example `VFDB` or `CARD`


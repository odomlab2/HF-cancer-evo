# Hair follicle cancer evolution - manuscript analysis code

This repository is a manuscript-reproduction snapshot for the submitted
analysis. It contains publication-safe processed tables, analysis code, and a
single offline validation entry point. It does not contain primary BAM files,
generated RDS objects, generated PDFs, or private computing infrastructure.

## Reproduce the submitted statistics

Requirements: R 4.4.2, base R packages, and `sha256sum` on `PATH`. From a fresh
checkout, run exactly:

```sh
Rscript scripts/public/reproduce_submission_snapshot.R public-output/manuscript
```

The command uses only tracked CSV/TXT inputs, performs no network request, and
writes ten validation files beneath the ignored `public-output/` directory. It
refuses to overwrite an existing output directory.

It reproduces the statistical tests used for the manuscript.

## Scope

Primary targeted/exome sequencing data are available through ArrayExpress
accession `E-MTAB-17446`. A full upstream rerun additionally requires a
compatible GRCm38/mm10-plus-PhiX FASTA, matching BAMs and QC exports, GATK,
samtools, Ensembl VEP cache 102 for GRCm38, the reference CDS object used by
dNdScv, and these R packages at the versions recorded in the validated R 4.4.2
environment: `babelgene` 22.9, `cbioportalR` 1.1.1, `dndscv` 0.0.1.0, `dplyr`
1.1.4, `emmeans` 2.0.0, `ggalluvial` 0.12.6, `ggbeeswarm` 0.7.3, `ggh4x`
0.3.1, `gghalves` 0.1.4, `ggpattern` 1.3.1, `ggplot2` 4.0.3, `ggrepel` 0.9.6,
`glmmTMB` 1.1.14, `gridExtra` 2.3, `httr` 1.4.7, `jsonlite` 2.0.0, `MASS`
7.3.65, `msigdbr` 26.1.0, `patchwork` 1.3.1, `purrr` 1.0.2, `readr` 2.1.5,
`sandwich` 3.1.1, `scales` 1.4.0, `stringr` 1.5.1, `survival` 3.8.6,
`survminer` 0.5.2, `tibble` 3.2.1, `tidyr` 1.3.1, `UpSetR` 1.4.0, and `vcfR`
1.15.0. Those inputs and environments are outside the public execution
contract.

## Scientific snapshot status

The frozen external-data tables were used for submission and do not require a
live cBioPortal or orthologue service. Their provenance is summarized in
`docs/public-inputs-and-provenance.md`.

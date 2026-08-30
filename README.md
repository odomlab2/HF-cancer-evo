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

It independently recalculates the Figure 2A Kruskal–Wallis tests, shared
mutation percentages, Figure 4F proportions, and Figure 5G permutation. It also
validates the frozen submitted Supplementary Figure 3C model results, Figure 5C
data invariant, and Figure 5E display aliases. The Supplementary Figure 3C
model is not refitted because its generated pipeline RDS inputs are deliberately
not distributed; its publication-critical fitted result is retained as an
approved frozen statistics table.

The frozen Supplementary Figure 3C table retains the submitted model's
false-convergence warning and convergence code `1` (`pdHess = TRUE`). The
public command validates the submitted odds ratio and BH-adjusted P-value but
does not suppress, reinterpret, or refit that model. The advertised command
itself completes without runtime warnings in the validated environment.

## Scope

The advertised command reproduces or validates the submitted-manuscript
statistics covered by the manuscript harness from included publication-safe
processed inputs. It does not claim BAM-to-manuscript reproducibility and does
not regenerate every manuscript PDF.

Primary targeted/exome sequencing data are available through ArrayExpress
accession `E-MTAB-17446`. A full upstream rerun additionally requires a
compatible GRCm38/mm10-plus-PhiX FASTA, matching BAMs and QC exports, GATK,
samtools, Ensembl VEP cache 102 for GRCm38, the reference CDS object used by
dNdScv, and the R/Python packages imported by individual legacy scripts. Those
inputs and environments are outside the public execution contract.

Optional upstream paths are configured with environment variables listed in
`config/public-inputs.env.example`. Missing required inputs fail explicitly;
no private fallback locations are embedded in the tree. Generated upstream and
figure files must be directed to `HF_SCC_OUTPUT_ROOT`.

## Scientific snapshot status

This tree intentionally preserves the submitted analysis behavior. Read
`docs/known-vcf-parser-limitation.md` before interpreting or extending the VCF
parsing code. A genotype-role correction is not included in this snapshot.

The frozen external-data tables were used for submission and do not require a
live cBioPortal or orthologue service. Their provenance is summarized in
`docs/public-inputs-and-provenance.md`.

## Citation and license

Citation metadata are provided in `CITATION.cff`. No software license was
present in the source repository, so no license has been invented for this
snapshot; the authors should select one before inviting reuse.

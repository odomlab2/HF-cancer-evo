# Release content policy

The public tree applies the following content classification:

| Class | Treatment |
|---|---|
| Private absolute and infrastructure paths | Removed or replaced by environment-variable interfaces |
| Private operational metadata exports | Excluded |
| Primary BAM/QC/reference inputs | External; not tracked |
| Generated RDS, PDF, PNG, and R Markdown reports | Excluded and ignored |
| Internal validation history, logs, and archived scheduler workflows | Excluded |
| Publication-critical statistics and plot-data tables | Retained |
| Public-safe sample annotations, panel definitions, and callable summaries | Retained |
| Credentials, secrets, and tokens | None detected in the release tree |

The R Markdown reports were excluded only after a dependency audit found no
publication-critical calculation unique to them. Some reports contained
exploratory alternatives with different thresholds or model specifications;
those were not migrated into the submitted-analysis path.

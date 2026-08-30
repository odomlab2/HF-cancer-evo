# Known VCF parser limitation

This release reproduces the analysis used for the submitted manuscript.

A post-submission audit found that the historical parser selected a genotype
column by position rather than by biological sample role. Correcting that
selection changes some downstream scientific values. The positional behavior
is therefore preserved here solely to reproduce the submitted analysis; it
must not be interpreted as biologically role-aware or biologically correct.

The role-aware correction will be handled as a separate, explicitly validated
scientific revision. It is intentionally not mixed into this reproduction
snapshot.

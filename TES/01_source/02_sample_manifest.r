tes_sample_stages <- c(
  "mutation_calling", "vep_annotation", "callable_bases", "sequencing_qc",
  "mutation_merge", "sample_annotation", "analysis_dataset",
  "technical_outliers")

tes_sample_manifest <- function(samples_path) {
  sample_ids <- readLines(samples_path, warn = FALSE)
  stages <- as.data.frame(setNames(
    replicate(length(tes_sample_stages), rep(TRUE, length(sample_ids)),
      simplify = FALSE),
    tes_sample_stages), check.names = FALSE)
  manifest <- cbind(data.frame(
    assay = rep("TES", length(sample_ids)),
    sample_name = sample_ids,
    role = rep("scientific_sample", length(sample_ids)),
    stringsAsFactors = FALSE), stages)
  validate_sample_manifest(manifest, tes_sample_stages)
}

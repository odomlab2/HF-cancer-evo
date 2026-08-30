wes_sample_stages <- c(
  "mutation_calling", "vep_annotation", "callable_bases", "sequencing_qc",
  "mutation_merge", "sample_annotation", "analysis_dataset",
  "technical_outliers")

wes_sample_manifest <- function(samples_path, controls_path) {
  sample_ids <- readLines(samples_path, warn = FALSE)
  controls <- utils::read.csv(controls_path, stringsAsFactors = FALSE)
  required <- c("sample_name", "role", "mutation_output_expected",
    "callable_expected")
  if (length(setdiff(required, names(controls))) ||
      anyDuplicated(controls$sample_name) ||
      !all(controls$sample_name %in% sample_ids) ||
      !all(controls$role == "matched_normal") ||
      any(controls$mutation_output_expected) ||
      !all(controls$callable_expected)) {
    stop("Invalid WES matched-normal control manifest.", call. = FALSE)
  }
  analysis <- !sample_ids %in% controls$sample_name
  stages <- data.frame(
    mutation_calling = analysis, vep_annotation = analysis,
    callable_bases = TRUE, sequencing_qc = TRUE,
    mutation_merge = analysis, sample_annotation = analysis,
    analysis_dataset = analysis, technical_outliers = analysis)
  manifest <- cbind(data.frame(
    assay = rep("WES", length(sample_ids)), sample_name = sample_ids,
    role = ifelse(analysis, "scientific_sample", "matched_normal"),
    stringsAsFactors = FALSE), stages)
  validate_sample_manifest(manifest, wes_sample_stages)
}

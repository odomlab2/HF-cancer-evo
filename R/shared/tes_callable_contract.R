# TES-specific machine-readable callable-bases provenance.
tes_callable_contract_fields <- c(
  "assay", "reference_genome", "min_depth", "mapq", "base_quality",
  "target_bed", "target_bed_sha256", "target_size_rule", "target_size",
  "depth_base_quality_option", "depth_mapping_quality_option",
  "depth_exclude_flags", "depth_suppress_overlaps", "depth_count_deletions",
  "on_target_depth_all_positions", "on_target_depth_scope",
  "on_target_reads_scope", "overall_depth_scope", "total_reads_scope",
  "view_exclude_flags", "view_require_flags", "sample_manifest",
  "sample_manifest_sha256", "expected_sample_count", "row_schema")

tes_callable_contract_error <- function(message) {
  stop(paste("Invalid TES callable provenance contract:", message),
    call. = FALSE)
}

read_tes_callable_contract <- function(path, role) {
  if (!file.exists(path)) {
    tes_callable_contract_error(paste("missing", role, "contract"))
  }
  lines <- readLines(path, warn = FALSE)
  if (length(lines) != 2L) {
    tes_callable_contract_error(paste(role,
      "contract must contain one header and one row"))
  }
  header <- strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]]
  values <- strsplit(lines[[2L]], "\t", fixed = TRUE)[[1L]]
  if (!identical(header, tes_callable_contract_fields)) {
    tes_callable_contract_error(paste(role, "contract schema mismatch"))
  }
  if (length(values) != length(header) || any(!nzchar(values))) {
    tes_callable_contract_error(paste(role, "contract has empty fields"))
  }
  names(values) <- header
  numeric <- c("min_depth", "mapq", "base_quality", "target_size",
    "expected_sample_count")
  if (any(!grepl("^[0-9]+$", values[numeric]))) {
    tes_callable_contract_error(paste(role,
      "numeric fields must be non-negative integers"))
  }
  values
}

validate_tes_callable_contract <- function(expected_path, declared_path,
                                            expected_ids, root) {
  expected <- read_tes_callable_contract(expected_path, "authoritative")
  declared <- read_tes_callable_contract(declared_path, "declared run")
  mismatched <- tes_callable_contract_fields[expected != declared]
  if (length(mismatched)) {
    details <- vapply(mismatched, function(field) sprintf(
      "%s expected %s, declared %s", field, expected[[field]],
      declared[[field]]), character(1))
    tes_callable_contract_error(paste("mismatch:", paste(details,
      collapse = "; ")))
  }
  if (expected[["assay"]] != "TES" ||
      as.integer(expected[["expected_sample_count"]]) != length(expected_ids)) {
    tes_callable_contract_error("assay or expected sample count mismatch")
  }
  manifest <- file.path(root, expected[["sample_manifest"]])
  if (!file.exists(manifest)) tes_callable_contract_error("missing sample manifest")
  manifest_hash <- strsplit(system2("sha256sum", shQuote(manifest), stdout = TRUE),
    " ", fixed = TRUE)[[1L]][1L]
  if (manifest_hash != expected[["sample_manifest_sha256"]]) {
    tes_callable_contract_error("sample manifest SHA-256 mismatch")
  }
  bed <- file.path(root, expected[["target_bed"]])
  if (!file.exists(bed)) tes_callable_contract_error("missing target BED")
  hash <- strsplit(system2("sha256sum", shQuote(bed), stdout = TRUE),
    " ", fixed = TRUE)[[1L]][1L]
  if (hash != expected[["target_bed_sha256"]]) {
    tes_callable_contract_error("target BED SHA-256 mismatch")
  }
  regions <- read.table(bed, header = FALSE, sep = "\t")
  size <- sum(regions[[3L]] - regions[[2L]] + 1)
  if (size != as.numeric(expected[["target_size"]])) {
    tes_callable_contract_error("target BED inclusive size mismatch")
  }
  expected
}

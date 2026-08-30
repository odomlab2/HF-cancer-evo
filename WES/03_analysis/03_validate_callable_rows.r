#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("Usage: 03_validate_callable_rows.r HF_SCC_ROOT ROWS_DIR",
    "AUTHORITATIVE_CONTRACT RUN_CONTRACT"), call. = FALSE)
}
root <- normalizePath(args[[1L]], mustWork = TRUE)
rows <- args[[2L]]
source(file.path(root, "R", "shared", "sample_manifest.R"))
source(file.path(root, "R", "shared", "callable_rows.R"))
source(file.path(root, "R", "shared", "callable_contract.R"))
source(file.path(root, "WES", "01_source", "03_sample_manifest.r"))

manifest <- wes_sample_manifest(
  file.path(root, "WES", "02_data", "00_raw", "samples.txt"),
  file.path(root, "WES", "01_source", "wes_matched_normal_controls.csv"))
expected <- select_manifest_stage(manifest, "callable_bases")$sample_name
if (length(expected) != 88L) {
  stop("WES callable manifest must contain exactly 88 samples.", call. = FALSE)
}
contract <- validate_callable_contract(args[[3L]], args[[4L]], expected)
validate_callable_row_directory(rows, expected,
  expected_mapq = as.integer(contract[["mapq"]]),
  expected_base_quality = as.integer(contract[["base_quality"]]))
cat(contract[["callable_column"]])

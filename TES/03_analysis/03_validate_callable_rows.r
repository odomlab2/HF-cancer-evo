#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("Usage: 03_validate_callable_rows.r HF_SCC_ROOT ROWS_DIR",
    "AUTHORITATIVE_CONTRACT RUN_CONTRACT"), call. = FALSE)
}
root <- normalizePath(args[[1L]], mustWork = TRUE)
source(file.path(root, "R", "shared", "sample_manifest.R"))
source(file.path(root, "R", "shared", "tes_callable_contract.R"))
source(file.path(root, "R", "shared", "tes_callable_rows.R"))
source(file.path(root, "TES", "01_source", "02_sample_manifest.r"))

manifest <- tes_sample_manifest(file.path(root, "TES", "02_data", "00_raw",
  "samples.txt"))
expected <- select_manifest_stage(manifest, "callable_bases")$sample_name
if (length(expected) != 131L || sum(grepl("^hairfollicle", expected)) != 99L ||
    sum(grepl("^skin", expected)) != 32L) {
  stop("TES callable manifest must contain 99 hair-follicle and 32 skin samples.",
    call. = FALSE)
}
contract <- validate_tes_callable_contract(args[[3L]], args[[4L]], expected, root)
schema <- strsplit(contract[["row_schema"]], ",", fixed = TRUE)[[1L]]
if (!identical(schema, tes_callable_row_fields)) {
  tes_callable_contract_error("row_schema does not match TES row validator")
}
validate_tes_callable_rows(args[[2L]], expected, contract)
cat(paste(schema, collapse = "\t"))

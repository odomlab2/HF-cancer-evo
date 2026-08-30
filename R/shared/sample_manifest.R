# Pure validation and selection helpers for assay-specific sample manifests.
validate_sample_manifest <- function(
    manifest,
    stage_columns,
    required_columns = c("assay", "sample_name", "role")) {
  if (!is.data.frame(manifest)) {
    stop("Sample manifest must be a data frame.", call. = FALSE)
  }
  required <- c(required_columns, stage_columns)
  missing_columns <- setdiff(required, names(manifest))
  if (length(missing_columns)) {
    stop("Sample manifest is missing columns: ",
      paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  for (column in required_columns) {
    value <- manifest[[column]]
    if (!is.character(value) || anyNA(value) || any(!nzchar(value))) {
      stop("Manifest column '", column,
        "' must contain non-empty character values.", call. = FALSE)
    }
  }
  duplicates <- unique(manifest$sample_name[duplicated(manifest$sample_name)])
  if (length(duplicates)) {
    stop("Duplicate sample IDs: ", paste(duplicates, collapse = ", "),
      call. = FALSE)
  }
  invalid_stages <- stage_columns[!vapply(manifest[stage_columns],
    function(x) is.logical(x) && !anyNA(x), logical(1))]
  if (length(invalid_stages)) {
    stop("Stage columns must be complete logical vectors: ",
      paste(invalid_stages, collapse = ", "), call. = FALSE)
  }
  manifest
}

select_manifest_stage <- function(manifest, stage) {
  validate_sample_manifest(manifest, stage)
  manifest[manifest[[stage]], , drop = FALSE]
}

compare_sample_ids <- function(expected, observed, order_sensitive = FALSE) {
  expected <- as.character(expected)
  observed <- as.character(observed)
  expected_duplicates <- unique(expected[duplicated(expected)])
  observed_duplicates <- unique(observed[duplicated(observed)])
  missing <- setdiff(expected, observed)
  extra <- setdiff(observed, expected)
  order_equal <- identical(expected, observed)
  list(
    pass = !length(expected_duplicates) && !length(observed_duplicates) &&
      !length(missing) && !length(extra) &&
      (!order_sensitive || order_equal),
    expected_duplicates = expected_duplicates,
    observed_duplicates = observed_duplicates,
    missing = missing,
    extra = extra,
    set_equal = !length(missing) && !length(extra),
    order_equal = order_equal)
}

assert_sample_ids <- function(
    expected, observed, order_sensitive = FALSE, label = "sample IDs") {
  comparison <- compare_sample_ids(expected, observed, order_sensitive)
  if (!comparison$pass) {
    stop(label, " mismatch; missing: ", paste(comparison$missing, collapse = ", "),
      "; extra: ", paste(comparison$extra, collapse = ", "),
      "; expected duplicates: ",
      paste(comparison$expected_duplicates, collapse = ", "),
      "; observed duplicates: ",
      paste(comparison$observed_duplicates, collapse = ", "),
      "; order equal: ", comparison$order_equal, call. = FALSE)
  }
  invisible(comparison)
}

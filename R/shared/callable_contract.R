# Machine-readable provenance contract for callable-bases row sets.
callable_contract_fields <- c(
  "assay", "min_depth", "mapq", "base_quality",
  "expected_sample_count", "callable_column")

callable_contract_error <- function(message) {
  stop(paste("Invalid WES callable provenance contract:", message),
    call. = FALSE)
}

read_callable_contract <- function(path, role) {
  if (!file.exists(path)) {
    callable_contract_error(paste("missing", role, "contract"))
  }
  lines <- readLines(path, warn = FALSE)
  if (length(lines) != 2L) {
    callable_contract_error(paste(role, "contract must contain one header and one row"))
  }
  header <- strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]]
  values <- strsplit(lines[[2L]], "\t", fixed = TRUE)[[1L]]
  if (!identical(header, callable_contract_fields)) {
    callable_contract_error(paste(role, "contract schema mismatch"))
  }
  if (length(values) != length(header) || any(!nzchar(values))) {
    callable_contract_error(paste(role, "contract must have six non-empty fields"))
  }
  names(values) <- header
  numeric_fields <- c("min_depth", "mapq", "base_quality", "expected_sample_count")
  if (any(!grepl("^[0-9]+$", values[numeric_fields]))) {
    callable_contract_error(paste(role, "numeric fields must be non-negative integers"))
  }
  values
}

validate_callable_contract <- function(expected_path, declared_path,
                                       expected_ids) {
  expected <- read_callable_contract(expected_path, "authoritative")
  declared <- read_callable_contract(declared_path, "declared run")
  mismatched <- callable_contract_fields[expected != declared]
  if (length(mismatched)) {
    details <- vapply(mismatched, function(field) sprintf(
      "%s expected %s, declared %s", field, expected[[field]], declared[[field]]),
      character(1))
    callable_contract_error(paste("mismatch:", paste(details, collapse = "; ")))
  }
  if (as.integer(expected[["expected_sample_count"]]) != length(expected_ids)) {
    callable_contract_error(sprintf(
      "expected_sample_count %s does not match manifest count %s",
      expected[["expected_sample_count"]], length(expected_ids)))
  }
  expected
}

publication_external_data_mode <- function() {
  mode <- tolower(Sys.getenv("HF_SCC_EXTERNAL_DATA_MODE", unset = "frozen"))
  if (!mode %in% c("frozen", "refresh")) {
    stop("HF_SCC_EXTERNAL_DATA_MODE must be 'frozen' or 'refresh'.")
  }
  mode
}

publication_require_live_refresh <- function() {
  if (!identical(publication_external_data_mode(), "refresh")) {
    stop("Live external-data access is disabled in publication mode.")
  }
  publication_output_root()
  invisible(TRUE)
}

publication_repo_input <- function(logical_path, expected_sha256) {
  path <- file.path(publication_repo_root(), logical_path)
  publication_verify_frozen_input(path, expected_sha256)
}

publication_sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(output, " +")[[1L]][[1L]]
}

publication_verify_frozen_input <- function(path, expected_sha256) {
  if (!file.exists(path)) stop("Missing frozen submission input: ", path)
  observed <- publication_sha256(path)
  if (!identical(observed, expected_sha256)) {
    stop("Frozen submission input hash mismatch for ", path,
      ": expected ", expected_sha256, ", observed ", observed)
  }
  invisible(path)
}

publication_copy_frozen <- function(source, logical_target, expected_sha256) {
  publication_verify_frozen_input(source, expected_sha256)
  target <- publication_output_path(logical_target)
  if (!file.copy(source, target, overwrite = TRUE, copy.mode = FALSE)) {
    stop("Failed to copy frozen submission artifact to ", target)
  }
  invisible(target)
}

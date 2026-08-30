# Validation helpers for headerless per-sample callable-bases rows.
callable_row_error <- function(problems) {
  message <- paste("Invalid callable row directory:", paste(problems, collapse = "; "))
  condition <- structure(list(message = message, call = NULL),
    class = c("callable_row_validation_error", "error", "condition"))
  stop(condition)
}

read_callable_row <- function(path) {
  size <- file.info(path)$size
  if (is.na(size) || size == 0) return(list(issue = "empty file", id = NA_character_))
  bytes <- readBin(path, what = "raw", n = size)
  if (tail(bytes, 1L) != as.raw(10) || sum(bytes == as.raw(10)) != 1L) {
    return(list(issue = "expected one newline-terminated row", id = NA_character_))
  }
  if (sum(bytes == as.raw(9)) != 4L) {
    return(list(issue = "schema mismatch: expected five tab-separated fields",
      id = NA_character_))
  }
  text <- tryCatch(rawToChar(head(bytes, -1L)), error = function(e) NA_character_)
  fields <- if (is.na(text)) character() else strsplit(text, "\t", fixed = TRUE)[[1L]]
  if (length(fields) != 5L || any(!nzchar(fields))) {
    return(list(issue = "schema mismatch: empty or undecodable field", id = NA_character_))
  }
  integers <- suppressWarnings(as.integer(fields[2:5]))
  if (any(!grepl("^[0-9]+$", fields[2:5])) || anyNA(integers)) {
    return(list(issue = "schema mismatch: fields 2-5 must be non-negative integers",
      id = fields[[1L]]))
  }
  list(issue = NULL, id = fields[[1L]], mapq = integers[[3L]],
    base_quality = integers[[4L]])
}

validate_callable_row_directory <- function(rows_dir, expected_ids,
                                            expected_mapq = NULL,
                                            expected_base_quality = NULL) {
  if (!dir.exists(rows_dir)) callable_row_error(paste("missing directory", rows_dir))
  entries <- list.files(rows_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  expected_files <- paste0(expected_ids, ".tsv")
  problems <- character()
  info <- file.info(entries)
  non_files <- basename(entries[is.na(info$isdir) | info$isdir])
  if (length(non_files)) problems <- c(problems,
    paste("unexpected non-file entries:", paste(non_files, collapse = ", ")))
  files <- entries[!is.na(info$isdir) & !info$isdir]
  names_found <- basename(files)
  missing_files <- setdiff(expected_files, names_found)
  unexpected_files <- setdiff(names_found, expected_files)
  if (length(missing_files)) problems <- c(problems,
    paste("missing expected files:", paste(missing_files, collapse = ", ")))
  if (length(unexpected_files)) problems <- c(problems,
    paste("unexpected files:", paste(unexpected_files, collapse = ", ")))
  rows <- lapply(files, read_callable_row)
  issues <- vapply(rows, function(x) if (is.null(x$issue)) "" else x$issue,
    character(1))
  if (any(nzchar(issues))) problems <- c(problems, paste0("malformed file ",
    names_found[nzchar(issues)], ": ", issues[nzchar(issues)]))
  observed <- vapply(rows, `[[`, character(1), "id")
  valid <- !is.na(observed)
  mismatched <- valid & names_found != paste0(observed, ".tsv")
  if (any(mismatched)) problems <- c(problems, paste0("filename/sample mismatch ",
    names_found[mismatched], " -> ", observed[mismatched]))
  comparison <- compare_sample_ids(expected_ids, observed[valid])
  if (length(comparison$missing)) problems <- c(problems,
    paste("missing sample IDs:", paste(comparison$missing, collapse = ", ")))
  if (length(comparison$extra)) problems <- c(problems,
    paste("unexpected sample IDs:", paste(comparison$extra, collapse = ", ")))
  for (id in comparison$observed_duplicates) problems <- c(problems,
    paste0("duplicate sample ID ", id, " in files: ",
      paste(names_found[valid][observed[valid] == id], collapse = ", ")))
  row_mapq <- vapply(rows, function(x) if (is.null(x$mapq)) NA_integer_ else x$mapq,
    integer(1))
  row_bq <- vapply(rows, function(x) {
    if (is.null(x$base_quality)) NA_integer_ else x$base_quality
  }, integer(1))
  mapq_mismatch <- valid & !is.na(row_mapq) & row_mapq != expected_mapq
  if (!is.null(expected_mapq) && any(mapq_mismatch)) {
    bad <- names_found[mapq_mismatch]
    problems <- c(problems, paste0("row MAPQ differs from declared contract in: ",
      paste(bad, collapse = ", ")))
  }
  bq_mismatch <- valid & !is.na(row_bq) & row_bq != expected_base_quality
  if (!is.null(expected_base_quality) && any(bq_mismatch)) {
    bad <- names_found[bq_mismatch]
    problems <- c(problems, paste0(
      "row base quality differs from declared contract in: ",
      paste(bad, collapse = ", ")))
  }
  if (length(problems)) callable_row_error(problems)
  invisible(data.frame(file = files, sample_name = observed,
    stringsAsFactors = FALSE))
}

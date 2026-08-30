# Validation helpers for headerless nine-field TES callable-bases rows.
tes_callable_row_fields <- c(
  "sample_name", "total_target_bases", "callable_bases_10x",
  "pct_callable_10x", "on_target_reads", "total_reads_mapq20",
  "overall_callable_10x", "mapq", "bq")

tes_callable_row_error <- function(problems) {
  stop(paste("Invalid TES callable row directory:",
    paste(problems, collapse = "; ")), call. = FALSE)
}

read_tes_callable_row <- function(path) {
  size <- file.info(path)$size
  if (is.na(size) || size == 0) return(list(issue = "empty file", id = NA_character_))
  bytes <- readBin(path, what = "raw", n = size)
  if (tail(bytes, 1L) != as.raw(10) || sum(bytes == as.raw(10)) != 1L) {
    return(list(issue = "expected one newline-terminated row", id = NA_character_))
  }
  if (sum(bytes == as.raw(9)) != 8L) {
    return(list(issue = "schema mismatch: expected nine tab-separated fields",
      id = NA_character_))
  }
  text <- tryCatch(rawToChar(head(bytes, -1L)), error = function(e) NA_character_)
  fields <- if (is.na(text)) character() else strsplit(text, "\t", fixed = TRUE)[[1L]]
  if (length(fields) != 9L || any(!nzchar(fields))) {
    return(list(issue = "schema mismatch: empty or undecodable field",
      id = NA_character_))
  }
  integers <- c(2L, 3L, 5L:9L)
  if (any(!grepl("^[0-9]+$", fields[integers]))) {
    return(list(issue = "schema mismatch: count and quality fields must be integers",
      id = fields[[1L]]))
  }
  if (!grepl("^[0-9]+[.][0-9]{4}$", fields[[4L]])) {
    return(list(issue = "schema mismatch: percentage must have four decimals",
      id = fields[[1L]]))
  }
  list(issue = NULL, id = fields[[1L]], values = fields)
}

validate_tes_callable_rows <- function(rows_dir, expected_ids, contract) {
  if (!dir.exists(rows_dir)) tes_callable_row_error(paste("missing directory", rows_dir))
  entries <- list.files(rows_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  info <- file.info(entries)
  files <- entries[!is.na(info$isdir) & !info$isdir]
  names_found <- basename(files)
  expected_files <- paste0(expected_ids, ".tsv")
  problems <- character()
  non_files <- basename(entries[is.na(info$isdir) | info$isdir])
  if (length(non_files)) problems <- c(problems,
    paste("unexpected non-file entries:", paste(non_files, collapse = ", ")))
  missing <- setdiff(expected_files, names_found)
  extra <- setdiff(names_found, expected_files)
  if (length(missing)) problems <- c(problems,
    paste("missing expected files:", paste(missing, collapse = ", ")))
  if (length(extra)) problems <- c(problems,
    paste("unexpected files:", paste(extra, collapse = ", ")))
  rows <- lapply(files, read_tes_callable_row)
  issues <- vapply(rows, function(x) if (is.null(x$issue)) "" else x$issue,
    character(1))
  if (any(nzchar(issues))) problems <- c(problems, paste0("malformed file ",
    names_found[nzchar(issues)], ": ", issues[nzchar(issues)]))
  ids <- vapply(rows, `[[`, character(1), "id")
  valid <- !is.na(ids) & !nzchar(issues)
  mismatch <- valid & names_found != paste0(ids, ".tsv")
  if (any(mismatch)) problems <- c(problems, paste0("filename/sample mismatch ",
    names_found[mismatch], " -> ", ids[mismatch]))
  comparison <- compare_sample_ids(expected_ids, ids[valid])
  if (length(comparison$missing)) problems <- c(problems,
    paste("missing sample IDs:", paste(comparison$missing, collapse = ", ")))
  if (length(comparison$extra)) problems <- c(problems,
    paste("unexpected sample IDs:", paste(comparison$extra, collapse = ", ")))
  for (id in comparison$observed_duplicates) problems <- c(problems,
    paste0("duplicate sample ID ", id, " in files: ",
      paste(names_found[valid][ids[valid] == id], collapse = ", ")))
  expected_size <- contract[["target_size"]]
  for (i in which(valid)) {
    fields <- rows[[i]]$values
    pct <- sprintf("%.4f", 100 * as.numeric(fields[[3L]]) / as.numeric(fields[[2L]]))
    if (fields[[2L]] != expected_size) problems <- c(problems,
      paste("target size differs from declared contract in", names_found[[i]]))
    if (fields[[4L]] != pct) problems <- c(problems,
      paste("callable percentage mismatch in", names_found[[i]]))
    if (fields[[8L]] != contract[["mapq"]]) problems <- c(problems,
      paste("row MAPQ differs from declared contract in", names_found[[i]]))
    if (fields[[9L]] != contract[["base_quality"]]) problems <- c(problems,
      paste("row base quality differs from declared contract in", names_found[[i]]))
  }
  if (length(problems)) tes_callable_row_error(problems)
  invisible(data.frame(file = files, sample_name = ids, stringsAsFactors = FALSE))
}

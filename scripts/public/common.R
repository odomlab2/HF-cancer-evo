assert_identical <- function(observed, expected, label) {
  if (!identical(observed, expected)) {
    stop(label, " differs: observed=", format(observed, digits = 17),
      ", expected=", format(expected, digits = 17), call. = FALSE)
  }
  invisible(observed)
}

read_publication_csv <- function(root, path) {
  read.csv(file.path(root, path), check.names = FALSE, stringsAsFactors = FALSE)
}

write_result <- function(data, output, name) {
  write.csv(data, file.path(output, name), row.names = FALSE, na = "")
}

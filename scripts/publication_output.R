# Publication writes require an explicit output root disjoint from the repository.
publication_repo_root <- function() {
  configured <- Sys.getenv("HF_SCC_ROOT", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))
  here <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(here, ".git"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) stop("Set HF_SCC_ROOT to the repository root.")
    here <- parent
  }
}

publication_is_within <- function(path, root) {
  identical(path, root) || startsWith(path, paste0(root, .Platform$file.sep))
}

publication_normalize_candidate <- function(path) {
  if (!grepl("^/", path)) path <- file.path(getwd(), path)
  normalizePath(path, mustWork = FALSE)
}

publication_output_root <- function() {
  value <- Sys.getenv("HF_SCC_OUTPUT_ROOT", unset = "")
  if (!nzchar(value)) stop("Set HF_SCC_OUTPUT_ROOT for every publication write.")
  if (!grepl("^/", value)) stop("HF_SCC_OUTPUT_ROOT must be an absolute path.")
  dir.create(value, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(value, mustWork = TRUE)
  repo <- publication_repo_root()
  if (publication_is_within(root, repo) || publication_is_within(repo, root)) {
    stop("HF_SCC_OUTPUT_ROOT must be disjoint from the repository and sealed outputs.")
  }
  root
}

publication_output_path <- function(path) {
  if (is.null(path) || !is.character(path) || length(path) != 1L || !nzchar(path)) {
    return(path)
  }
  repo <- publication_repo_root()
  root <- publication_output_root()
  source_path <- publication_normalize_candidate(path)
  if (publication_is_within(source_path, root)) {
    target <- source_path
  } else if (publication_is_within(source_path, repo)) {
    relative <- substring(source_path, nchar(repo) + 2L)
    target <- file.path(root, relative)
  } else {
    stop("Publication output escapes HF_SCC_ROOT/HF_SCC_OUTPUT_ROOT: ", path)
  }
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  parent <- normalizePath(dirname(target), mustWork = TRUE)
  if (!publication_is_within(parent, root)) stop("Publication output parent escapes root.")
  file.path(parent, basename(target))
}

publication_trace_writer <- function(name, namespace, argument) {
  where <- asNamespace(namespace)
  tracer <- bquote(.(as.name(argument)) <- publication_output_path(.(as.name(argument))))
  invisible(capture.output(trace(name, where = where, tracer = tracer, print = FALSE)))
}

capture.output <- function(..., file = NULL, append = FALSE,
    type = c("output", "message"), split = FALSE) {
  utils::capture.output(..., file = publication_output_path(file), append = append,
    type = type, split = split)
}

if (!isTRUE(getOption("hf_scc.publication_writers_guarded"))) {
  publication_trace_writer("saveRDS", "base", "file")
  publication_trace_writer("save", "base", "file")
  publication_trace_writer("writeLines", "base", "con")
  publication_trace_writer("write.csv", "utils", "file")
  publication_trace_writer("write.table", "utils", "file")
  publication_trace_writer("pdf", "grDevices", "file")
  publication_trace_writer("png", "grDevices", "filename")
  publication_trace_writer("jpeg", "grDevices", "filename")
  publication_trace_writer("tiff", "grDevices", "filename")
  publication_trace_writer("svg", "grDevices", "filename")
  publication_trace_writer("write_csv", "readr", "file")
  publication_trace_writer("write_tsv", "readr", "file")
  publication_trace_writer("ggsave", "ggplot2", "filename")
  if ("arrow" %in% loadedNamespaces()) {
    publication_trace_writer("write_parquet", "arrow", "sink")
  } else {
    setHook(packageEvent("arrow", "onLoad"), function(...) {
      publication_trace_writer("write_parquet", "arrow", "sink")
    }, action = "append")
  }
  options(hf_scc.publication_writers_guarded = TRUE)
}

sha256_file <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  sub(" .*", "", output)
}

validate_submitted_panels <- function(root) {
  figure4f_path <- file.path("Figures", "Figure 3",
    "Figure3_mutations_over_time_space_data.csv")
  data <- read_publication_csv(root, figure4f_path)
  expected <- data.frame(
    plot_timepoint = rep(c("Week 8", "Week 14", "Week 17"), each = 3),
    sharing_class = rep(c("Shared with SK", "Shared with later HF",
      "Shared with earlier HF"), 3),
    numerator = c(82, 195, 0, 126, 138, 79, 107, 0, 160),
    denominator = rep(c(277, 343, 267), each = 3))
  observed <- merge(expected, data,
    by = c("plot_timepoint", "sharing_class"), sort = FALSE)
  keys <- paste(expected$plot_timepoint, expected$sharing_class)
  observed <- observed[match(keys,
    paste(observed$plot_timepoint, observed$sharing_class)), ]
  stopifnot(nrow(data) == 9L)
  stopifnot(identical(observed$pct_mutations,
    expected$numerator / expected$denominator))
  csv_hash <- sha256_file(file.path(root, figure4f_path))
  stopifnot(identical(csv_hash,
    "7460bdc3eb39e50ec8059dd0337d3e49db1c0783d1a74b532a4598cf94be25a2"))
  observed$csv_sha256 <- csv_hash

  figure5c_path <- file.path("Figures", "Figure 4",
    "Figure4_skin_gene_VAF_data.csv")
  figure5c <- read_publication_csv(root, figure5c_path)
  stopifnot(nrow(figure5c) == 30L)
  stopifnot(all(figure5c$n_mutations == figure5c$n_samples_mutated))
  stopifnot(max(figure5c$n_samples_mutated) == 2L)
  figure5c_hash <- sha256_file(file.path(root, figure5c_path))
  stopifnot(identical(figure5c_hash,
    "1a6c22ef5197b31baf126c07e504c9f8e8243a106352990e0588d05bbbaada42"))
  proof <- data.frame(panel = "Figure 5C", rows = nrow(figure5c),
    one_mutation_per_gene_sample = TRUE, maximum_mutated_samples = 2L,
    before_after_csv_sha256 = figure5c_hash,
    active_statistic = "median of per-sample mean VAF",
    compatibility_column_name = "mean_vaf", values_exactly_unchanged = TRUE)
  list(figure4f = observed, figure5c = proof)
}

reproduce_figure5g <- function(root) {
  dir <- file.path("Figures", "Figure 4")
  incidence <- read.delim(file.path(root, dir, "cBioPortal_Mutated_Genes.txt"),
    check.names = FALSE, stringsAsFactors = FALSE)
  incidence$incidence_pct <- as.numeric(sub("%", "", incidence$Freq,
    fixed = TRUE))
  incidence <- incidence[!is.na(incidence$Gene) &
    !is.na(incidence$incidence_pct), ]
  incidence <- incidence[!duplicated(incidence$Gene), ]
  candidates <- read_publication_csv(root,
    file.path(dir, "Figure4_candidate_genes_data.csv"))
  candidates <- unique(candidates[c("gene_symbol", "human_symbols")])
  values <- lapply(strsplit(candidates$human_symbols, ";", fixed = TRUE),
    function(symbols) incidence$incidence_pct[match(trimws(symbols), incidence$Gene)])
  candidate_incidence <- vapply(values, mean, numeric(1), na.rm = TRUE)
  stopifnot(length(candidate_incidence) == 30L, !anyNA(candidate_incidence))
  set.seed(20260516L)
  null <- vapply(seq_len(10000L), function(i) median(sample(
    incidence$incidence_pct, length(candidate_incidence), replace = FALSE)),
    numeric(1))
  tracked <- read_publication_csv(root, file.path(dir,
    "Figure4_cbioportal_incidence_permutation_long_data.csv"))
  stopifnot(identical(null, tracked$median_incidence_pct))
  observed <- median(candidate_incidence)
  empirical <- (sum(null >= observed) + 1) / (length(null) + 1)
  assert_identical(observed, 13.95, "Figure 5G observed median")
  assert_identical(empirical, 9.999000099990002e-5,
    "Figure 5G empirical P")
  data.frame(panel = "Figure 5G", statistic = "candidate human incidence permutation",
    summary = "median", n_candidates = length(candidate_incidence),
    universe_genes = nrow(incidence), n_iterations = length(null),
    seed = 20260516L, observed_median_incidence_pct = observed,
    null_mean = mean(null), null_median = median(null),
    null_q2.5 = unname(quantile(null, 0.025)),
    null_q97.5 = unname(quantile(null, 0.975)),
    empirical_p_upper = empirical, live_query = FALSE,
    source = file.path(dir, "cBioPortal_Mutated_Genes.txt"))
}

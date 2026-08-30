validate_frozen_statistics <- function(root) {
  sharing <- read_publication_csv(root,
    "publication_inputs/shared_mutation_statistics.csv")
  calculated <- 100 * sharing$numerator / sharing$denominator
  stopifnot(identical(sharing$numerator, c(214L, 115L)))
  stopifnot(identical(sharing$denominator, c(8133L, 953L)))
  stopifnot(identical(calculated, sharing$exact_percentage))
  stopifnot(identical(round(calculated, c(2L, 1L)),
    sharing$rounded_percentage))

  supp3c <- read_publication_csv(root,
    "publication_inputs/suppfigure3c_beta_binomial.csv")
  target <- supp3c[supp3c$contrast == "TES Acetone / TES DT", ]
  stopifnot(nrow(target) == 1L)
  assert_identical(target$odds_ratio[[1]], 1.8121921461362716,
    "Supplementary Figure 3C odds ratio")
  assert_identical(target$p.value[[1]], 0.00063761766723724066,
    "Supplementary Figure 3C BH P")
  stopifnot(identical(target$model_rows[[1]], 13686L))
  stopifnot(identical(target$model_samples[[1]], 45L))
  list(sharing = sharing, supp3c = supp3c)
}

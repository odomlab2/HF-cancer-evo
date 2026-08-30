reproduce_figure2a <- function(root) {
  files <- c(TES = "Figure2_TES_muts_mbp_data.csv",
    WES = "Figure2_WES_muts_mbp_data.csv")
  result <- do.call(rbind, lapply(names(files), function(assay) {
    path <- file.path("Figures", "Figure 2", files[[assay]])
    data <- read_publication_csv(root, path)
    data <- data[data$treatment == "DT", , drop = FALSE]
    test <- kruskal.test(muts_per_mbp ~ tissue, data = data)
    counts <- table(factor(data$tissue,
      levels = c("Hair follicle", "Skin")))
    data.frame(panel = "Figure 2A", assay = assay,
      test = "Kruskal-Wallis", hf_n = unname(counts[[1]]),
      skin_n = unname(counts[[2]]), statistic = unname(test$statistic),
      df = unname(test$parameter), p_value = test$p.value,
      adjustment = "none stated", source = path)
  }))
  stopifnot(identical(result$hf_n, c(87L, 23L)))
  stopifnot(identical(result$skin_n, c(28L, 22L)))
  assert_identical(result$p_value[[1]], 1.1097650834976858e-12,
    "Figure 2A TES P")
  assert_identical(result$p_value[[2]], 0.00077689382336869459,
    "Figure 2A WES P")
  result
}

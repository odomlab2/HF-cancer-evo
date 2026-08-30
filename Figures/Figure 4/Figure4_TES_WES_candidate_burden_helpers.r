run_candidate_burden_spearman <- function(data, y, analysis_id) {
  result <- suppressWarnings(cor.test(
    data$n_mutated_candidate_genes, data[[y]], method = "spearman",
    alternative = "two.sided", exact = FALSE))
  tibble(
    analysis_id = analysis_id,
    population = "All non-outlier TES and last-cycle WES samples",
    x_variable = paste(
      "Number of distinct genes among all 30 final candidates with at least",
      "one configured nonsynonymous mutation"),
    y_variable = y,
    test_method = "Two-sided Spearman rank correlation",
    n_samples = nrow(data),
    n_TES = sum(data$dataset == "TES"),
    n_WES = sum(data$dataset == "WES"),
    n_candidate_genes_considered = first(data$n_candidate_genes_considered),
    n_tes_targeted_candidate_genes =
      first(data$n_tes_targeted_candidate_genes),
    n_zero_candidate_genes = sum(data$n_mutated_candidate_genes == 0L),
    n_zero_y = sum(data[[y]] == 0),
    spearman_rho = unname(result$estimate),
    p_value = result$p.value,
    multiple_testing_adjustment =
      "None: two prespecified descriptive correlations",
    candidate_profiling_caveat =
      first(data$candidate_profiling_caveat))
}

plot_candidate_burden_spearman <- function(data, y, y_label, test_data) {
  p_text <- ifelse(
    test_data$p_value < 0.001,
    "P <0.001",
    paste0("P = ", formatC(test_data$p_value, digits = 3, format = "f")))
  ggplot(data, aes(x = n_mutated_candidate_genes, y = .data[[y]])) +
    geom_point(
      position = position_jitter(
        width = 0.08, height = 0, seed = 20260729),
      shape = 21, size = 2, stroke = 0.3,
      colour = "black", fill = "#D79AAA") +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5)) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
      breaks = c(0, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000),
      expand = expansion(mult = c(0.03, 0.05))) +
    labs(
      title = sprintf(
        "Spearman rho = %.2f\n%s", test_data$spearman_rho, p_text),
      x = "# of mutated candidate genes",
      y = y_label) +
    my_theme +
    theme(plot.title = element_text(
      face = "plain", hjust = 0, size = 8, margin = margin(b = 2)))
}

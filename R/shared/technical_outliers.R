# Identify technical outliers from a per-sample mutation-count table.
identify_technical_outliers <- function(mutation_counts) {
  `%>%` <- dplyr::`%>%`
  coalesce <- dplyr::coalesce
  filter <- dplyr::filter
  group_by <- dplyr::group_by
  if_else <- dplyr::if_else
  left_join <- dplyr::left_join
  mutate <- dplyr::mutate
  summarise <- dplyr::summarise
  quantile <- stats::quantile
  IQR <- stats::IQR
  log_pseudocount <- 0.05

  mutation_counts <- mutation_counts %>%
    mutate(
      mutations = coalesce(mutations, 0),
      mutations_per_mbp = if_else(
        !is.na(callable_mbp) & callable_mbp > 0,
        mutations / callable_mbp,
        NA_real_),
      log_mutations_per_mbp = if_else(
        !is.na(callable_mbp) & callable_mbp > 0,
        log((mutations + log_pseudocount) / callable_mbp),
        NA_real_))

  hypermutator_thresholds <- mutation_counts %>%
    filter(!is.na(callable_mbp), callable_mbp > 0) %>%
    group_by(tissue) %>%
    summarise(
      log_mutations_per_mbp_q3 = quantile(log_mutations_per_mbp, 0.75),
      log_mutations_per_mbp_iqr = IQR(log_mutations_per_mbp),
      log_mutations_per_mbp_1_5iqr_cutoff =
        log_mutations_per_mbp_q3 + (1.5 * log_mutations_per_mbp_iqr),
      .groups = "drop")

  hypermutator_samples <- mutation_counts %>%
    left_join(hypermutator_thresholds, by = "tissue") %>%
    mutate(
      outlier = log_mutations_per_mbp >= log_mutations_per_mbp_1_5iqr_cutoff,
      outlier = replace(
        outlier,
        is.na(callable_mbp) | callable_mbp <= 0,
        TRUE))

  technical_outliers <- hypermutator_samples %>% filter(outlier == TRUE)

  list(
    mutation_counts = mutation_counts,
    hypermutator_thresholds = hypermutator_thresholds,
    hypermutator_samples = hypermutator_samples,
    technical_outliers = technical_outliers)
}

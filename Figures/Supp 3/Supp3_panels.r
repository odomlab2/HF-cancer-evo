# 28/04/2026
# Author: Yoav Avi-Guy
# Purpose: Supplementary Figure 3 panels - Deep characterization of shared
# mutations

#-------------------------------------------------------------------------------
# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(UpSetR)
library(emmeans)
library(glmmTMB)
repo_root <- normalizePath(Sys.getenv("HF_SCC_ROOT", unset = "."), mustWork = TRUE)
output_root <- Sys.getenv("HF_SCC_OUTPUT_ROOT")
if (!nzchar(output_root)) stop("Set HF_SCC_OUTPUT_ROOT to an external output directory.")

# Colours
source(file.path(repo_root, "TES", "01_source", "01_plotting.r"))

# Loading files
tes_metadata <- readRDS(file.path(repo_root, "TES/02_data/01_metadata/sample_metadata.rds"))
wes_metadata <- readRDS(file.path(repo_root, "WES/02_data/01_metadata/sample_metadata.rds"))
tes_unique <- readRDS(file.path(repo_root, "TES/02_data/02_processed/mutations_unique_dp.rds"))
wes_unique <- readRDS(file.path(repo_root, "WES/02_data/02_processed/mutations_unique_dp.rds"))
wes_shared_hf_skin <- readRDS(file.path(repo_root, "WES/02_data/02_processed/shared_muts_hf_skin.rds"))
tes_shared_hf_skin <- readRDS(file.path(repo_root, "TES/02_data/02_processed/shared_muts_hf_skin.rds"))
tes_outliers <- readRDS(file.path(repo_root, "TES/02_data/02_processed/technical_outliers.rds"))
wes_outliers <- readRDS(file.path(repo_root, "WES/02_data/02_processed/technical_outliers.rds"))
refcds_path <- Sys.getenv("HF_SCC_REFCDS")
if (!nzchar(refcds_path)) stop("Set HF_SCC_REFCDS to the reference CDS RDA file.")
load(refcds_path)

supp3_dir <- file.path(output_root, "Figures", "Supp 3")
dir.create(supp3_dir, recursive = TRUE, showWarnings = FALSE)

tes_metadata <- tes_metadata %>% filter(!str_detect(sample_name, "ds"))

min_vaf <- 0.01
variant_levels <- c("SNV", "deletion", "insertion", "substitution")
mutation_nature_levels <- c(
  "Stop gain/loss",
  "Missense",
  "Insertion",
  "Deletion",
  "DBS",
  "Synonymous",
  "Gene-flanking")
shared_mutation_nature_levels <- c(
  "Synonymous",
  "Gene-flanking",
  "Insertion",
  "Deletion",
  "DBS",
  "Missense",
  "Stop gain/loss")
shared_mutation_group_levels <- c(
  "TES\nShared with skin",
  "TES\nShared with hair follicles\nfrom different weeks",
  "WES\nShared with skin")
ttoa_distribution_group_levels <- c(
  "TES\nShared HF-SK",
  "TES\nShared HF-HF\nat different timepoints",
  "TES\nNon-shared",
  "WES\nShared HF-SK",
  "WES\nNon-shared")
ttoa_proportion_group_levels <- c(
  "WES\nAdjacent mutations",
  "WES\nDistal mutations",
  "TES\nWeek 8 shared mutations",
  "TES\nWeek 14 shared mutations",
  "TES\nWeek 17 shared mutations")
mutation_nature_colours <- c(
  "Stop gain/loss" = "#9E2F2F",
  Missense = "#3D1B1B",
  Insertion = "#7EA3CF",
  Deletion = "#436503",
  DBS = "#2F9A8B",
  Synonymous = "#B2585E",
  `Gene-flanking` = "#C7A439")
mutation_effect_levels <- c(
  "Nonsynonymous",
  "Synonymous",
  "Gene-flanking")
gene_family_levels <- c(
  "Keratin-associated proteins",
  "Genome integrity",
  "Tumour suppressor genes",
  "Epigenetic regulators",
  "Cell adhesion and cytoskeleton",
  "Ras",
  "Metabolism",
  "Other")
KNOWN_CSCC_GENES <- c(
  "Ajuba", "Arid2", "Asxl1", "Casp8", "Card11", "Ccnd1", "Cdkn2a", "Chuk",
  "Crebbp", "Egfr", "Ep300", "Erbb2", "Erbb3", "Erbb4", "Ezh2", "Fat1",
  "Hras", "Iqgap1", "Irf6", "Kmt2a", "Kmt2c", "Kmt2d", "Knstrn", "Kras",
  "Mtor", "Myh9", "Ncor1", "Nfe2l2", "Notch1", "Notch2", "Notch3", "Nrp1",
  "Pbrm1", "Pik3ca", "Ptch1", "Pten", "Rhbdf2", "Ros1", "Rras2", "Runx1",
  "Setd2", "Tert", "Tp53", "Tp63", "Tsc1", "Usp28", "Vegfa")

filter_last_cycle <- function(df) {
  if ("last_cycle" %in% names(df)) {
    df <- df %>% filter(last_cycle == TRUE)
  }
  df
}

recode_morphology <- function(category) {
  recode(
    as.character(category),
    "Visually normal" = "Morphologically normal")
}

prepare_mutation_base <- function(mutations, outliers, last_cycle = FALSE) {
  if (last_cycle) mutations <- filter_last_cycle(mutations)

  base <- mutations %>%
    plot_levels() %>%
    mutate(gt_AF = suppressWarnings(as.numeric(gt_AF))) %>%
    filter(
      !is.na(sample_name),
      !is.na(SYMBOL),
      !sample_name %in% outliers$sample_name,
      gt_AF >= min_vaf)

  if ("time" %in% names(base)) {
    base <- base %>%
      filter(!(tissue == "Hair follicle" & as.character(time) == "Week 19"))
  }

  base
}

prepare_mutation_base_with_synonymous <- function(
  mutations,
  outliers,
  last_cycle = FALSE,
  min_vaf_filter = min_vaf) {

  if (last_cycle) mutations <- filter_last_cycle(mutations)

  base <- mutations %>%
    plot_levels() %>%
    mutate(gt_AF = suppressWarnings(as.numeric(gt_AF))) %>%
    filter(
      !is.na(sample_name),
      !is.na(SYMBOL),
      !sample_name %in% outliers$sample_name)

  if (!is.null(min_vaf_filter)) {
    base <- base %>% filter(gt_AF >= min_vaf_filter)
  }

  base
}

save_ggplot_data <- function(
  plot,
  data,
  file_name,
  cols = NULL,
  width,
  height,
  units = "mm") {

  out <- if (is.null(cols)) data else dplyr::select(data, any_of(cols))

  readr::write_csv(out, paste0(file_name, "_data.csv"))
  ggsave(
    filename = paste0(file_name, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = units)

  invisible(out)
}

make_shared_mutation_groups <- function(
  mutations,
  tissue_labels,
  category_levels,
  group_levels,
  id_cols = c("CHROM", "POS", "REF", "ALT")) {

  shared_long <- mutations %>%
    filter(
      treatment == "DT",
      tissue %in% names(tissue_labels),
      category %in% category_levels) %>%
    mutate(
      tissue_label = recode(as.character(tissue), !!!tissue_labels),
      morphology = recode_morphology(category),
      group = paste(tissue_label, morphology),
      mutation_id = do.call(paste, c(across(all_of(id_cols)), sep = ":")),
      mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
    distinct(
      mutation_label,
      mutation_id,
      CHROM,
      POS,
      REF,
      ALT,
      mutation,
      SYMBOL,
      VARIANT_CLASS,
      Consequence,
      group)

  shared_matrix <- shared_long %>%
    distinct(mutation_label, group) %>%
    mutate(present = 1L) %>%
    pivot_wider(
      names_from = group,
      values_from = present,
      values_fill = 0L)

  shared_matrix[setdiff(group_levels, names(shared_matrix))] <- 0L

  shared_export <- shared_matrix %>%
    dplyr::select(mutation_label, all_of(group_levels)) %>%
    as.data.frame()
  rownames(shared_export) <- shared_export$mutation_label

  shared_summary <- shared_long %>%
    group_by(mutation_label) %>%
    summarise(
      mutation_id = dplyr::first(mutation_id),
      CHROM = dplyr::first(CHROM),
      POS = dplyr::first(POS),
      REF = dplyr::first(REF),
      ALT = dplyr::first(ALT),
      mutation = dplyr::first(mutation),
      SYMBOL = dplyr::first(SYMBOL),
      VARIANT_CLASS = dplyr::first(VARIANT_CLASS),
      Consequence = dplyr::first(Consequence),
      n_groups = n_distinct(group),
      groups = paste(sort(unique(group)), collapse = "; "),
      .groups = "drop") %>%
    filter(n_groups >= 2)

  list(
    matrix = shared_export %>%
      dplyr::select(all_of(group_levels)) %>%
      as.data.frame(),
    export = shared_export,
    long = shared_long,
    shared = shared_summary)
}

save_upset <- function(shared_groups, group_levels, file_stub) {
  readr::write_csv(
    shared_groups$export,
    paste0(file_stub, "_data.csv"))

  upset_plot <- UpSetR::upset(
      shared_groups$matrix,
      sets = group_levels,
      nsets = length(group_levels),
      nintersects = 2 ^ length(group_levels),
      keep.order = TRUE,
      order.by = c("degree", "freq"),
      mb.ratio = c(0.65, 0.35),
      mainbar.y.label = "# mutations across group intersections",
      sets.x.label = "# mutations in group",
      main.bar.color = "black",
      sets.bar.color = "black",
      matrix.color = "black",
      shade.color = "white",
      shade.alpha = 0,
      text.scale = c(1.5, 1.4, 1.2, 1.2, 1.4, 1.2),
      point.size = 3,
      line.size = 1)

  pdf(paste0(file_stub, ".pdf"), width = 6, height = 5)
  print(upset_plot)
  dev.off()

  if (interactive()) print(upset_plot)
  invisible(upset_plot)
}

summarise_shared_variant_classes <- function(shared_mutations) {
  shared_mutations %>%
    count(VARIANT_CLASS, name = "n") %>%
    mutate(
      VARIANT_CLASS = factor(VARIANT_CLASS, levels = variant_levels)) %>%
    tidyr::complete(
      VARIANT_CLASS = factor(variant_levels, levels = variant_levels),
      fill = list(n = 0L)) %>%
    mutate(
      total = sum(n),
      pct = if_else(total > 0, n / total, 0),
      label = if_else(n > 0, paste0(n, "\n", scales::percent(pct, 1)), ""))
}

plot_shared_variant_classes <- function(df) {
  ggplot(df, aes(x = "", y = n, fill = VARIANT_CLASS)) +
    geom_col(width = 1, colour = "white") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 2.5,
      lineheight = 0.85) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = col_palette$mutation_type, drop = FALSE) +
    labs(x = "", y = "") +
    my_theme +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      legend.title = element_blank())
}

summarise_shared_ta <- function(shared_mutations) {
  shared_mutations %>%
    mutate(ta_class = if_else(!is.na(mutation) & mutation == "T>A", "T>A", "Other")) %>%
    count(ta_class, name = "n") %>%
    tidyr::complete(
      ta_class = c("T>A", "Other"),
      fill = list(n = 0L)) %>%
    mutate(
      ta_class = factor(ta_class, levels = c("T>A", "Other")),
      total = sum(n),
      pct = if_else(total > 0, n / total, 0),
      label = if_else(n > 0, paste0(n, "\n", scales::percent(pct, 1)), "")) %>%
    arrange(ta_class)
}

plot_shared_ta <- function(df) {
  ggplot(df, aes(x = "", y = n, fill = ta_class)) +
    geom_col(width = 1, colour = "white") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 2.5,
      lineheight = 0.85) +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = c(
        "T>A" = col_palette$mutation_signature[["T>A"]],
        Other = "#D9D9D9"),
      drop = FALSE) +
    labs(x = "", y = "") +
    my_theme +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      legend.title = element_blank())
}

plot_class_with_synonymous <- function(consequence, variant_class, ref, alt) {
  consequence <- coalesce(as.character(consequence), "")
  variant_class <- coalesce(as.character(variant_class), "")
  ref <- coalesce(as.character(ref), "")
  alt <- coalesce(as.character(alt), "")

  case_when(
    str_detect(consequence, "stop_gained|stop_lost") ~ "Stop gain/loss",
    str_detect(consequence, "missense_variant") ~ "Missense",
    str_detect(consequence, "synonymous_variant") ~ "Synonymous",
    str_detect(consequence, "intron_variant|UTR_variant") ~ "Gene-flanking",
    variant_class == "insertion" | nchar(ref) < nchar(alt) ~ "Insertion",
    variant_class == "deletion" | nchar(ref) > nchar(alt) ~ "Deletion",
    variant_class %in% c("substitution", "MNV", "DNV") |
      (nchar(ref) == 2 & nchar(alt) == 2) ~ "DBS",
    TRUE ~ "Gene-flanking")
}

plot_effect_with_synonymous <- function(consequence) {
  consequence <- coalesce(as.character(consequence), "")

  case_when(
    str_detect(
      consequence,
      "intron_variant|UTR_variant|upstream_gene_variant|downstream_gene_variant|non_coding_transcript") ~
      "Gene-flanking",
    str_detect(consequence, "synonymous_variant") ~ "Synonymous",
    TRUE ~ "Nonsynonymous")
}

summarise_shared_mutation_nature <- function(shared_mutations) {
  shared_mutations %>%
    mutate(
      shared_group = factor(shared_group, levels = shared_mutation_group_levels),
      mutation_nature = factor(
        mutation_nature,
        levels = shared_mutation_nature_levels)) %>%
    filter(!is.na(shared_group), !is.na(mutation_nature)) %>%
    count(shared_group, mutation_nature, name = "n") %>%
    tidyr::complete(
      shared_group = factor(
        shared_mutation_group_levels,
        levels = shared_mutation_group_levels),
      mutation_nature = factor(
        shared_mutation_nature_levels,
        levels = shared_mutation_nature_levels),
      fill = list(n = 0L)) %>%
    group_by(shared_group) %>%
    mutate(
      total = sum(n),
      pct = if_else(total > 0, n / total, 0)) %>%
    ungroup()
}

plot_shared_mutation_nature <- function(df) {
  ggplot(df, aes(x = shared_group, y = pct, fill = mutation_nature)) +
    geom_col(width = 0.72, colour = "black") +
    scale_fill_manual(
      values = mutation_nature_colours,
      breaks = shared_mutation_nature_levels,
      drop = FALSE) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "",
      y = "Proportion of shared mutations (%)",
      fill = "") +
    my_theme +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(lineheight = 0.85))
}

summarise_mutation_nature_composition <- function(
  mutations,
  group_levels,
  group_col = "composition_group") {

  mutations %>%
    mutate(
      composition_group = factor(.data[[group_col]], levels = group_levels),
      mutation_nature = factor(
        mutation_nature,
        levels = shared_mutation_nature_levels)) %>%
    filter(!is.na(composition_group), !is.na(mutation_nature)) %>%
    count(composition_group, mutation_nature, name = "n") %>%
    tidyr::complete(
      composition_group = factor(group_levels, levels = group_levels),
      mutation_nature = factor(
        shared_mutation_nature_levels,
        levels = shared_mutation_nature_levels),
      fill = list(n = 0L)) %>%
    group_by(composition_group) %>%
    mutate(
      total = sum(n),
      pct = if_else(total > 0, n / total, 0)) %>%
    ungroup()
}

plot_mutation_nature_composition <- function(df) {
  ggplot(df, aes(x = composition_group, y = pct, fill = mutation_nature)) +
    geom_col(width = 0.72, colour = "black") +
    scale_fill_manual(
      values = mutation_nature_colours,
      breaks = shared_mutation_nature_levels,
      drop = FALSE) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "",
      y = "Proportion of shared mutations (%)",
      fill = "") +
    my_theme +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(lineheight = 0.85))
}

summarise_ta_distribution <- function(mutations) {
  mutations %>%
    mutate(
      shared_group = factor(shared_group, levels = ttoa_distribution_group_levels),
      ta_class = if_else(
        !is.na(mutation) & mutation == "T>A",
        "T>A",
        "Other")) %>%
    filter(!is.na(shared_group)) %>%
    count(shared_group, ta_class, name = "n") %>%
    tidyr::complete(
      shared_group = factor(
        ttoa_distribution_group_levels,
        levels = ttoa_distribution_group_levels),
      ta_class = c("T>A", "Other"),
      fill = list(n = 0L)) %>%
    group_by(shared_group) %>%
    mutate(
      ta_class = factor(ta_class, levels = c("T>A", "Other")),
      total = sum(n),
      pct = if_else(total > 0, n / total, 0),
      label = if_else(n > 0, paste0(n, "\n", scales::percent(pct, 1)), "")) %>%
    ungroup() %>%
    arrange(shared_group, ta_class)
}

plot_ta_doughnuts_by_group <- function(df) {
  ggplot(df, aes(x = 2, y = n, fill = ta_class)) +
    geom_col(width = 1, colour = NA) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 8 / ggplot2::.pt,
      lineheight = 0.85) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    facet_wrap(~ shared_group, nrow = 1, scales = "free") +
    scale_fill_manual(
      values = c(
        "T>A" = col_palette$mutation_signature[["T>A"]],
        Other = "#D9D9D9"),
      drop = FALSE) +
    labs(x = "", y = "", fill = "") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      strip.text.x = element_text(size = 8, lineheight = 0.85),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.position = "bottom")
}

summarise_ta_proportions <- function(mutations) {
  mutations %>%
    mutate(
      ttoa_group = factor(ttoa_group, levels = ttoa_proportion_group_levels),
      ta_class = if_else(
        !is.na(mutation) & mutation == "T>A",
        "T>A",
        "Other")) %>%
    filter(!is.na(ttoa_group)) %>%
    distinct(ttoa_group, mutation_label, ta_class) %>%
    count(ttoa_group, ta_class, name = "n") %>%
    tidyr::complete(
      ttoa_group = factor(
        ttoa_proportion_group_levels,
        levels = ttoa_proportion_group_levels),
      ta_class = c("T>A", "Other"),
      fill = list(n = 0L)) %>%
    group_by(ttoa_group) %>%
    mutate(
      ta_class = factor(ta_class, levels = c("Other", "T>A")),
      total = sum(n),
      pct = if_else(total > 0, n / total, 0)) %>%
    ungroup() %>%
    arrange(ttoa_group, ta_class)
}

plot_ta_proportions <- function(df) {
  df <- df %>%
    mutate(ttoa_group = droplevels(ttoa_group))

  ggplot(df, aes(x = ttoa_group, y = pct, fill = ta_class)) +
    geom_col(width = 0.72, colour = "black") +
    scale_fill_manual(
      values = c(
        Other = "#D9D9D9",
        "T>A" = col_palette$mutation_signature[["T>A"]]),
      breaks = c("T>A", "Other"),
      drop = FALSE) +
    scale_y_continuous(
      breaks = c(0, 1),
      labels = c("0", "100"),
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "",
      y = "Proportion of shared mutations",
      fill = "") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(size = 8, lineheight = 0.85),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.position = "bottom")
}

run_ta_proportion_tests <- function(df, tes_trend_data) {
  counts <- df %>%
    mutate(
      ttoa_group = as.character(ttoa_group),
      ta_class = as.character(ta_class)) %>%
    dplyr::select(ttoa_group, ta_class, n) %>%
    tidyr::pivot_wider(
      names_from = ta_class,
      values_from = n,
      values_fill = 0L) %>%
    mutate(
      total = `T>A` + Other,
      ttoa_proportion = if_else(total > 0, `T>A` / total, NA_real_))

  wes_counts <- counts %>%
    filter(
      ttoa_group %in% c(
        "WES\nAdjacent mutations",
        "WES\nDistal mutations")) %>%
    arrange(factor(
      ttoa_group,
      levels = c("WES\nAdjacent mutations", "WES\nDistal mutations")))

  wes_table <- as.matrix(wes_counts %>% dplyr::select(`T>A`, Other))
  rownames(wes_table) <- wes_counts$ttoa_group

  wes_fisher_greater <- fisher.test(wes_table, alternative = "greater")
  wes_fisher_two_sided <- fisher.test(wes_table, alternative = "two.sided")

  tes_trend_data <- tes_trend_data %>%
    mutate(
      log_callable = log(callable_mbp)) %>%
    filter(
      is.finite(log_callable),
      total > 0)

  tes_trend_fit <- glm(
    cbind(TA, Other) ~ week + log_callable,
    family = binomial(),
    data = tes_trend_data)
  tes_reduced_fit <- glm(
    cbind(TA, Other) ~ log_callable,
    family = binomial(),
    data = tes_trend_data)
  tes_trend_coef <- summary(tes_trend_fit)$coefficients["week", ]
  tes_trend_ci <- confint.default(tes_trend_fit)["week", ]
  tes_trend_lrt <- anova(
    tes_reduced_fit,
    tes_trend_fit,
    test = "Chisq")

  tests <- bind_rows(
    tibble(
      analysis = "WES Adjacent vs Distal",
      test = "Fisher's exact test",
      alternative = "Adjacent T>A proportion > Distal T>A proportion",
      estimate = unname(wes_fisher_greater$estimate),
      estimate_type = "odds_ratio",
      conf_low = wes_fisher_greater$conf.int[[1]],
      conf_high = wes_fisher_greater$conf.int[[2]],
      p_value = wes_fisher_greater$p.value),
    tibble(
      analysis = "WES Adjacent vs Distal",
      test = "Fisher's exact test",
      alternative = "two-sided",
      estimate = unname(wes_fisher_two_sided$estimate),
      estimate_type = "odds_ratio",
      conf_low = wes_fisher_two_sided$conf.int[[1]],
      conf_high = wes_fisher_two_sided$conf.int[[2]],
      p_value = wes_fisher_two_sided$p.value),
    tibble(
      analysis = "TES Week 8/14/17 trend",
      test = "Logistic regression",
      alternative = "continuous week effect adjusted for log(callable_mbp)",
      estimate = unname(tes_trend_coef[["Estimate"]]),
      estimate_type = "log_odds_slope_per_week",
      conf_low = unname(tes_trend_ci[[1]]),
      conf_high = unname(tes_trend_ci[[2]]),
      p_value = unname(tes_trend_coef[["Pr(>|z|)"]])),
    tibble(
      analysis = "TES Week 8/14/17 trend",
      test = "Likelihood ratio test",
      alternative = "continuous week effect adjusted for log(callable_mbp)",
      estimate = unname(tes_trend_lrt$Deviance[[2]]),
      estimate_type = "deviance_reduction",
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = unname(tes_trend_lrt$`Pr(>Chi)`[[2]]))) %>%
    mutate(
      odds_ratio = case_when(
        estimate_type == "log_odds_slope_per_week" ~ exp(estimate),
        estimate_type == "odds_ratio" ~ estimate,
        TRUE ~ NA_real_),
      odds_ratio_conf_low = case_when(
        estimate_type == "log_odds_slope_per_week" ~ exp(conf_low),
        estimate_type == "odds_ratio" ~ conf_low,
        TRUE ~ NA_real_),
      odds_ratio_conf_high = case_when(
        estimate_type == "log_odds_slope_per_week" ~ exp(conf_high),
        estimate_type == "odds_ratio" ~ conf_high,
        TRUE ~ NA_real_))

  list(
    tests = tests,
    counts = counts %>%
      mutate(
        assay = if_else(str_detect(ttoa_group, "^WES"), "WES", "TES"),
        week = if_else(
          str_detect(ttoa_group, "^TES\\nWeek "),
          as.numeric(str_match(ttoa_group, "Week ([0-9]+)")[, 2]),
          NA_real_)) %>%
      arrange(factor(ttoa_group, levels = ttoa_proportion_group_levels)),
    tes_trend_model_data = tes_trend_data)
}

fit_ta_betabinomial <- function(data, assay_label) {
  data <- data %>%
    mutate(
      assay = assay_label,
      shared_group = droplevels(shared_group),
      log_callable = log(min_callable_mbp)) %>%
    filter(is.finite(log_callable))

  fit_data <- data %>%
    group_by(shared_group, sample_name, min_callable_mbp) %>%
    summarise(
      TA = sum(TA),
      Other = sum(Other),
      log_callable = log(dplyr::first(min_callable_mbp)),
      .groups = "drop") %>%
    transmute(
      TA = as.integer(TA),
      Other = as.integer(Other),
      shared_group = as.character(shared_group),
      log_callable = as.numeric(log_callable)) %>%
    as.data.frame()
  fit_data <- readr::read_csv(
    I(readr::format_csv(fit_data)),
    show_col_types = FALSE)
  fit_data$shared_group <- factor(
    fit_data$shared_group,
    levels = levels(data$shared_group))

  model <- glmmTMB(
    cbind(TA, Other) ~ shared_group + log_callable,
    family = betabinomial(link = "logit"),
    control = glmmTMBControl(
      optCtrl = list(iter.max = 1e4, eval.max = 1e4)),
    data = fit_data)

  fixed_effects <- as.data.frame(summary(model)$coefficients$cond) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(
      assay = assay_label,
      pdHess = model$sdr$pdHess,
      convergence = model$fit$convergence,
      .before = 1)

  null_model <- glmmTMB(
    cbind(TA, Other) ~ log_callable,
    family = betabinomial(link = "logit"),
    control = glmmTMBControl(
      optCtrl = list(iter.max = 1e4, eval.max = 1e4)),
    data = fit_data)

  lrt_model <- glmmTMB(
    cbind(TA, Other) ~ shared_group + log_callable,
    family = betabinomial(link = "logit"),
    control = glmmTMBControl(
      optCtrl = list(iter.max = 1e4, eval.max = 1e4)),
    data = fit_data)

  lrt <- as.data.frame(anova(null_model, lrt_model)) %>%
    tibble::rownames_to_column("model") %>%
    as_tibble() %>%
    mutate(assay = assay_label, .before = 1)

  emm <- emmeans(
    model,
    ~ shared_group,
    at = list(log_callable = 0),
    type = "response")

  list(
    data = data,
    model = model,
    fixed_effects = fixed_effects,
    lrt = lrt,
    emmeans = as_tibble(as.data.frame(summary(emm))) %>%
      mutate(assay = assay_label, .before = 1),
    pairs = as_tibble(
      as.data.frame(summary(pairs(emm), adjust = "BH", infer = TRUE))) %>%
      mutate(assay = assay_label, .before = 1))
}

gene_family_label <- function(gene_symbol, human_symbols) {
  human_symbols <- coalesce(as.character(human_symbols), "")

  case_when(
    str_detect(gene_symbol, "^Krtap") |
      str_detect(human_symbols, "(^|;)KRTAP") ~
      "Keratin-associated proteins",
    gene_symbol %in% c("Akap9", "Bbs9", "Cdk8", "Cnot1", "Fbxo21",
                       "Pcnx2", "R3hdm1", "Smg1", "Tbc1d5", "Tex2",
                       "Ttll11") ~
      "Genome integrity",
    gene_symbol %in% c("Fhit", "Nf1", "Nfkb1", "Phlpp1", "Pten") ~
      "Tumour suppressor genes",
    gene_symbol %in% c("Ash1l", "Atrx", "Cmtr1", "Kdm5b", "Kdm6a",
                       "Kmt2c", "Kmt2d", "Setd2", "Smarca4") ~
      "Epigenetic regulators",
    gene_symbol %in% c("Anxa3", "Chl1", "Cobll1", "Col7a1", "Dnm2",
                       "Dst", "Fbn2", "Gpc5", "Hmcn1", "Krt6a",
                       "Lrp1", "Mag", "Myh9", "Ntng1", "Pard3",
                       "Rcsd1", "Ryr2", "Sfi1", "Syne1", "Syne2",
                       "Tnr", "Tspan10", "Zan") ~
      "Cell adhesion and cytoskeleton",
    gene_symbol %in% c("Hras", "Kras", "Rras2") ~
      "Ras",
    gene_symbol %in% c("Atp2a2", "Atp2b2", "Hcn1", "Slc6a5", "Wnk1") ~
      "Metabolism",
    TRUE ~ "Other")
}

make_wes_oncoplot_group <- function(tissue, morphology) {
  morphology <- if_else(
    morphology == "Morphologically normal",
    "Morphologically\nnormal",
    morphology)
  paste(tissue, morphology, sep = "\n")
}

tes_base <- prepare_mutation_base(
  mutations = tes_unique,
  outliers = tes_outliers)
wes_base <- prepare_mutation_base(
  mutations = wes_unique,
  outliers = wes_outliers,
  last_cycle = TRUE)
wes_base_with_synonymous <- prepare_mutation_base_with_synonymous(
  mutations = wes_unique,
  outliers = wes_outliers,
  last_cycle = TRUE)
wes_metadata_last <- filter_last_cycle(wes_metadata)

tes_group_levels <- c(
  "HF SCC",
  "HF Papilloma",
  "HF Morphologically normal",
  "SK SCC",
  "SK Papilloma",
  "SK Morphologically normal")
wes_group_levels <- c(
  "HF Morphologically normal",
  "HF Papilloma",
  "SK Morphologically normal",
  "SK Papilloma")

tes_shared_groups <- make_shared_mutation_groups(
  mutations = tes_base,
  tissue_labels = c("Hair follicle" = "HF", "Skin" = "SK"),
  category_levels = c("SCC", "Papilloma", "Visually normal"),
  group_levels = tes_group_levels,
  id_cols = c("CHROM", "POS", "REF", "ALT"))
wes_shared_groups <- make_shared_mutation_groups(
  mutations = wes_base,
  tissue_labels = c("Hair follicle" = "HF", "Skin" = "SK"),
  category_levels = c("Papilloma", "Visually normal"),
  group_levels = wes_group_levels,
  id_cols = c("CHROM", "POS", "REF", "ALT"))

figure3_tes_mutation_labels <- tes_base %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("SCC", "Papilloma")) %>%
  mutate(
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
  distinct(mutation_label) %>%
  pull(mutation_label)

supp3_tes_tumour_mutation_labels <- tes_shared_groups$export %>%
  as_tibble() %>%
  filter(
    rowSums(
      dplyr::select(
        .,
        any_of(c("HF SCC", "HF Papilloma", "SK SCC", "SK Papilloma")))) > 0) %>%
  pull(mutation_label)

tes_upset_missing_from_supp3 <- setdiff(
  figure3_tes_mutation_labels,
  supp3_tes_tumour_mutation_labels)
tes_upset_extra_non_normal <- setdiff(
  supp3_tes_tumour_mutation_labels,
  figure3_tes_mutation_labels)

if (length(tes_upset_missing_from_supp3) > 0 ||
    length(tes_upset_extra_non_normal) > 0) {
  stop(
    "Supp 3 TES UpSet does not match the Figure 3 mutation universe outside ",
    "morphologically normal-only mutations.")
}

timepoint_levels <- c("Week 8", "Week 14", "Week 17", "Skin")
time_space_plot_timepoints <- c("Week 8", "Week 14", "Week 17")

tes_panel3_presence <- tes_base %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "SCC")) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_"),
    mutation_nature = plot_class_with_synonymous(
      Consequence,
      VARIANT_CLASS,
      REF,
      ALT),
    mutation_effect = plot_effect_with_synonymous(Consequence)) %>%
  filter(!is.na(plot_timepoint)) %>%
  filter(!(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
  mutate(
    plot_timepoint = factor(plot_timepoint, levels = timepoint_levels),
    timepoint_rank = as.integer(plot_timepoint))

tes_skin_shared_mutation_labels <- tes_panel3_presence %>%
  filter(
    plot_timepoint == "Skin",
    category %in% c("SCC", "Papilloma")) %>%
  distinct(mutation_label) %>%
  pull(mutation_label)

tes_hf_rank_summary <- tes_panel3_presence %>%
  filter(tissue == "Hair follicle") %>%
  distinct(mutation_label, timepoint_rank) %>%
  group_by(mutation_label) %>%
  summarise(
    first_hf_rank = min(timepoint_rank, na.rm = TRUE),
    last_hf_rank = max(timepoint_rank, na.rm = TRUE),
    .groups = "drop")

tes_panel3_shared_mutations <- tes_panel3_presence %>%
  filter(
    tissue == "Hair follicle",
    plot_timepoint %in% time_space_plot_timepoints) %>%
  left_join(tes_hf_rank_summary, by = "mutation_label") %>%
  mutate(
    shared_skin = mutation_label %in% tes_skin_shared_mutation_labels,
    shared_different_hf_week =
      first_hf_rank < timepoint_rank | last_hf_rank > timepoint_rank,
    shared_group = case_when(
      shared_skin ~ "TES\nShared with skin",
      shared_different_hf_week ~
        "TES\nShared with hair follicles\nfrom different weeks",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(shared_group)) %>%
  distinct(
    shared_group,
    mutation_label,
    mutation_id,
    CHROM,
    POS,
    REF,
    ALT,
    mutation,
    SYMBOL,
    VARIANT_CLASS,
    Consequence,
    mutation_nature,
    mutation_effect,
    .keep_all = TRUE)

wes_shared_skin_mutations <- wes_base_with_synonymous %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "Visually normal")) %>%
  mutate(
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_"),
    mutation_nature = plot_class_with_synonymous(
      Consequence,
      VARIANT_CLASS,
      REF,
      ALT),
    mutation_effect = plot_effect_with_synonymous(Consequence)) %>%
  group_by(mutation_label) %>%
  filter(any(tissue == "Hair follicle") & any(tissue == "Skin")) %>%
  ungroup() %>%
  arrange(mutation_label, tissue, sample_name) %>%
  distinct(mutation_label, .keep_all = TRUE) %>%
  mutate(shared_group = "WES\nShared with skin")

shared_mutations_panel3 <- bind_rows(
  tes_panel3_shared_mutations,
  wes_shared_skin_mutations) %>%
  mutate(
    shared_group = factor(shared_group, levels = shared_mutation_group_levels),
    mutation_nature = factor(
      mutation_nature,
      levels = shared_mutation_nature_levels),
    mutation_effect = factor(
      mutation_effect,
      levels = mutation_effect_levels)) %>%
  filter(!is.na(shared_group), !is.na(mutation_nature))

shared_mutation_nature_df <- summarise_shared_mutation_nature(
  shared_mutations_panel3)

figure3_wes_hf_skin_max <- wes_shared_hf_skin %>%
  filter(
    hf_treatment == "DT",
    sk_treatment == "DT",
    hf_coverage >= 5) %>%
  group_by(hf_sample_name, distance) %>%
  slice_max(order_by = n_shared, with_ties = FALSE) %>%
  ungroup()

figure3_wes_hf_skin_plot_df <- figure3_wes_hf_skin_max %>%
  filter(!is.na(adjacency)) %>%
  group_by(hf_sample_name, adjacency) %>%
  slice_max(order_by = n_shared, with_ties = FALSE) %>%
  ungroup()

wes_adjacency_mutation_base <- wes_base_with_synonymous %>%
  mutate(
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_"),
    mutation_nature = plot_class_with_synonymous(
      Consequence,
      VARIANT_CLASS,
      REF,
      ALT)) %>%
  dplyr::select(
    sample_name,
    mutation_label,
    mutation_id,
    CHROM,
    POS,
    REF,
    ALT,
    mutation,
    SYMBOL,
    VARIANT_CLASS,
    Consequence,
    mutation_nature) %>%
  distinct()

figure3_wes_adjacency_shared_mutations <- figure3_wes_hf_skin_plot_df %>%
  dplyr::select(hf_sample_name, sk_sample_name, adjacency) %>%
  left_join(
    wes_adjacency_mutation_base %>%
      rename(hf_sample_name = sample_name),
    by = "hf_sample_name",
    relationship = "many-to-many") %>%
  inner_join(
    wes_adjacency_mutation_base %>%
      distinct(sample_name, mutation_label) %>%
      rename(sk_sample_name = sample_name),
    by = c("sk_sample_name", "mutation_label"),
    relationship = "many-to-many") %>%
  mutate(composition_group = as.character(adjacency)) %>%
  distinct(
    composition_group,
    mutation_label,
    mutation_id,
    CHROM,
    POS,
    REF,
    ALT,
    mutation,
    SYMBOL,
    VARIANT_CLASS,
    Consequence,
    mutation_nature)

figure3_wes_adjacency_mutation_nature_df <-
  summarise_mutation_nature_composition(
    figure3_wes_adjacency_shared_mutations,
    group_levels = c("Adjacent", "Distal"))

figure3_time_space_presence <- tes_base %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "SCC")) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_"),
    mutation_nature = plot_class_with_synonymous(
      Consequence,
      VARIANT_CLASS,
      REF,
      ALT)) %>%
  filter(!is.na(plot_timepoint)) %>%
  filter(!(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
  distinct(
    mutation_label,
    plot_timepoint,
    tissue,
    category,
    mutation,
    mutation_nature) %>%
  mutate(
    plot_timepoint = factor(plot_timepoint, levels = timepoint_levels),
    timepoint_rank = as.integer(plot_timepoint))

figure3_time_space_shared_skin_labels <- figure3_time_space_presence %>%
  filter(
    plot_timepoint == "Skin",
    category %in% c("SCC", "Papilloma")) %>%
  distinct(mutation_label) %>%
  pull(mutation_label)

figure3_time_space_flags <- tes_base %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "SCC")) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_"),
    mutation_nature = plot_class_with_synonymous(
      Consequence,
      VARIANT_CLASS,
      REF,
      ALT)) %>%
  filter(!is.na(plot_timepoint)) %>%
  filter(!(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
  distinct(
    sample_name,
    mutation_label,
    mutation,
    callable_mbp,
    plot_timepoint,
    mutation_nature) %>%
  mutate(
    plot_timepoint = factor(plot_timepoint, levels = timepoint_levels),
    timepoint_rank = as.integer(plot_timepoint)) %>%
  filter(plot_timepoint %in% time_space_plot_timepoints) %>%
  mutate(
    shared_skin = mutation_label %in% figure3_time_space_shared_skin_labels) %>%
  rowwise() %>%
  mutate(
    shared_later_hf = any(
      figure3_time_space_presence$mutation_label == mutation_label &
        figure3_time_space_presence$tissue == "Hair follicle" &
        figure3_time_space_presence$timepoint_rank > timepoint_rank),
    shared_earlier_hf = any(
      figure3_time_space_presence$mutation_label == mutation_label &
        figure3_time_space_presence$tissue == "Hair follicle" &
        figure3_time_space_presence$timepoint_rank < timepoint_rank)) %>%
  ungroup()

figure3_time_space_shared_mutations <- figure3_time_space_flags %>%
  filter(shared_skin | shared_later_hf | shared_earlier_hf) %>%
  mutate(composition_group = as.character(plot_timepoint)) %>%
  distinct(composition_group, mutation_label, mutation, mutation_nature)

figure3_time_space_shared_mutation_sample_counts <- figure3_time_space_flags %>%
  filter(shared_skin | shared_later_hf | shared_earlier_hf) %>%
  distinct(sample_name, mutation_label, mutation, callable_mbp, plot_timepoint) %>%
  mutate(
    week = as.numeric(str_match(as.character(plot_timepoint), "Week ([0-9]+)")[, 2]),
    TA = as.integer(!is.na(mutation) & mutation == "T>A"),
    Other = 1L - TA) %>%
  filter(!is.na(week)) %>%
  group_by(sample_name, plot_timepoint, week, callable_mbp) %>%
  summarise(
    TA = sum(TA),
    Other = sum(Other),
    total = TA + Other,
    .groups = "drop")

figure3_time_space_mutation_nature_df <-
  summarise_mutation_nature_composition(
    figure3_time_space_shared_mutations,
    group_levels = time_space_plot_timepoints)

#-------------------------------------------------------------------------------
### TES raw mutation sharing across time and skin ###
#-------------------------------------------------------------------------------

time_space_raw_count_levels <- c(
  "Shared with SK",
  "Shared with later HF",
  "Shared with earlier HF",
  "Non-shared")

tes_time_space_mutation_flags <- tes_panel3_presence %>%
  filter(
    tissue == "Hair follicle",
    plot_timepoint %in% time_space_plot_timepoints) %>%
  distinct(
    sample_name,
    mutation_label,
    mutation,
    callable_mbp,
    plot_timepoint,
    timepoint_rank) %>%
  mutate(
    shared_skin = mutation_label %in% tes_skin_shared_mutation_labels) %>%
  rowwise() %>%
  mutate(
    shared_later_hf = any(
      tes_panel3_presence$mutation_label == mutation_label &
        tes_panel3_presence$tissue == "Hair follicle" &
        tes_panel3_presence$timepoint_rank > timepoint_rank),
    shared_earlier_hf = any(
      tes_panel3_presence$mutation_label == mutation_label &
        tes_panel3_presence$tissue == "Hair follicle" &
        tes_panel3_presence$timepoint_rank < timepoint_rank)) %>%
  ungroup() %>%
  mutate(
    plot_timepoint = factor(
      as.character(plot_timepoint),
      levels = time_space_plot_timepoints),
    unique_to_week = !shared_skin & !shared_later_hf & !shared_earlier_hf)

tes_ta_skin_callable <- tes_panel3_presence %>%
  filter(plot_timepoint == "Skin", category %in% c("SCC", "Papilloma")) %>%
  distinct(
    mutation_label,
    shared_sample_name = sample_name,
    shared_callable_mbp = callable_mbp)

tes_ta_hf_callable <- tes_panel3_presence %>%
  filter(tissue == "Hair follicle") %>%
  distinct(
    mutation_label,
    shared_sample_name = sample_name,
    shared_callable_mbp = callable_mbp,
    shared_timepoint_rank = timepoint_rank)

tes_ta_distribution_mutations <- tes_time_space_mutation_flags %>%
  pivot_longer(
    cols = c(
      shared_skin,
      shared_later_hf,
      shared_earlier_hf,
      unique_to_week),
    names_to = "sharing_class",
    values_to = "included") %>%
  filter(included) %>%
  mutate(
    shared_group = case_when(
      sharing_class == "shared_skin" ~ "TES\nShared HF-SK",
      sharing_class %in% c("shared_later_hf", "shared_earlier_hf") ~
        "TES\nShared HF-HF\nat different timepoints",
      sharing_class == "unique_to_week" ~ "TES\nNon-shared")) %>%
  dplyr::select(
    sample_name,
    mutation_label,
    mutation,
    callable_mbp,
    plot_timepoint,
    timepoint_rank,
    sharing_class,
    shared_group)

tes_ta_shared_skin_model_data <- tes_ta_distribution_mutations %>%
  filter(sharing_class == "shared_skin") %>%
  left_join(
    tes_ta_skin_callable,
    by = "mutation_label",
    relationship = "many-to-many") %>%
  group_by(
    sample_name,
    mutation_label,
    mutation,
    callable_mbp,
    plot_timepoint,
    timepoint_rank,
    sharing_class,
    shared_group) %>%
  summarise(
    shared_sample_names = paste(sort(unique(shared_sample_name)), collapse = "; "),
    min_callable_mbp = min(
      c(dplyr::first(callable_mbp), shared_callable_mbp),
      na.rm = TRUE),
    .groups = "drop")

tes_ta_shared_hf_model_data <- tes_ta_distribution_mutations %>%
  filter(sharing_class %in% c("shared_later_hf", "shared_earlier_hf")) %>%
  left_join(
    tes_ta_hf_callable,
    by = "mutation_label",
    relationship = "many-to-many") %>%
  filter(
    (sharing_class == "shared_later_hf" &
       shared_timepoint_rank > timepoint_rank) |
      (sharing_class == "shared_earlier_hf" &
         shared_timepoint_rank < timepoint_rank)) %>%
  group_by(
    sample_name,
    mutation_label,
    mutation,
    callable_mbp,
    plot_timepoint,
    timepoint_rank,
    sharing_class,
    shared_group) %>%
  summarise(
    shared_sample_names = paste(sort(unique(shared_sample_name)), collapse = "; "),
    min_callable_mbp = min(
      c(dplyr::first(callable_mbp), shared_callable_mbp),
      na.rm = TRUE),
    .groups = "drop")

tes_ta_nonshared_model_data <- tes_ta_distribution_mutations %>%
  filter(sharing_class == "unique_to_week") %>%
  mutate(
    shared_sample_names = NA_character_,
    min_callable_mbp = callable_mbp)

wes_ta_distribution_mutations <- wes_base_with_synonymous %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "Visually normal")) %>%
  mutate(
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_"))

wes_ta_distribution_skin_samples <- wes_ta_distribution_mutations %>%
  filter(tissue == "Skin") %>%
  distinct(
    mutation_label,
    shared_sample_name = sample_name,
    shared_callable_mbp = callable_mbp)

wes_ta_distribution_skin_labels <- wes_ta_distribution_skin_samples %>%
  distinct(mutation_label) %>%
  pull(mutation_label)

wes_ta_distribution_mutations <- wes_ta_distribution_mutations %>%
  filter(tissue == "Hair follicle") %>%
  distinct(sample_name, mutation_label, mutation, callable_mbp) %>%
  mutate(
    shared_group = if_else(
      mutation_label %in% wes_ta_distribution_skin_labels,
      "WES\nShared HF-SK",
      "WES\nNon-shared")) %>%
  dplyr::select(
    sample_name,
    mutation_label,
    mutation,
    callable_mbp,
    shared_group)

wes_ta_shared_skin_model_data <- wes_ta_distribution_mutations %>%
  filter(shared_group == "WES\nShared HF-SK") %>%
  left_join(
    wes_ta_distribution_skin_samples,
    by = "mutation_label",
    relationship = "many-to-many") %>%
  group_by(sample_name, mutation_label, mutation, callable_mbp, shared_group) %>%
  summarise(
    shared_sample_names = paste(sort(unique(shared_sample_name)), collapse = "; "),
    min_callable_mbp = min(
      c(dplyr::first(callable_mbp), shared_callable_mbp),
      na.rm = TRUE),
    .groups = "drop")

wes_ta_nonshared_model_data <- wes_ta_distribution_mutations %>%
  filter(shared_group == "WES\nNon-shared") %>%
  mutate(
    shared_sample_names = NA_character_,
    min_callable_mbp = callable_mbp)

shared_ta <- bind_rows(
  tes_ta_distribution_mutations,
  wes_ta_distribution_mutations) %>%
  summarise_ta_distribution()

shared_ta_proportions <- bind_rows(
  figure3_wes_adjacency_shared_mutations %>%
    transmute(
      ttoa_group = recode(
        as.character(composition_group),
        Adjacent = "WES\nAdjacent mutations",
        Distal = "WES\nDistal mutations"),
      mutation_label,
      mutation),
  figure3_time_space_shared_mutations %>%
    transmute(
      ttoa_group = recode(
        as.character(composition_group),
        "Week 8" = "TES\nWeek 8 shared mutations",
        "Week 14" = "TES\nWeek 14 shared mutations",
        "Week 17" = "TES\nWeek 17 shared mutations"),
      mutation_label,
      mutation)) %>%
  summarise_ta_proportions()

shared_ta_proportions_wes <- shared_ta_proportions %>%
  filter(
    ttoa_group %in% c(
      "WES\nAdjacent mutations",
      "WES\nDistal mutations")) %>%
  mutate(
    ttoa_group = factor(
      ttoa_group,
      levels = c(
        "WES\nAdjacent mutations",
        "WES\nDistal mutations")))

shared_ta_proportions_tes <- shared_ta_proportions %>%
  filter(
    ttoa_group %in% c(
      "TES\nWeek 8 shared mutations",
      "TES\nWeek 14 shared mutations",
      "TES\nWeek 17 shared mutations")) %>%
  mutate(
    ttoa_group = factor(
      ttoa_group,
      levels = c(
        "TES\nWeek 8 shared mutations",
        "TES\nWeek 14 shared mutations",
        "TES\nWeek 17 shared mutations")))

shared_ta_proportion_tests <- run_ta_proportion_tests(
  shared_ta_proportions,
  figure3_time_space_shared_mutation_sample_counts)
shared_ta_test_results <- shared_ta_proportion_tests$tests
shared_ta_test_counts <- shared_ta_proportion_tests$counts
shared_ta_tes_trend_model_data <-
  shared_ta_proportion_tests$tes_trend_model_data

tes_ta_model_data <- bind_rows(
  tes_ta_shared_skin_model_data,
  tes_ta_shared_hf_model_data,
  tes_ta_nonshared_model_data) %>%
  mutate(
    assay = "TES",
    shared_group = factor(
      shared_group,
      levels = ttoa_distribution_group_levels[1:3]),
    TA = as.integer(!is.na(mutation) & mutation == "T>A"),
    Other = 1L - TA)

wes_ta_model_data <- bind_rows(
  wes_ta_shared_skin_model_data,
  wes_ta_nonshared_model_data) %>%
  mutate(
    assay = "WES",
    shared_group = factor(
      shared_group,
      levels = ttoa_distribution_group_levels[4:5]),
    TA = as.integer(!is.na(mutation) & mutation == "T>A"),
    Other = 1L - TA)

tes_ta_betabinomial <- fit_ta_betabinomial(tes_ta_model_data, "TES")
wes_ta_betabinomial <- fit_ta_betabinomial(wes_ta_model_data, "WES")

ttoa_betabinomial_model_data <- bind_rows(
  tes_ta_betabinomial$data,
  wes_ta_betabinomial$data)
ttoa_betabinomial_fixed_effects <- bind_rows(
  tes_ta_betabinomial$fixed_effects,
  wes_ta_betabinomial$fixed_effects)
ttoa_betabinomial_lrt <- bind_rows(
  tes_ta_betabinomial$lrt,
  wes_ta_betabinomial$lrt)
ttoa_betabinomial_emmeans <- bind_rows(
  tes_ta_betabinomial$emmeans,
  wes_ta_betabinomial$emmeans)
ttoa_betabinomial_pairs <- bind_rows(
  tes_ta_betabinomial$pairs,
  wes_ta_betabinomial$pairs)

tes_time_space_raw_counts <- tes_time_space_mutation_flags %>%
  pivot_longer(
    cols = c(
      shared_skin,
      shared_later_hf,
      shared_earlier_hf,
      unique_to_week),
    names_to = "sharing_class",
    values_to = "included") %>%
  filter(included) %>%
  mutate(
    sharing_class = recode(
      sharing_class,
      shared_skin = "Shared with SK",
      shared_later_hf = "Shared with later HF",
      shared_earlier_hf = "Shared with earlier HF",
      unique_to_week = "Non-shared"),
    sharing_class = factor(
      sharing_class,
      levels = time_space_raw_count_levels)) %>%
  count(plot_timepoint, sharing_class, name = "n_mutations") %>%
  tidyr::complete(
    plot_timepoint = factor(
      time_space_plot_timepoints,
      levels = time_space_plot_timepoints),
    sharing_class = factor(
      time_space_raw_count_levels,
      levels = time_space_raw_count_levels),
    fill = list(n_mutations = 0L)) %>%
  left_join(
    tes_time_space_mutation_flags %>%
      count(plot_timepoint, name = "total_week_mutations"),
    by = "plot_timepoint") %>%
  mutate(
    total_week_mutations = coalesce(total_week_mutations, 0L),
    pct_week_mutations = if_else(
      total_week_mutations > 0,
      n_mutations / total_week_mutations,
      0))

TES_mutations_over_time_space_raw_counts <- ggplot(
  tes_time_space_raw_counts,
  aes(x = plot_timepoint, y = n_mutations, fill = sharing_class)) +
  geom_col(
    width = 0.72,
    colour = "black") +
  geom_text(
    aes(
      label = if_else(n_mutations > 0, as.character(n_mutations), ""),
      group = sharing_class),
    position = position_stack(vjust = 0.5),
    size = 8 / ggplot2::.pt) +
  scale_fill_manual(
    values = c(
      "Shared with SK" = "#f0b981",
      "Shared with later HF" = "#300358",
      "Shared with earlier HF" = "#9780aa",
      "Non-shared" = "#5B6770"),
    drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "",
    y = "# mutations",
    fill = "") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = grid::unit(0.5, "cm"),
    legend.position = "bottom")

ggplot2::set_last_plot(TES_mutations_over_time_space_raw_counts)
export_plot_data(
  data = tes_time_space_raw_counts,
  file_name = file.path(supp3_dir, "Supp3_TES_mutations_over_time_space_raw_counts"),
  cols = c(
    "plot_timepoint",
    "sharing_class",
    "n_mutations",
    "total_week_mutations",
    "pct_week_mutations"),
  width = 105,
  height = 65)
readr::write_csv(
  tes_time_space_mutation_flags,
  file.path(
    supp3_dir,
    "Supp3_TES_mutations_over_time_space_raw_count_flags_data.csv"))
if (interactive()) TES_mutations_over_time_space_raw_counts

#-------------------------------------------------------------------------------
### WES per-sample shared mutation counts ###
#-------------------------------------------------------------------------------

wes_hf_skin_adjacency <- wes_shared_hf_skin %>%
  filter(
    hf_treatment == "DT",
    sk_treatment == "DT",
    same_mouse) %>%
  group_by(hf_sample_name, adjacency) %>%
  slice_max(order_by = n_shared, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(adjacency = factor(adjacency, levels = c("Distal", "Adjacent")))

wes_paired_shared <- wes_hf_skin_adjacency %>%
  dplyr::select(hf_sample_name, adjacency, n_shared) %>%
  distinct() %>%
  pivot_wider(
    names_from = adjacency,
    values_from = n_shared) %>%
  filter(!is.na(Adjacent), !is.na(Distal)) %>%
  mutate(delta = case_when(
    Adjacent > Distal ~ "Adj",
    Distal > Adjacent ~ "Dist",
    TRUE ~ "Tie"))

wes_paired_shared_long <- wes_paired_shared %>%
  pivot_longer(
    cols = c(Adjacent, Distal),
    names_to = "adjacency",
    values_to = "n_shared") %>%
  mutate(adjacency = factor(adjacency, levels = c("Distal", "Adjacent")))

WES_shared_mutations_per_sample <- ggplot(
  wes_paired_shared_long,
  aes(x = adjacency, y = n_shared, group = hf_sample_name, colour = delta)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  facet_wrap(~ delta) +
  scale_colour_manual(
    values = c(
      Adj = "#26547c",
      Dist = "#ef476f",
      Tie = "#ffd166")) +
  scale_y_continuous(
    breaks = c(0, 130),
    limits = c(0, 130),
    expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "# of shared mutations",
    colour = NULL) +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 90),
    axis.title = element_text(size = 8),
    legend.position = "none")

export_plot_data(
  data = wes_paired_shared_long,
  file_name = file.path(supp3_dir, "Supp3_WES_shared_mutations_per_sample"),
  cols = c("hf_sample_name", "adjacency", "n_shared", "delta"),
  panel_width = 9.667)
if (interactive()) WES_shared_mutations_per_sample

set.seed(20260504)
n_distal_resample <- 10000L

hf_distal_wr_input <- wes_shared_hf_skin %>%
  filter(
    hf_treatment == "DT",
    sk_treatment == "DT",
    hf_coverage >= 5,
    same_mouse,
    adjacency %in% c("Adjacent", "Distal")) %>%
  group_by(hf_sample_name) %>%
  summarise(
    n_adj = sum(adjacency == "Adjacent"),
    n_distal = sum(adjacency == "Distal"),
    observed_adj_max = if (any(adjacency == "Adjacent")) {
      max(n_shared[adjacency == "Adjacent"])
    } else {
      NA_real_
    },
    distal_shared = list(n_shared[adjacency == "Distal"]),
    .groups = "drop") %>%
  mutate(
    included_in_distal_wr = !is.na(observed_adj_max) &
      n_adj > 0 &
      n_distal > 0)

hf_distal_wr_excluded <- hf_distal_wr_input %>%
  filter(!included_in_distal_wr) %>%
  dplyr::select(hf_sample_name, n_adj, n_distal)

if (nrow(hf_distal_wr_excluded) > 0) {
  message(
    "Excluded HF samples without both observed adjacent and distal samples: ",
    paste(hf_distal_wr_excluded$hf_sample_name, collapse = ", "))
}

hf_distal_wr_test_input <- hf_distal_wr_input %>%
  filter(included_in_distal_wr)

stopifnot(nrow(hf_distal_wr_test_input) > 0)

distal_wr_max_by_hf <- vapply(seq_len(n_distal_resample), function(i) {
  vapply(seq_len(nrow(hf_distal_wr_test_input)), function(j) {
    distal_shared <- hf_distal_wr_test_input$distal_shared[[j]]
    sampled_distal <- distal_shared[
      sample.int(
        length(distal_shared),
        size = hf_distal_wr_test_input$n_adj[[j]],
        replace = TRUE)]
    max(sampled_distal)
  }, numeric(1))
}, numeric(nrow(hf_distal_wr_test_input)))

rownames(distal_wr_max_by_hf) <- hf_distal_wr_test_input$hf_sample_name

T_obs_distal_wr <- sum(hf_distal_wr_test_input$observed_adj_max)
T_null_distal_wr <- colSums(distal_wr_max_by_hf)
adjacent_minus_distal_wr <- T_obs_distal_wr - T_null_distal_wr
adjacent_gt_distal_pct <- mean(adjacent_minus_distal_wr > 0) * 100
adjacent_gt_distal_title <- paste0(
  "Adjacent shared mutations > distal shared mutations in ",
  scales::number(adjacent_gt_distal_pct, accuracy = 1),
  " % iterations")
round_axis_limit <- function(x) {
  if (!is.finite(x) || x <= 0) return(0)
  scale <- 10 ^ max(0, floor(log10(x)) - 1)
  ceiling(x / scale) * scale
}
adjacent_minus_distal_xmax <- round_axis_limit(max(adjacent_minus_distal_wr))

global_distal_wr_summary <- data.frame(
  n_hf_samples = nrow(hf_distal_wr_test_input),
  n_resamples = n_distal_resample,
  T_obs = T_obs_distal_wr,
  T_null_median = stats::median(T_null_distal_wr),
  T_null_q025 = as.numeric(stats::quantile(T_null_distal_wr, 0.025)),
  T_null_q975 = as.numeric(stats::quantile(T_null_distal_wr, 0.975)),
  effect_vs_null_median = T_obs_distal_wr - stats::median(T_null_distal_wr),
  p_one_sided = (sum(T_null_distal_wr >= T_obs_distal_wr) + 1) /
    (n_distal_resample + 1))

distal_wr_null_global <- data.frame(
  T_null = T_null_distal_wr,
  T_obs = T_obs_distal_wr,
  adjacent_minus_distal = adjacent_minus_distal_wr,
  adjacent_gt_distal_pct = adjacent_gt_distal_pct,
  p_one_sided = global_distal_wr_summary$p_one_sided)
adjacent_minus_distal_ymax <- ggplot_build(
  ggplot(distal_wr_null_global, aes(x = adjacent_minus_distal)) +
    geom_histogram(bins = 40))$data[[1]]$count %>%
  max() %>%
  round_axis_limit()

WES_adjacent_replacement_distal_null_global <- ggplot(
  distal_wr_null_global,
  aes(x = adjacent_minus_distal)) +
  geom_histogram(bins = 40, color = "white", fill = "grey50") +
  annotate(
    "text",
    x = adjacent_minus_distal_xmax,
    y = adjacent_minus_distal_ymax,
    label = "P = 1e-04",
    hjust = 1,
    vjust = 1,
    size = 8 / ggplot2::.pt,
    color = "black") +
  scale_x_continuous(
    breaks = c(0, adjacent_minus_distal_xmax),
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    breaks = c(0, adjacent_minus_distal_ymax),
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(
    xlim = c(0, adjacent_minus_distal_xmax),
    ylim = c(0, adjacent_minus_distal_ymax)) +
  labs(
    title = adjacent_gt_distal_title,
    x = "Adjacent shared mutations - distal shared mutations, equal opportunities",
    y = "# of iterations") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5))

ggplot2::set_last_plot(WES_adjacent_replacement_distal_null_global)
export_plot_data(
  data = distal_wr_null_global,
  file_name = file.path(
    supp3_dir,
    "Supp3_WES_adjacent_max_replacement_distal_null_global"),
  cols = c(
    "T_null",
    "T_obs",
    "adjacent_minus_distal",
    "adjacent_gt_distal_pct",
    "p_one_sided"))
if (interactive()) WES_adjacent_replacement_distal_null_global

T_obs_distal_wr_median <- stats::median(hf_distal_wr_test_input$observed_adj_max)
T_null_distal_wr_median <- apply(distal_wr_max_by_hf, 2, stats::median)
adjacent_minus_distal_wr_median <- T_obs_distal_wr_median -
  T_null_distal_wr_median
adjacent_median_gt_distal_pct <- mean(adjacent_minus_distal_wr_median > 0) *
  100
adjacent_median_gt_distal_title <- paste0(
  scales::number(adjacent_median_gt_distal_pct, accuracy = 1),
  " % iterations: adjacent median > distal median")
round_axis_lower_limit <- function(x) {
  if (!is.finite(x) || x >= 0) return(0)
  -round_axis_limit(abs(x))
}
adjacent_minus_distal_median_xmin <- round_axis_lower_limit(
  min(adjacent_minus_distal_wr_median))
adjacent_minus_distal_median_xmax <- round_axis_limit(
  max(adjacent_minus_distal_wr_median))
adjacent_minus_distal_median_xbreaks <- unique(c(
  adjacent_minus_distal_median_xmin,
  0,
  adjacent_minus_distal_median_xmax))

global_distal_wr_median_summary <- data.frame(
  n_hf_samples = nrow(hf_distal_wr_test_input),
  n_resamples = n_distal_resample,
  T_obs_median = T_obs_distal_wr_median,
  T_null_median_of_medians = stats::median(T_null_distal_wr_median),
  T_null_median_q025 = as.numeric(
    stats::quantile(T_null_distal_wr_median, 0.025)),
  T_null_median_q975 = as.numeric(
    stats::quantile(T_null_distal_wr_median, 0.975)),
  effect_vs_null_median = T_obs_distal_wr_median -
    stats::median(T_null_distal_wr_median),
  p_one_sided = (sum(T_null_distal_wr_median >= T_obs_distal_wr_median) + 1) /
    (n_distal_resample + 1))

distal_wr_median_null_global <- data.frame(
  T_null_median = T_null_distal_wr_median,
  T_obs_median = T_obs_distal_wr_median,
  adjacent_minus_distal_median = adjacent_minus_distal_wr_median,
  adjacent_median_gt_distal_pct = adjacent_median_gt_distal_pct,
  p_one_sided = global_distal_wr_median_summary$p_one_sided)
p_one_sided_median_label <- paste0(
  "P = ",
  scales::number(global_distal_wr_median_summary$p_one_sided, accuracy = 1e-4))
adjacent_minus_distal_median_ymax <- ggplot_build(
  ggplot(
    distal_wr_median_null_global,
    aes(x = adjacent_minus_distal_median)) +
    geom_histogram(bins = 40))$data[[1]]$count %>%
  max() %>%
  round_axis_limit()

WES_adjacent_median_replacement_distal_null_global <- ggplot(
  distal_wr_median_null_global,
  aes(x = adjacent_minus_distal_median)) +
  geom_histogram(bins = 40, color = "white", fill = "grey50") +
  annotate(
    "text",
    x = adjacent_minus_distal_median_xmax,
    y = adjacent_minus_distal_median_ymax,
    label = p_one_sided_median_label,
    hjust = 1,
    vjust = 1,
    size = 8 / ggplot2::.pt,
    color = "black") +
  scale_x_continuous(
    breaks = adjacent_minus_distal_median_xbreaks,
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    breaks = c(0, adjacent_minus_distal_median_ymax),
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(
    xlim = c(
      adjacent_minus_distal_median_xmin,
      adjacent_minus_distal_median_xmax),
    ylim = c(0, adjacent_minus_distal_median_ymax)) +
  labs(
    title = adjacent_median_gt_distal_title,
    x = "Adjacent shared mutations - distal shared mutations,\nmedian",
    y = "# of iterations") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5))

ggplot2::set_last_plot(WES_adjacent_median_replacement_distal_null_global)
export_plot_data(
  data = distal_wr_median_null_global,
  file_name = file.path(
    supp3_dir,
    "Supp3_WES_adjacent_median_max_replacement_distal_null_global"),
  cols = c(
    "T_null_median",
    "T_obs_median",
    "adjacent_minus_distal_median",
    "adjacent_median_gt_distal_pct",
    "p_one_sided"))
if (interactive()) WES_adjacent_median_replacement_distal_null_global

wes_distal_median_distribution_xmax <- round_axis_limit(max(c(
  distal_wr_median_null_global$T_null_median,
  T_obs_distal_wr_median)))
wes_distal_median_distribution_ymax <- ggplot_build(
  ggplot(distal_wr_median_null_global, aes(x = T_null_median)) +
    geom_histogram(bins = 40))$data[[1]]$count %>%
  max() %>%
  round_axis_limit()

WES_distal_median_replacement_null_adjacent_median <- ggplot(
  distal_wr_median_null_global,
  aes(x = T_null_median)) +
  geom_histogram(bins = 40, color = "white", fill = "grey50") +
  geom_vline(
    xintercept = T_obs_distal_wr_median,
    color = "black",
    linewidth = 0.4) +
  scale_x_continuous(
    breaks = c(0, T_obs_distal_wr_median, wes_distal_median_distribution_xmax),
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    breaks = c(0, wes_distal_median_distribution_ymax),
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(
    xlim = c(0, wes_distal_median_distribution_xmax),
    ylim = c(0, wes_distal_median_distribution_ymax)) +
  labs(
    x = "Distal shared mutations, median",
    y = "# of iterations") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8))

ggplot2::set_last_plot(WES_distal_median_replacement_null_adjacent_median)
export_plot_data(
  data = distal_wr_median_null_global,
  file_name = file.path(
    supp3_dir,
    "Supp3_WES_distal_median_max_replacement_null_adjacent_median"),
  cols = c(
    "T_null_median",
    "T_obs_median",
    "adjacent_minus_distal_median",
    "adjacent_median_gt_distal_pct",
    "p_one_sided"))
if (interactive()) WES_distal_median_replacement_null_adjacent_median

set.seed(20260504)

tes_hf_distal_wr_input <- tes_shared_hf_skin %>%
  filter(
    hf_treatment == "DT",
    sk_treatment == "DT",
    hf_coverage >= 15,
    hf_mouse == sk_mouse,
    adjacency %in% c("Adjacent", "Distal")) %>%
  group_by(hf_sample_name) %>%
  summarise(
    n_adj = sum(adjacency == "Adjacent"),
    n_distal = sum(adjacency == "Distal"),
    observed_adj_max = if (any(adjacency == "Adjacent")) {
      max(n_shared[adjacency == "Adjacent"])
    } else {
      NA_real_
    },
    distal_shared = list(n_shared[adjacency == "Distal"]),
    .groups = "drop") %>%
  mutate(
    included_in_distal_wr = !is.na(observed_adj_max) &
      n_adj > 0 &
      n_distal > 0)

tes_hf_distal_wr_excluded <- tes_hf_distal_wr_input %>%
  filter(!included_in_distal_wr) %>%
  dplyr::select(hf_sample_name, n_adj, n_distal)

if (nrow(tes_hf_distal_wr_excluded) > 0) {
  message(
    "Excluded TES HF samples without both observed adjacent and distal samples: ",
    paste(tes_hf_distal_wr_excluded$hf_sample_name, collapse = ", "))
}

tes_hf_distal_wr_test_input <- tes_hf_distal_wr_input %>%
  filter(included_in_distal_wr)

stopifnot(nrow(tes_hf_distal_wr_test_input) > 0)

tes_distal_wr_max_by_hf <- vapply(seq_len(n_distal_resample), function(i) {
  vapply(seq_len(nrow(tes_hf_distal_wr_test_input)), function(j) {
    distal_shared <- tes_hf_distal_wr_test_input$distal_shared[[j]]
    sampled_distal <- distal_shared[
      sample.int(
        length(distal_shared),
        size = tes_hf_distal_wr_test_input$n_adj[[j]],
        replace = TRUE)]
    max(sampled_distal)
  }, numeric(1))
}, numeric(nrow(tes_hf_distal_wr_test_input)))

rownames(tes_distal_wr_max_by_hf) <-
  tes_hf_distal_wr_test_input$hf_sample_name

tes_T_obs_distal_wr <- sum(tes_hf_distal_wr_test_input$observed_adj_max)
tes_T_null_distal_wr <- colSums(tes_distal_wr_max_by_hf)
tes_adjacent_minus_distal_wr <- tes_T_obs_distal_wr - tes_T_null_distal_wr
tes_adjacent_gt_distal_pct <- mean(tes_adjacent_minus_distal_wr > 0) * 100
tes_adjacent_gt_distal_title <- paste0(
  "Adjacent shared mutations > distal shared mutations in ",
  scales::number(tes_adjacent_gt_distal_pct, accuracy = 1),
  " % iterations")
tes_adjacent_minus_distal_xmin <- round_axis_lower_limit(
  min(tes_adjacent_minus_distal_wr))
tes_adjacent_minus_distal_xmax <- round_axis_limit(
  max(tes_adjacent_minus_distal_wr))
tes_adjacent_minus_distal_xbreaks <- unique(c(
  tes_adjacent_minus_distal_xmin,
  0,
  tes_adjacent_minus_distal_xmax))

tes_global_distal_wr_summary <- data.frame(
  n_hf_samples = nrow(tes_hf_distal_wr_test_input),
  n_resamples = n_distal_resample,
  T_obs = tes_T_obs_distal_wr,
  T_null_median = stats::median(tes_T_null_distal_wr),
  T_null_q025 = as.numeric(stats::quantile(tes_T_null_distal_wr, 0.025)),
  T_null_q975 = as.numeric(stats::quantile(tes_T_null_distal_wr, 0.975)),
  effect_vs_null_median = tes_T_obs_distal_wr -
    stats::median(tes_T_null_distal_wr),
  p_one_sided = (sum(tes_T_null_distal_wr >= tes_T_obs_distal_wr) + 1) /
    (n_distal_resample + 1))

tes_distal_wr_null_global <- data.frame(
  T_null = tes_T_null_distal_wr,
  T_obs = tes_T_obs_distal_wr,
  adjacent_minus_distal = tes_adjacent_minus_distal_wr,
  adjacent_gt_distal_pct = tes_adjacent_gt_distal_pct,
  p_one_sided = tes_global_distal_wr_summary$p_one_sided)
tes_p_one_sided_label <- paste0(
  "P = ",
  scales::number(tes_global_distal_wr_summary$p_one_sided, accuracy = 1e-4))
tes_adjacent_minus_distal_ymax <- ggplot_build(
  ggplot(tes_distal_wr_null_global, aes(x = adjacent_minus_distal)) +
    geom_histogram(bins = 40))$data[[1]]$count %>%
  max() %>%
  round_axis_limit()

TES_adjacent_replacement_distal_null_global <- ggplot(
  tes_distal_wr_null_global,
  aes(x = adjacent_minus_distal)) +
  geom_histogram(bins = 40, color = "white", fill = "grey50") +
  annotate(
    "text",
    x = tes_adjacent_minus_distal_xmax,
    y = tes_adjacent_minus_distal_ymax,
    label = tes_p_one_sided_label,
    hjust = 1,
    vjust = 1,
    size = 8 / ggplot2::.pt,
    color = "black") +
  scale_x_continuous(
    breaks = tes_adjacent_minus_distal_xbreaks,
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    breaks = c(0, tes_adjacent_minus_distal_ymax),
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(
    xlim = c(tes_adjacent_minus_distal_xmin, tes_adjacent_minus_distal_xmax),
    ylim = c(0, tes_adjacent_minus_distal_ymax)) +
  labs(
    title = tes_adjacent_gt_distal_title,
    x = "Adjacent shared mutations - distal shared mutations, equal opportunities",
    y = "# of iterations") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5))

ggplot2::set_last_plot(TES_adjacent_replacement_distal_null_global)
export_plot_data(
  data = tes_distal_wr_null_global,
  file_name = file.path(
    supp3_dir,
    "Supp3_TES_adjacent_max_replacement_distal_null_global"),
  cols = c(
    "T_null",
    "T_obs",
    "adjacent_minus_distal",
    "adjacent_gt_distal_pct",
    "p_one_sided"))
if (interactive()) TES_adjacent_replacement_distal_null_global

tes_T_obs_distal_wr_median <- stats::median(
  tes_hf_distal_wr_test_input$observed_adj_max)
tes_T_null_distal_wr_median <- apply(
  tes_distal_wr_max_by_hf,
  2,
  stats::median)
tes_adjacent_minus_distal_wr_median <- tes_T_obs_distal_wr_median -
  tes_T_null_distal_wr_median
tes_adjacent_median_gt_distal_pct <- mean(
  tes_adjacent_minus_distal_wr_median > 0) * 100
tes_adjacent_median_gt_distal_title <- paste0(
  scales::number(tes_adjacent_median_gt_distal_pct, accuracy = 1),
  " % iterations: adjacent median > distal median")
tes_adjacent_minus_distal_median_xmin <- round_axis_lower_limit(
  min(tes_adjacent_minus_distal_wr_median))
tes_adjacent_minus_distal_median_xmax <- round_axis_limit(
  max(tes_adjacent_minus_distal_wr_median))
tes_adjacent_minus_distal_median_xbreaks <- unique(c(
  tes_adjacent_minus_distal_median_xmin,
  0,
  tes_adjacent_minus_distal_median_xmax))

tes_global_distal_wr_median_summary <- data.frame(
  n_hf_samples = nrow(tes_hf_distal_wr_test_input),
  n_resamples = n_distal_resample,
  T_obs_median = tes_T_obs_distal_wr_median,
  T_null_median_of_medians = stats::median(tes_T_null_distal_wr_median),
  T_null_median_q025 = as.numeric(
    stats::quantile(tes_T_null_distal_wr_median, 0.025)),
  T_null_median_q975 = as.numeric(
    stats::quantile(tes_T_null_distal_wr_median, 0.975)),
  effect_vs_null_median = tes_T_obs_distal_wr_median -
    stats::median(tes_T_null_distal_wr_median),
  p_one_sided = (
    sum(tes_T_null_distal_wr_median >= tes_T_obs_distal_wr_median) + 1) /
    (n_distal_resample + 1))

tes_distal_wr_median_null_global <- data.frame(
  T_null_median = tes_T_null_distal_wr_median,
  T_obs_median = tes_T_obs_distal_wr_median,
  adjacent_minus_distal_median = tes_adjacent_minus_distal_wr_median,
  adjacent_median_gt_distal_pct = tes_adjacent_median_gt_distal_pct,
  p_one_sided = tes_global_distal_wr_median_summary$p_one_sided)
tes_p_one_sided_median_label <- paste0(
  "P = ",
  scales::number(
    tes_global_distal_wr_median_summary$p_one_sided,
    accuracy = 1e-4))
tes_adjacent_minus_distal_median_ymax <- ggplot_build(
  ggplot(
    tes_distal_wr_median_null_global,
    aes(x = adjacent_minus_distal_median)) +
    geom_histogram(bins = 40))$data[[1]]$count %>%
  max() %>%
  round_axis_limit()

TES_adjacent_median_replacement_distal_null_global <- ggplot(
  tes_distal_wr_median_null_global,
  aes(x = adjacent_minus_distal_median)) +
  geom_histogram(bins = 40, color = "white", fill = "grey50") +
  annotate(
    "text",
    x = tes_adjacent_minus_distal_median_xmax,
    y = tes_adjacent_minus_distal_median_ymax,
    label = tes_p_one_sided_median_label,
    hjust = 1,
    vjust = 1,
    size = 8 / ggplot2::.pt,
    color = "black") +
  scale_x_continuous(
    breaks = tes_adjacent_minus_distal_median_xbreaks,
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    breaks = c(0, tes_adjacent_minus_distal_median_ymax),
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(
    xlim = c(
      tes_adjacent_minus_distal_median_xmin,
      tes_adjacent_minus_distal_median_xmax),
    ylim = c(0, tes_adjacent_minus_distal_median_ymax)) +
  labs(
    title = tes_adjacent_median_gt_distal_title,
    x = "Adjacent shared mutations - distal shared mutations,\nmedian",
    y = "# of iterations") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5))

ggplot2::set_last_plot(TES_adjacent_median_replacement_distal_null_global)
export_plot_data(
  data = tes_distal_wr_median_null_global,
  file_name = file.path(
    supp3_dir,
    "Supp3_TES_adjacent_median_max_replacement_distal_null_global"),
  cols = c(
    "T_null_median",
    "T_obs_median",
    "adjacent_minus_distal_median",
    "adjacent_median_gt_distal_pct",
    "p_one_sided"))
if (interactive()) TES_adjacent_median_replacement_distal_null_global

tes_distal_median_distribution_xmax <- round_axis_limit(max(c(
  tes_distal_wr_median_null_global$T_null_median,
  tes_T_obs_distal_wr_median)))
tes_distal_median_distribution_ymax <- ggplot_build(
  ggplot(tes_distal_wr_median_null_global, aes(x = T_null_median)) +
    geom_histogram(bins = 40))$data[[1]]$count %>%
  max() %>%
  round_axis_limit()

TES_distal_median_replacement_null_adjacent_median <- ggplot(
  tes_distal_wr_median_null_global,
  aes(x = T_null_median)) +
  geom_histogram(bins = 40, color = "white", fill = "grey50") +
  geom_vline(
    xintercept = tes_T_obs_distal_wr_median,
    color = "black",
    linewidth = 0.4) +
  scale_x_continuous(
    breaks = c(0, tes_T_obs_distal_wr_median, tes_distal_median_distribution_xmax),
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    breaks = c(0, tes_distal_median_distribution_ymax),
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(
    xlim = c(0, tes_distal_median_distribution_xmax),
    ylim = c(0, tes_distal_median_distribution_ymax)) +
  labs(
    x = "Distal shared mutations, median",
    y = "# of iterations") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8))

ggplot2::set_last_plot(TES_distal_median_replacement_null_adjacent_median)
export_plot_data(
  data = tes_distal_wr_median_null_global,
  file_name = file.path(
    supp3_dir,
    "Supp3_TES_distal_median_max_replacement_null_adjacent_median"),
  cols = c(
    "T_null_median",
    "T_obs_median",
    "adjacent_minus_distal_median",
    "adjacent_median_gt_distal_pct",
    "p_one_sided"))
if (interactive()) TES_distal_median_replacement_null_adjacent_median

#-------------------------------------------------------------------------------
### Shared mutation UpSet plots ###
#-------------------------------------------------------------------------------

TES_shared_mutations_morphologically_normal <- save_upset(
  shared_groups = tes_shared_groups,
  group_levels = tes_group_levels,
  file_stub = file.path(
    supp3_dir,
    "Supp3_TES_shared_mutations_morphologically_normal"))

WES_shared_mutations <- save_upset(
  shared_groups = wes_shared_groups,
  group_levels = wes_group_levels,
  file_stub = file.path(
    supp3_dir,
    "Supp3_WES_shared_mutations"))

#-------------------------------------------------------------------------------
### TES pairwise shared mutation heatmap ###
#-------------------------------------------------------------------------------

tes_pairwise_group_levels <- c(
  "Hair follicle SCC",
  "Hair follicle Papilloma",
  "Hair follicle Morphologically normal",
  "Hair follicle Acetone",
  "Skin SCC",
  "Skin Papilloma",
  "Skin Morphologically normal",
  "Skin Acetone")

tes_pairwise_base <- prepare_mutation_base_with_synonymous(
  mutations = tes_unique,
  outliers = tes_outliers,
  min_vaf_filter = min_vaf) %>%
  filter(
    !is.na(SYMBOL),
    !(tissue == "Hair follicle" &
        as.character(time) == "Week 19" &
        treatment != "Acetone"))

tes_pairwise_presence <- tes_pairwise_base %>%
  filter(
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("SCC", "Papilloma", "Visually normal", "Acetone")) %>%
  mutate(
    morphology = recode_morphology(category),
    group = paste(as.character(tissue), morphology),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
  filter(group %in% tes_pairwise_group_levels) %>%
  distinct(
    mutation_label,
    mutation_id,
    CHROM,
    POS,
    REF,
    ALT,
    mutation,
    SYMBOL,
    VARIANT_CLASS,
    Consequence,
    group)

tes_pairwise_shared_mutations <- tes_pairwise_presence %>%
  group_by(mutation_label) %>%
  summarise(
    mutation_id = dplyr::first(mutation_id),
    CHROM = dplyr::first(CHROM),
    POS = dplyr::first(POS),
    REF = dplyr::first(REF),
    ALT = dplyr::first(ALT),
    mutation = dplyr::first(mutation),
    SYMBOL = dplyr::first(SYMBOL),
    VARIANT_CLASS = dplyr::first(VARIANT_CLASS),
    Consequence = dplyr::first(Consequence),
    n_groups = n_distinct(group),
    groups = paste(sort(unique(group)), collapse = "; "),
    .groups = "drop") %>%
  filter(n_groups >= 2)

tes_pairwise_presence_matrix <- tes_pairwise_presence %>%
  semi_join(
    tes_pairwise_shared_mutations %>% dplyr::select(mutation_label),
    by = "mutation_label") %>%
  distinct(mutation_label, group) %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = group,
    values_from = present,
    values_fill = 0L)

tes_pairwise_presence_matrix[
  setdiff(tes_pairwise_group_levels, names(tes_pairwise_presence_matrix))
] <- 0L

tes_pairwise_matrix <- tes_pairwise_presence_matrix %>%
  dplyr::select(all_of(tes_pairwise_group_levels)) %>%
  as.matrix()

tes_pairwise_jaccard <- outer(
  tes_pairwise_group_levels,
  tes_pairwise_group_levels,
  Vectorize(function(group_a, group_b) {
    present_a <- tes_pairwise_matrix[, group_a] == 1L
    present_b <- tes_pairwise_matrix[, group_b] == 1L
    union_n <- sum(present_a | present_b)

    if (union_n == 0) {
      0
    } else {
      1 - (sum(present_a & present_b) / union_n)
    }
  }))
dimnames(tes_pairwise_jaccard) <- list(
  tes_pairwise_group_levels,
  tes_pairwise_group_levels)

tes_pairwise_cluster <- hclust(
  as.dist(tes_pairwise_jaccard),
  method = "average")
tes_pairwise_cluster_order <- tes_pairwise_cluster$labels[
  tes_pairwise_cluster$order]

tes_pairwise_counts <- as.data.frame(as.table(crossprod(tes_pairwise_matrix))) %>%
  setNames(c("group_x", "group_y", "n_shared")) %>%
  mutate(
    group_x = factor(group_x, levels = tes_pairwise_cluster_order),
    group_y = factor(group_y, levels = tes_pairwise_cluster_order),
    comparison_type = if_else(group_x == group_y, "Within group", "Between groups"))

tes_pairwise_group_totals <- tes_pairwise_counts %>%
  filter(group_x == group_y) %>%
  transmute(
    group_y = group_x,
    group_total_shared = n_shared)

tes_pairwise_counts <- tes_pairwise_counts %>%
  left_join(tes_pairwise_group_totals, by = "group_y") %>%
  mutate(
    pct_shared = if_else(
      group_total_shared > 0,
      n_shared / group_total_shared,
      NA_real_),
    pct_label = if_else(
      is.na(pct_shared),
      "NA",
      scales::percent(pct_shared, accuracy = 1)))

tes_pairwise_total_labels <- tes_pairwise_group_totals %>%
  mutate(
    group_y = factor(group_y, levels = tes_pairwise_cluster_order),
    group_total_label = paste0("n=", group_total_shared))

tes_pairwise_unique_pairs <- tes_pairwise_counts %>%
  filter(as.integer(group_x) < as.integer(group_y))

tes_pairwise_summary <- tes_pairwise_counts %>%
  summarise(
    n_mutations_shared_by_two_or_more_groups =
      nrow(tes_pairwise_shared_mutations),
    max_pairwise_shared =
      max(tes_pairwise_unique_pairs$n_shared, na.rm = TRUE),
    max_pct_shared =
      max(pct_shared[group_x != group_y], na.rm = TRUE),
    max_pairwise_groups = paste(
      paste(
        as.character(tes_pairwise_unique_pairs$group_x)[
          tes_pairwise_unique_pairs$n_shared ==
            max(tes_pairwise_unique_pairs$n_shared, na.rm = TRUE)
        ],
        as.character(tes_pairwise_unique_pairs$group_y)[
          tes_pairwise_unique_pairs$n_shared ==
            max(tes_pairwise_unique_pairs$n_shared, na.rm = TRUE)
        ],
        sep = " / "),
      collapse = "; "))

TES_pairwise_sharing_heatmap <- ggplot(
  tes_pairwise_counts,
  aes(x = group_x, y = group_y, fill = pct_shared)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(
    data = tes_pairwise_total_labels,
    aes(x = "Total", y = group_y, label = group_total_label),
    inherit.aes = FALSE,
    hjust = 0.5,
    size = 8 / ggplot2::.pt) +
  coord_fixed() +
  scale_x_discrete(
    limits = c(tes_pairwise_cluster_order, "Total"),
    labels = function(x) {
      x <- str_replace(x, "Morphologically normal", "Morph. normal")
      if_else(x == "Total", "n", str_replace_all(x, " ", "\n"))
    }) +
  scale_y_discrete(
    labels = function(x) {
      x <- str_replace(x, "Morphologically normal", "Morph. normal")
      str_replace_all(x, " ", "\n")
    }) +
  scale_fill_gradient(
    low = "#F2F2F2",
    high = "#1A3A63",
    breaks = c(0, 0.5, 1),
    labels = c("0", "50", "100"),
    limits = c(0, 1)) +
  labs(
    x = "",
    y = "",
    fill = "% of row\nshared mutations") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right")

ggplot2::set_last_plot(TES_pairwise_sharing_heatmap)
export_plot_data(
  data = tes_pairwise_counts,
  file_name = file.path(supp3_dir, "Supp3_TES_pairwise_sharing_heatmap"),
  cols = c(
    "group_x",
    "group_y",
    "n_shared",
    "group_total_shared",
    "pct_shared",
    "comparison_type")
  # width = 130,
  # height = 105,
  # panel_width = 62,
  # panel_height = 62
)
readr::write_csv(
  tes_pairwise_shared_mutations,
  file.path(supp3_dir, "Supp3_TES_pairwise_sharing_heatmap_mutations_data.csv"))
readr::write_csv(
  tes_pairwise_summary,
  file.path(supp3_dir, "Supp3_TES_pairwise_sharing_heatmap_summary_data.csv"))
readr::write_csv(
  tibble(
    group = tes_pairwise_cluster_order,
    cluster_rank = seq_along(tes_pairwise_cluster_order)),
  file.path(
    supp3_dir,
    "Supp3_TES_pairwise_sharing_heatmap_cluster_order_data.csv"))
if (interactive()) TES_pairwise_sharing_heatmap

tes_pairwise_group_totals_jaccard <- tes_pairwise_counts %>%
  filter(group_x == group_y) %>%
  transmute(
    group = as.character(group_x),
    group_total_shared = n_shared)

tes_pairwise_jaccard_counts <- tes_pairwise_counts %>%
  dplyr::select(group_x, group_y, n_shared, comparison_type) %>%
  mutate(
    group_x = as.character(group_x),
    group_y = as.character(group_y)) %>%
  left_join(
    tes_pairwise_group_totals_jaccard %>%
      rename(group_x = group, group_total_shared_x = group_total_shared),
    by = "group_x") %>%
  left_join(
    tes_pairwise_group_totals_jaccard %>%
      rename(group_y = group, group_total_shared_y = group_total_shared),
    by = "group_y") %>%
  mutate(
    union_shared = group_total_shared_x + group_total_shared_y - n_shared,
    jaccard_similarity = if_else(
      union_shared > 0,
      n_shared / union_shared,
      NA_real_),
    group_x = factor(group_x, levels = tes_pairwise_cluster_order),
    group_y = factor(group_y, levels = tes_pairwise_cluster_order))

TES_pairwise_sharing_heatmap_jaccard <- ggplot(
  tes_pairwise_jaccard_counts,
  aes(x = group_x, y = group_y, fill = jaccard_similarity)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  coord_fixed() +
  scale_x_discrete(
    limits = tes_pairwise_cluster_order,
    labels = function(x) {
      x <- str_replace(x, "Morphologically normal", "Morph. normal")
      str_replace_all(x, " ", "\\n")
    }) +
  scale_y_discrete(
    labels = function(x) {
      x <- str_replace(x, "Morphologically normal", "Morph. normal")
      str_replace_all(x, " ", "\\n")
    }) +
  scale_fill_gradient(
    low = "#F2F2F2",
    high = "#1A3A63",
    breaks = c(0, 0.5, 1),
    labels = c("0", "50", "100"),
    limits = c(0, 1),
    na.value = "grey80") +
  labs(
    x = "",
    y = "",
    fill = "Jaccard similarity (%)") +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right")

ggplot2::set_last_plot(TES_pairwise_sharing_heatmap_jaccard)
export_plot_data(
  data = tes_pairwise_jaccard_counts,
  file_name = file.path(
    supp3_dir,
    "Supp3_TES_pairwise_sharing_heatmap_jaccard"),
  cols = c(
    "group_x",
    "group_y",
    "n_shared",
    "group_total_shared_x",
    "group_total_shared_y",
    "union_shared",
    "jaccard_similarity",
    "comparison_type"))
if (interactive()) TES_pairwise_sharing_heatmap_jaccard

#-------------------------------------------------------------------------------
### WES shared mutations per sample count ###
#-------------------------------------------------------------------------------

wes_shared_mutation_sample_counts <- wes_base_with_synonymous %>%
  filter(
    tissue %in% c("Hair follicle", "Skin"),
    treatment == "DT",
    category %in% c("Papilloma", "Visually normal")) %>%
  mutate(
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
  group_by(mutation_label) %>%
  filter(any(tissue == "Hair follicle") & any(tissue == "Skin")) %>%
  ungroup() %>%
  distinct(mutation_label, sample_name) %>%
  count(mutation_label, name = "n_samples") %>%
  count(n_samples, name = "n_shared_mutations") %>%
  arrange(n_samples)

tes_shared_mutations_10plus_samples <- tes_pairwise_base %>%
  filter(
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("SCC", "Papilloma", "Visually normal", "Acetone")) %>%
  mutate(
    morphology = recode_morphology(category),
    group = paste(as.character(tissue), morphology),
    mutation_ID = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_ID, SYMBOL, sep = "_")) %>%
  filter(group %in% tes_pairwise_group_levels) %>%
  semi_join(
    tes_pairwise_shared_mutations %>% dplyr::select(mutation_label),
    by = "mutation_label") %>%
  distinct(
    mutation_label,
    sample_name,
    CHROM,
    POS,
    mutation,
    mutation_ID) %>%
  group_by(mutation_label) %>%
  summarise(
    CHROM = dplyr::first(CHROM),
    POS = dplyr::first(POS),
    mutation = dplyr::first(mutation),
    mutation_ID = dplyr::first(mutation_ID),
    shared_by = n_distinct(sample_name),
    .groups = "drop") %>%
  filter(shared_by >= 10) %>%
  arrange(desc(shared_by), CHROM, POS) %>%
  dplyr::select(CHROM, POS, mutation, mutation_ID, shared_by)

WES_shared_mutations_per_sample_count <- ggplot(
  wes_shared_mutation_sample_counts,
  aes(x = n_samples, y = n_shared_mutations)) +
  geom_col(fill = "#BDBDBD", colour = "black", width = 0.8) +
  scale_x_continuous(
    breaks = wes_shared_mutation_sample_counts$n_samples,
    expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "# samples sharing mutation",
    y = "# shared mutations") +
  my_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

export_plot_data(
  data = wes_shared_mutation_sample_counts,
  file_name = file.path(supp3_dir, "Supp3_WES_shared_mutations_per_sample_count"),
  cols = c("n_samples", "n_shared_mutations"))
if (interactive()) WES_shared_mutations_per_sample_count

readr::write_csv(
  tes_shared_mutations_10plus_samples,
  file.path(supp3_dir, "Supp3_TES_shared_mutations_10plus_samples_data.csv"))

#-------------------------------------------------------------------------------
### Shared mutation functional composition ###
#-------------------------------------------------------------------------------

Shared_mutation_nature_composition <- plot_shared_mutation_nature(
  shared_mutation_nature_df)
export_plot_data(
  data = shared_mutation_nature_df,
  file_name = file.path(supp3_dir, "Supp3_shared_mutation_nature_composition"),
  cols = c("shared_group", "mutation_nature", "n", "total", "pct"))
if (interactive()) Shared_mutation_nature_composition

Figure3_WES_adjacency_mutation_nature_composition <-
  plot_mutation_nature_composition(figure3_wes_adjacency_mutation_nature_df)
ggplot2::set_last_plot(Figure3_WES_adjacency_mutation_nature_composition)
export_plot_data(
  data = figure3_wes_adjacency_mutation_nature_df,
  file_name = file.path(
    supp3_dir,
    "Supp3_Figure3_wes_adjacency_mutation_nature_composition"),
  cols = c("composition_group", "mutation_nature", "n", "total", "pct"))
if (interactive()) Figure3_WES_adjacency_mutation_nature_composition

Figure3_mutations_over_time_space_mutation_nature_composition <-
  plot_mutation_nature_composition(figure3_time_space_mutation_nature_df)
ggplot2::set_last_plot(
  Figure3_mutations_over_time_space_mutation_nature_composition)
export_plot_data(
  data = figure3_time_space_mutation_nature_df,
  file_name = file.path(
    supp3_dir,
    "Supp3_Figure3_mutations_over_time_space_weighted_mutation_nature_composition"),
  cols = c("composition_group", "mutation_nature", "n", "total", "pct"))
if (interactive()) Figure3_mutations_over_time_space_mutation_nature_composition

#-------------------------------------------------------------------------------
### WES shared mutation oncoplot ###
#-------------------------------------------------------------------------------

wes_oncoplot_group_levels <- c(
  make_wes_oncoplot_group("Hair follicle", "Morphologically normal"),
  make_wes_oncoplot_group("Skin", "Morphologically normal"),
  make_wes_oncoplot_group("Hair follicle", "Papilloma"),
  make_wes_oncoplot_group("Skin", "Papilloma"))
plot_text_size <- 8
plot_text_size_mm <- plot_text_size / ggplot2::.pt

wes_shared_mutation_ids <- wes_base_with_synonymous %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "Visually normal")) %>%
  mutate(
    morphology = recode_morphology(category),
    group = paste(as.character(tissue), morphology),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
  distinct(mutation_id, group) %>%
  count(mutation_id, name = "n_groups") %>%
  filter(n_groups >= 2)

wes_oncoplot_base <- wes_base_with_synonymous %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "Visually normal")) %>%
  mutate(
    morphology = recode_morphology(category),
    biopsy_group = make_wes_oncoplot_group(tissue, morphology),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_nature = plot_class_with_synonymous(
      Consequence,
      VARIANT_CLASS,
      REF,
      ALT),
    mutation_effect = plot_effect_with_synonymous(Consequence)) %>%
  semi_join(wes_shared_mutation_ids, by = "mutation_id") %>%
  distinct(
    sample_name,
    biopsy_group,
    mutation_id,
    SYMBOL,
    mutation_nature,
    mutation_effect,
    .keep_all = TRUE)

wes_biopsy_order <- wes_metadata_last %>%
  plot_levels() %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "Visually normal"),
    !sample_name %in% wes_outliers$sample_name) %>%
  mutate(
    morphology = recode_morphology(category),
    biopsy_group = make_wes_oncoplot_group(tissue, morphology),
    biopsy_group = factor(biopsy_group, levels = wes_oncoplot_group_levels)) %>%
  semi_join(wes_oncoplot_base %>% distinct(sample_name), by = "sample_name") %>%
  arrange(biopsy_group, sample_name) %>%
  mutate(biopsy_rank = row_number()) %>%
  dplyr::select(sample_name, biopsy_group, biopsy_rank)

mouse_human_orthologues <- babelgene::orthologs(
  genes = sort(unique(na.omit(wes_oncoplot_base$SYMBOL))),
  species = "mouse",
  human = FALSE) %>%
  filter(!is.na(human_symbol), human_symbol != "") %>%
  distinct(symbol, human_symbol)

gene_family_assignments <- mouse_human_orthologues %>%
  group_by(symbol) %>%
  summarise(
    human_symbols = paste(sort(unique(human_symbol)), collapse = ";"),
    .groups = "drop") %>%
  mutate(
    gene_symbol = symbol,
    gene_family = gene_family_label(gene_symbol, human_symbols)) %>%
  dplyr::select(gene_symbol, human_symbols, gene_family)

wes_gene_summary <- wes_oncoplot_base %>%
  semi_join(wes_biopsy_order, by = "sample_name") %>%
  group_by(gene_symbol = SYMBOL) %>%
  summarise(
    n_mutated_biopsies = n_distinct(sample_name),
    n_shared_mutations = n_distinct(mutation_id),
    first_biopsy_rank = min(
      wes_biopsy_order$biopsy_rank[
        match(sample_name, wes_biopsy_order$sample_name)
      ],
      na.rm = TRUE),
    .groups = "drop") %>%
  left_join(gene_family_assignments, by = "gene_symbol") %>%
  mutate(
    gene_family = coalesce(gene_family, "Other"),
    gene_family = factor(gene_family, levels = gene_family_levels),
    known_cscc_gene = gene_symbol %in% KNOWN_CSCC_GENES) %>%
  arrange(
    desc(n_mutated_biopsies),
    desc(n_shared_mutations),
    first_biopsy_rank,
    gene_symbol) %>%
  mutate(
    gene_rank = row_number(),
    plot_y = n() - gene_rank + 1,
    gene_label = if_else(known_cscc_gene, paste0(gene_symbol, " *"), gene_symbol))

wes_plot_y_limits <- range(wes_gene_summary$plot_y) + c(-0.5, 0.5)

wes_gene_label_df <- wes_gene_summary %>%
  mutate(gene_label_strip = factor("Skin\nPapilloma"))

wes_gene_mutation_counts <- wes_oncoplot_base %>%
  filter(SYMBOL %in% wes_gene_summary$gene_symbol) %>%
  distinct(SYMBOL, mutation_id, mutation_effect) %>%
  group_by(gene_symbol = SYMBOL, mutation_id) %>%
  summarise(
    mutation_effect = case_when(
      any(mutation_effect == "Nonsynonymous") ~ "Nonsynonymous",
      any(mutation_effect == "Gene-flanking") ~ "Gene-flanking",
      TRUE ~ "Synonymous"),
    .groups = "drop") %>%
  count(gene_symbol, mutation_effect, name = "n_mutations")

wes_gene_mutation_bar_df <- expand_grid(
  gene_symbol = wes_gene_summary$gene_symbol,
  mutation_effect = mutation_effect_levels) %>%
  left_join(
    wes_gene_mutation_counts,
    by = c("gene_symbol", "mutation_effect")) %>%
  left_join(
    wes_gene_summary %>% dplyr::select(gene_symbol, plot_y),
    by = "gene_symbol") %>%
  mutate(
    n_mutations = coalesce(n_mutations, 0L),
    bar_strip = factor("Skin\nPapilloma"),
    mutation_effect = factor(mutation_effect, levels = mutation_effect_levels))

wes_gene_mutation_totals <- wes_gene_mutation_bar_df %>%
  pivot_wider(
    names_from = mutation_effect,
    values_from = n_mutations,
    values_fill = 0L) %>%
  mutate(total_mutations = Nonsynonymous + Synonymous + `Gene-flanking`) %>%
  dplyr::select(
    gene_symbol,
    Nonsynonymous,
    Synonymous,
    `Gene-flanking`,
    total_mutations)

wes_gene_mutation_bar_df <- wes_gene_mutation_bar_df %>%
  left_join(wes_gene_mutation_totals, by = "gene_symbol") %>%
  mutate(pct_mutations = if_else(total_mutations > 0, n_mutations / total_mutations, 0))

wes_gene_hits <- wes_oncoplot_base %>%
  semi_join(wes_gene_summary, by = c("SYMBOL" = "gene_symbol")) %>%
  mutate(
    mutation_nature = factor(
      mutation_nature,
      levels = mutation_nature_levels)) %>%
  group_by(gene_symbol = SYMBOL, sample_name) %>%
  summarise(
    mutation_classes = paste(sort(unique(as.character(mutation_nature))), collapse = "; "),
    mutation_ids = paste(sort(unique(mutation_id)), collapse = "; "),
    n_cell_mutations = n_distinct(mutation_id),
    mutation_nature = as.character(
      mutation_nature[which.min(as.integer(mutation_nature))]),
    .groups = "drop")

WES_shared_mutations_oncoplot_df <- expand_grid(
  gene_symbol = wes_gene_summary$gene_symbol,
  sample_name = wes_biopsy_order$sample_name) %>%
  left_join(wes_biopsy_order, by = "sample_name") %>%
  left_join(wes_gene_hits, by = c("gene_symbol", "sample_name")) %>%
  left_join(wes_gene_mutation_totals, by = "gene_symbol") %>%
  left_join(wes_gene_summary, by = "gene_symbol") %>%
  mutate(
    sample_name = factor(sample_name, levels = wes_biopsy_order$sample_name),
    mutation_nature = factor(mutation_nature, levels = mutation_nature_levels),
    mutation_nature_plot = factor(
      coalesce(as.character(mutation_nature), "Not detected"),
      levels = c(mutation_nature_levels, "Not detected")))

WES_shared_mutations_oncoplot_bar <- ggplot(
  wes_gene_mutation_bar_df,
  aes(x = pct_mutations, y = plot_y, fill = mutation_effect)) +
  geom_col(width = 0.85, orientation = "y") +
  facet_grid(. ~ bar_strip) +
  scale_x_reverse(
    limits = c(1, 0),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = wes_plot_y_limits,
    breaks = wes_gene_summary$plot_y,
    labels = NULL,
    expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(
    values = c(
      Nonsynonymous = "#496B80",
      Synonymous = "#B2585E",
      `Gene-flanking` = "#C7A439"),
    drop = FALSE) +
  labs(x = "% of mutations", y = "", fill = "") +
  my_theme +
  theme(
    text = element_text(size = plot_text_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = plot_text_size),
    legend.key.size = grid::unit(0.35, "cm"),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = plot_text_size,
      angle = 90,
      hjust = 0,
      colour = "white"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)) +
  guides(fill = guide_legend(ncol = 1))

WES_shared_mutations_oncoplot_gene_labels <- ggplot(
  wes_gene_label_df,
  aes(x = 1, y = plot_y, label = gene_symbol, colour = known_cscc_gene)) +
  geom_text(hjust = 1, size = plot_text_size_mm) +
  facet_grid(. ~ gene_label_strip) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    limits = wes_plot_y_limits,
    breaks = wes_gene_summary$plot_y,
    labels = NULL,
    expand = expansion(mult = c(0, 0))) +
  scale_colour_manual(
    values = c(
      `FALSE` = "black",
      `TRUE` = "#D63B2A"),
    guide = "none") +
  labs(x = "", y = "") +
  my_theme +
  theme(
    text = element_text(size = plot_text_size),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = plot_text_size,
      angle = 90,
      hjust = 0,
      colour = "white"),
    plot.margin = margin(5.5, 2, 5.5, 0))

WES_shared_mutations_oncoplot_heatmap <- ggplot(
  WES_shared_mutations_oncoplot_df,
  aes(x = sample_name, y = plot_y, fill = mutation_nature_plot)) +
  geom_tile(colour = "white", linewidth = 0.2, height = 0.85) +
  facet_grid(
    . ~ biopsy_group,
    scales = "free_x",
    space = "free_x") +
  scale_y_continuous(
    limits = wes_plot_y_limits,
    breaks = wes_gene_summary$plot_y,
    labels = NULL,
    expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(
    values = c(
      "Stop gain/loss" = "#9E2F2F",
      Missense = "#3D1B1B",
      Insertion = "#7EA3CF",
      Deletion = "#436503",
      DBS = "#2F9A8B",
      Synonymous = "#B2585E",
      `Gene-flanking` = "#C7A439",
      "Not detected" = "#D9D9D9"),
    limits = c(mutation_nature_levels, "Not detected"),
    breaks = mutation_nature_levels,
    drop = FALSE) +
  labs(x = "", y = "", fill = "") +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      lineheight = 0.85),
    legend.position = "bottom",
    legend.key.size = grid::unit(0.35, "cm"),
    panel.spacing.x = grid::unit(0.08, "lines")) +
  guides(fill = guide_legend(nrow = 1))

WES_shared_mutations_oncoplot <- grid::grid.grabExpr(
  {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(
        nrow = 1,
        ncol = 3,
        widths = grid::unit(c(1.4, 0.9, 6), "null"))))
    print(
      WES_shared_mutations_oncoplot_bar,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(
      WES_shared_mutations_oncoplot_gene_labels,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(
      WES_shared_mutations_oncoplot_heatmap,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  },
  wrap.grobs = TRUE)

readr::write_csv(
  WES_shared_mutations_oncoplot_df,
  file.path(supp3_dir, "Supp3_WES_shared_mutations_oncoplot_data.csv"))
readr::write_csv(
  wes_gene_summary,
  file.path(supp3_dir, "Supp3_WES_shared_mutations_gene_summary_data.csv"))

pdf(
  file.path(supp3_dir, "Supp3_WES_shared_mutations_oncoplot.pdf"),
  width = max(14, min(24, 7 + 0.16 * nrow(wes_biopsy_order))),
  height = max(5, 2 + 0.16 * nrow(wes_gene_summary)))
grid::grid.draw(WES_shared_mutations_oncoplot)
dev.off()
if (interactive()) grid::grid.draw(WES_shared_mutations_oncoplot)

save_wes_shared_oncoplot <- function(wes_base_input, file_stub) {
  shared_mutation_ids <- wes_base_input %>%
    filter(
      treatment == "DT",
      tissue %in% c("Hair follicle", "Skin"),
      category %in% c("Papilloma", "Visually normal")) %>%
    mutate(
      morphology = recode_morphology(category),
      group = paste(as.character(tissue), morphology),
      mutation_id = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
    distinct(mutation_id, group) %>%
    count(mutation_id, name = "n_groups") %>%
    filter(n_groups >= 2)

  oncoplot_base <- wes_base_input %>%
    filter(
      treatment == "DT",
      tissue %in% c("Hair follicle", "Skin"),
      category %in% c("Papilloma", "Visually normal")) %>%
    mutate(
      morphology = recode_morphology(category),
      biopsy_group = make_wes_oncoplot_group(tissue, morphology),
      mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
      mutation_nature = plot_class_with_synonymous(
        Consequence,
        VARIANT_CLASS,
        REF,
        ALT),
      mutation_effect = plot_effect_with_synonymous(Consequence)) %>%
    semi_join(shared_mutation_ids, by = "mutation_id") %>%
    distinct(
      sample_name,
      biopsy_group,
      mutation_id,
      SYMBOL,
      mutation_nature,
      mutation_effect,
      .keep_all = TRUE)

  biopsy_order <- wes_metadata_last %>%
    plot_levels() %>%
    filter(
      treatment == "DT",
      tissue %in% c("Hair follicle", "Skin"),
      category %in% c("Papilloma", "Visually normal"),
      !sample_name %in% wes_outliers$sample_name) %>%
    mutate(
      morphology = recode_morphology(category),
      biopsy_group = make_wes_oncoplot_group(tissue, morphology),
      biopsy_group = factor(biopsy_group, levels = wes_oncoplot_group_levels)) %>%
    semi_join(oncoplot_base %>% distinct(sample_name), by = "sample_name") %>%
    arrange(biopsy_group, sample_name) %>%
    mutate(biopsy_rank = row_number()) %>%
    dplyr::select(sample_name, biopsy_group, biopsy_rank)

  orthologues <- babelgene::orthologs(
    genes = sort(unique(na.omit(oncoplot_base$SYMBOL))),
    species = "mouse",
    human = FALSE) %>%
    filter(!is.na(human_symbol), human_symbol != "") %>%
    distinct(symbol, human_symbol)

  family_assignments <- orthologues %>%
    group_by(symbol) %>%
    summarise(
      human_symbols = paste(sort(unique(human_symbol)), collapse = ";"),
      .groups = "drop") %>%
    mutate(
      gene_symbol = symbol,
      gene_family = gene_family_label(gene_symbol, human_symbols)) %>%
    dplyr::select(gene_symbol, human_symbols, gene_family)

  gene_summary <- oncoplot_base %>%
    semi_join(biopsy_order, by = "sample_name") %>%
    group_by(gene_symbol = SYMBOL) %>%
    summarise(
      n_mutated_biopsies = n_distinct(sample_name),
      n_shared_mutations = n_distinct(mutation_id),
      first_biopsy_rank = min(
        biopsy_order$biopsy_rank[
          match(sample_name, biopsy_order$sample_name)
        ],
        na.rm = TRUE),
      .groups = "drop") %>%
    left_join(family_assignments, by = "gene_symbol") %>%
    mutate(
      gene_family = coalesce(gene_family, "Other"),
      gene_family = factor(gene_family, levels = gene_family_levels),
      known_cscc_gene = gene_symbol %in% KNOWN_CSCC_GENES) %>%
    group_by(gene_family) %>%
    mutate(
      family_mutated_biopsies = sum(n_mutated_biopsies),
      family_shared_mutations = sum(n_shared_mutations),
      family_n_genes = n()) %>%
    ungroup() %>%
    arrange(
      gene_family == "Other",
      desc(family_mutated_biopsies),
      desc(family_shared_mutations),
      gene_family,
      desc(n_mutated_biopsies),
      desc(n_shared_mutations),
      first_biopsy_rank,
      gene_symbol) %>%
    mutate(
      gene_rank = row_number(),
      family_rank = as.integer(factor(gene_family, levels = unique(gene_family))),
      n_gene_families = n_distinct(gene_family),
      plot_y = n() - gene_rank + 1 + (n_gene_families - family_rank) * 0.35,
      gene_label = if_else(known_cscc_gene, paste0(gene_symbol, " *"), gene_symbol))

  plot_y_limits <- range(gene_summary$plot_y) + c(-0.5, 0.5)

  gene_group_label_df <- gene_summary %>%
    group_by(gene_family) %>%
    summarise(
      plot_y = mean(range(plot_y)),
      family_mutated_biopsies = dplyr::first(family_mutated_biopsies),
      family_shared_mutations = dplyr::first(family_shared_mutations),
      .groups = "drop") %>%
    arrange(
      gene_family == "Other",
      desc(family_mutated_biopsies),
      desc(family_shared_mutations),
      gene_family)

  gene_label_df <- gene_summary %>%
    mutate(gene_label_strip = factor("Skin\nPapilloma"))

  gene_mutation_counts <- oncoplot_base %>%
    filter(SYMBOL %in% gene_summary$gene_symbol) %>%
    distinct(SYMBOL, mutation_id, mutation_effect) %>%
    group_by(gene_symbol = SYMBOL, mutation_id) %>%
    summarise(
      mutation_effect = case_when(
        any(mutation_effect == "Nonsynonymous") ~ "Nonsynonymous",
        any(mutation_effect == "Gene-flanking") ~ "Gene-flanking",
        TRUE ~ "Synonymous"),
      .groups = "drop") %>%
    count(gene_symbol, mutation_effect, name = "n_mutations")

  gene_mutation_bar_df <- expand_grid(
    gene_symbol = gene_summary$gene_symbol,
    mutation_effect = mutation_effect_levels) %>%
    left_join(gene_mutation_counts, by = c("gene_symbol", "mutation_effect")) %>%
    left_join(gene_summary %>% dplyr::select(gene_symbol, plot_y), by = "gene_symbol") %>%
    mutate(
      n_mutations = coalesce(n_mutations, 0L),
      bar_strip = factor("Skin\nPapilloma"),
      mutation_effect = factor(mutation_effect, levels = mutation_effect_levels))

  gene_mutation_totals <- gene_mutation_bar_df %>%
    pivot_wider(
      names_from = mutation_effect,
      values_from = n_mutations,
      values_fill = 0L) %>%
    mutate(total_mutations = Nonsynonymous + Synonymous + `Gene-flanking`) %>%
    dplyr::select(
      gene_symbol,
      Nonsynonymous,
      Synonymous,
      `Gene-flanking`,
      total_mutations)

  gene_mutation_bar_df <- gene_mutation_bar_df %>%
    left_join(gene_mutation_totals, by = "gene_symbol") %>%
    mutate(pct_mutations = if_else(total_mutations > 0, n_mutations / total_mutations, 0))

  gene_hits <- oncoplot_base %>%
    semi_join(gene_summary, by = c("SYMBOL" = "gene_symbol")) %>%
    mutate(
      mutation_nature = factor(mutation_nature, levels = mutation_nature_levels)) %>%
    group_by(gene_symbol = SYMBOL, sample_name) %>%
    summarise(
      mutation_classes = paste(sort(unique(as.character(mutation_nature))), collapse = "; "),
      mutation_ids = paste(sort(unique(mutation_id)), collapse = "; "),
      n_cell_mutations = n_distinct(mutation_id),
      mutation_nature = as.character(
        mutation_nature[which.min(as.integer(mutation_nature))]),
      .groups = "drop")

  oncoplot_df <- expand_grid(
    gene_symbol = gene_summary$gene_symbol,
    sample_name = biopsy_order$sample_name) %>%
    left_join(biopsy_order, by = "sample_name") %>%
    left_join(gene_hits, by = c("gene_symbol", "sample_name")) %>%
    left_join(gene_mutation_totals, by = "gene_symbol") %>%
    left_join(gene_summary, by = "gene_symbol") %>%
    mutate(
      sample_name = factor(sample_name, levels = biopsy_order$sample_name),
      mutation_nature = factor(mutation_nature, levels = mutation_nature_levels),
      mutation_nature_plot = factor(
        coalesce(as.character(mutation_nature), "Not detected"),
        levels = c(mutation_nature_levels, "Not detected")))

  oncoplot_bar <- ggplot(
    gene_mutation_bar_df,
    aes(x = pct_mutations, y = plot_y, fill = mutation_effect)) +
    geom_col(width = 0.85, orientation = "y") +
    facet_grid(. ~ bar_strip) +
    scale_x_reverse(
      limits = c(1, 0),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = plot_y_limits,
      breaks = gene_summary$plot_y,
      labels = NULL,
      expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(
      values = c(
        Nonsynonymous = "#496B80",
        Synonymous = "#B2585E",
        `Gene-flanking` = "#C7A439"),
      drop = FALSE) +
    labs(x = "% of mutations", y = "", fill = "") +
    my_theme +
    theme(
      text = element_text(size = plot_text_size),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = plot_text_size),
      legend.key.size = grid::unit(0.35, "cm"),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(
        size = plot_text_size,
        angle = 90,
        hjust = 0,
        colour = "white"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)) +
    guides(fill = guide_legend(ncol = 1))

  oncoplot_gene_labels <- ggplot(
    gene_label_df,
    aes(x = 1, y = plot_y, label = gene_symbol, colour = known_cscc_gene)) +
    geom_text(hjust = 1, size = plot_text_size_mm) +
    facet_grid(. ~ gene_label_strip) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      limits = plot_y_limits,
      breaks = gene_summary$plot_y,
      labels = NULL,
      expand = expansion(mult = c(0, 0))) +
    scale_colour_manual(
      values = c(`FALSE` = "black", `TRUE` = "#D63B2A"),
      guide = "none") +
    labs(x = "", y = "") +
    my_theme +
    theme(
      text = element_text(size = plot_text_size),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(
        size = plot_text_size,
        angle = 90,
        hjust = 0,
        colour = "white"),
      plot.margin = margin(5.5, 2, 5.5, 0))

  oncoplot_heatmap <- ggplot(
    oncoplot_df,
    aes(x = sample_name, y = plot_y, fill = mutation_nature_plot)) +
    geom_tile(colour = "white", linewidth = 0.2, height = 0.85) +
    facet_grid(. ~ biopsy_group, scales = "free_x", space = "free_x") +
    scale_y_continuous(
      limits = plot_y_limits,
      breaks = gene_summary$plot_y,
      labels = NULL,
      expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(
      values = c(
        "Stop gain/loss" = "#9E2F2F",
        Missense = "#3D1B1B",
        Insertion = "#7EA3CF",
        Deletion = "#436503",
        DBS = "#2F9A8B",
        Synonymous = "#B2585E",
        `Gene-flanking` = "#C7A439",
        "Not detected" = "#D9D9D9"),
      limits = c(mutation_nature_levels, "Not detected"),
      breaks = mutation_nature_levels,
      drop = FALSE) +
    labs(x = "", y = "", fill = "") +
    my_theme +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.85),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.35, "cm"),
      panel.spacing.x = grid::unit(0.08, "lines")) +
    guides(fill = guide_legend(nrow = 1))

  oncoplot_gene_family <- ggplot(gene_summary, aes(x = 0, y = plot_y)) +
    geom_text(
      data = gene_group_label_df,
      aes(x = 0, y = plot_y, label = gene_family),
      inherit.aes = FALSE,
      hjust = 0,
      size = plot_text_size_mm) +
    facet_grid(. ~ factor("Skin\nPapilloma")) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      limits = plot_y_limits,
      breaks = gene_summary$plot_y,
      labels = NULL,
      expand = expansion(mult = c(0, 0))) +
    labs(x = "", y = "") +
    my_theme +
    theme(
      text = element_text(size = plot_text_size),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(
        size = plot_text_size,
        angle = 90,
        hjust = 0,
        colour = "white"),
      plot.margin = margin(5.5, 5.5, 5.5, 2))

  oncoplot <- grid::grid.grabExpr(
    {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(
        layout = grid::grid.layout(
          nrow = 1,
          ncol = 4,
          widths = grid::unit(c(1.4, 0.9, 6, 2.2), "null"))))
      print(oncoplot_bar, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
      print(oncoplot_gene_labels, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
      print(oncoplot_heatmap, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
      print(oncoplot_gene_family, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 4))
    },
    wrap.grobs = TRUE)

  readr::write_csv(oncoplot_df, paste0(file_stub, "_data.csv"))
  readr::write_csv(gene_summary, paste0(file_stub, "_gene_summary_data.csv"))

  pdf(
    paste0(file_stub, ".pdf"),
    width = max(14, min(24, 7 + 0.16 * nrow(biopsy_order))),
    height = max(5, 2 + 0.16 * nrow(gene_summary)))
  grid::grid.draw(oncoplot)
  dev.off()

  invisible(oncoplot)
}

# The main WES oncoplot above uses last-cycle WES calls after outlier removal
# and the same gt_AF >= 0.01 threshold as the other WES panels.

#-------------------------------------------------------------------------------
### T>A mutation distribution ###
#-------------------------------------------------------------------------------

Shared_TtoA_distribution_WES <- plot_ta_proportions(shared_ta_proportions_wes)
ggplot2::set_last_plot(Shared_TtoA_distribution_WES)
export_plot_data(
  data = shared_ta_proportions_wes,
  file_name = file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_WES_adjacency"),
  cols = c("ttoa_group", "ta_class", "n", "total", "pct"))

Shared_TtoA_distribution_TES <- plot_ta_proportions(shared_ta_proportions_tes)
ggplot2::set_last_plot(Shared_TtoA_distribution_TES)
export_plot_data(
  data = shared_ta_proportions_tes,
  file_name = file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_TES_timepoints"),
  cols = c("ttoa_group", "ta_class", "n", "total", "pct"))
readr::write_csv(
  shared_ta_test_results,
  file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_tests.csv"))
readr::write_csv(
  shared_ta_test_counts,
  file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_test_counts.csv"))
readr::write_csv(
  shared_ta_tes_trend_model_data,
  file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_tes_trend_model_data.csv"))
readr::write_csv(
  ttoa_betabinomial_model_data,
  file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_betabinomial_model_data.csv"))
readr::write_csv(
  ttoa_betabinomial_fixed_effects,
  file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_betabinomial_fixed_effects.csv"))
readr::write_csv(
  ttoa_betabinomial_lrt,
  file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_betabinomial_lrt.csv"))
readr::write_csv(
  ttoa_betabinomial_emmeans,
  file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_betabinomial_emmeans.csv"))
readr::write_csv(
  ttoa_betabinomial_pairs,
  file.path(
    supp3_dir,
    "Supp3_shared_TtoA_distribution_betabinomial_pairs.csv"))
if (interactive()) Shared_TtoA_distribution_WES
if (interactive()) Shared_TtoA_distribution_TES

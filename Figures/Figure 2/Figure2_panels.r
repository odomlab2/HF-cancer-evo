# 10/04/2026
# Author: Yoav Avi-Guy
# Purpose: Figure 2 panels - Longitudinal profiling reveals elevated mutational 
# burden with a pecualiar signature in hair follicles and its reduction during 
# tumor development 
#-------------------------------------------------------------------------------
# Libraries
library(dplyr)
library(tidyr)
library(ggalluvial)
library(stringr)
library(emmeans)
library(glmmTMB)
repo_root <- normalizePath(Sys.getenv("HF_SCC_ROOT", unset = "."), mustWork = TRUE)
output_dir <- Sys.getenv("HF_SCC_OUTPUT_ROOT")
if (!nzchar(output_dir)) stop("Set HF_SCC_OUTPUT_ROOT to an external output directory.")
output_dir <- file.path(output_dir, "Figures", "Figure 2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Plotting parameters
source(file.path(repo_root, "TES", "01_source", "01_plotting.r"))

# Loading files
tes_metadata <- readRDS(file.path(repo_root, "TES/02_data/01_metadata/sample_metadata.rds"))
tes_mutations_df <- readRDS(file.path(repo_root, "TES/02_data/02_processed/mutations_unique_dp.rds"))
tes_outliers <- readRDS(file.path(repo_root, "TES/02_data/02_processed/technical_outliers.rds"))
wes_metadata <- readRDS(file.path(repo_root, "WES/02_data/01_metadata/sample_metadata.rds"))
wes_mutations_df <- readRDS(file.path(repo_root, "WES/02_data/02_processed/mutations_unique_dp.rds"))
wes_outliers <- readRDS(file.path(repo_root, "WES/02_data/02_processed/technical_outliers.rds"))
tes_snvs <- readRDS(file.path(repo_root, "TES/02_data/02_processed/snvs_dp.rds"))

#-------------------------------------------------------------------------------
# General parameters and functions
#-------------------------------------------------------------------------------
mutation_levels <- c("C>A","C>G","C>T","T>A","T>C","T>G")

## emmeans tidy table output
emm_pairs_tbl <- function(model, specs, by = NULL, adjust = "BH") {
  emm <- emmeans(
    model,
    specs = specs,
    by    = by,
    at    = list(log_callable = 0),   # predictions per 1 Mb
    type  = "response")

  list(
    emmeans = as_tibble(
      as.data.frame(summary(emm))),
    pairs   = as_tibble(
      as.data.frame(summary(pairs(emm), adjust = adjust, infer = TRUE))))
}

## Mutation count per sample table, for targeted exome data
tes_mutations_count <- tes_metadata %>%
  filter(
    !sample_name %in% tes_outliers$sample_name,
    !str_detect(sample_name, "ds")) %>%
  left_join(
    tes_mutations_df %>%
    filter(
    !sample_name %in% tes_outliers$sample_name,
    !str_detect(sample_name, "ds"),
    gt_AF >= 0.01) %>%
    group_by(sample_name) %>%
    summarise(
      mutations = n(),
      .groups = "drop"),
      by = "sample_name") %>%
    mutate(
      mutations    = dplyr::coalesce(mutations, 0L),
      callable_mbp = dplyr::coalesce(callable_mbp, 0),
      threshold    = callable_mbp < 0.1,
      muts_per_mbp = dplyr::if_else(
        callable_mbp > 0, mutations / callable_mbp, 0))

## Mutation count per sample table, for whole exome data
wes_mutations_count <- wes_metadata %>%
  filter(
    !sample_name %in% wes_outliers$sample_name,
    last_cycle == TRUE) %>%
  left_join(
    wes_mutations_df %>%
    filter(
      gt_AF >= 0.01,
      !sample_name %in% wes_outliers$sample_name) %>%
  group_by(sample_name) %>%
  summarise(
    mutations = n(),
    .groups = "drop"),
    by = "sample_name") %>%
  mutate(
    mutations    = dplyr::coalesce(mutations, 0L),
    callable_mbp = dplyr::coalesce(callable_mbp, 0),
    threshold    = callable_mbp < 0.1,
    muts_per_mbp = dplyr::if_else(
      callable_mbp > 0, mutations / callable_mbp, 0),
    log_coverage = log(coverage),
    log_callable = log(callable_mbp))

## 
tes_sig_df <- tes_snvs %>%
  filter(!sample_name %in% tes_outliers$sample_name, gt_AF >= 0.01) %>%
  count(
    sample_name, 
    tissue, 
    treatment, 
    time, 
    category, 
    callable_mbp, 
    mutation) %>%
  group_by(sample_name, tissue, treatment, time, category, callable_mbp) %>%
  complete(mutation = c("C>A","C>G","C>T","T>A","T>C","T>G"), fill = list(n = 0)) %>%
  mutate(total = sum(n), ratio = if_else(total > 0, n / total, 0)) %>%
  ungroup() %>%
  group_by(tissue, mutation) %>%
  mutate(
      acetone_ratio = if (any(treatment == "Acetone"))
      mean(ratio[treatment == "Acetone"], na.rm = TRUE)
      else NA_real_) %>%
  ungroup() %>%
  filter(treatment == "DT") %>%
  mutate(ratio_diff = ratio - acetone_ratio)

tes_mutations_count <- plot_levels(tes_mutations_count)
wes_mutations_count <- plot_levels(wes_mutations_count)
tes_sig_df <- plot_levels(tes_sig_df)

tes_mutations_count_stats <- tes_mutations_count %>%
  filter(
    !(tissue == "Hair follicle" & time == "Week 19" & treatment == "DT")) %>%
  mutate(
    log_callable = log(callable_mbp),
    log_coverage = log(coverage),
    weeks = case_when(
      time == "Week 8" ~ 8,
      time == "Week 14" ~ 14,
      time == "Week 17" ~ 17,
      TRUE ~ 19),
    category = droplevels(factor(category, ordered = FALSE)),
    tissue = droplevels(factor(tissue, ordered = FALSE)),
    treatment = droplevels(factor(treatment, ordered = FALSE))) %>%
  filter(is.finite(log_callable))

tes_alluvial_counts <- tes_snvs %>%
  filter(
    !sample_name %in% tes_outliers$sample_name,
    gt_AF >= 0.01,
    !(tissue == "Hair follicle" & time == "Week 19"),
    treatment == "DT") %>%
  count(sample_name, mutation) %>%
  complete(
    sample_name, 
    mutation = mutation_levels,
    fill = list(n = 0)) %>%
  left_join(
    tes_metadata,
    by = "sample_name") %>%
  pivot_wider(names_from = mutation, values_from = n) %>%
  mutate(
    category = factor(category, ordered = FALSE),
    category = droplevels(category),
    time = factor(time, ordered = FALSE),
    tissue = factor(tissue, ordered = FALSE),
    total_mutations = rowSums(across(all_of(mutation_levels))),
    log_reads = log(reads),
    log_reads_scaled = scale(log_reads),
    log_callable = log(callable_mbp),
    log_coverage = log(coverage),
    log_coverage_scaled = scale(log_coverage)) %>%
  filter(total_mutations >= 10)

TA_counts <- tes_alluvial_counts %>%
  filter(tissue == "Hair follicle", treatment == "DT") %>%
  mutate(
    time_num = as.numeric(gsub("\\D+", "", as.character(time))),
    TA = .data[["T>A"]],
    total_snvs = rowSums(across(all_of(mutation_levels))),
    Other = total_snvs - TA) 

#-------------------------------------------------------------------------------
### Mutations per Mb ###
#-------------------------------------------------------------------------------
# TES DT
tes_muts_mbp <- ggplot(
  tes_mutations_count %>% filter(treatment == "DT"),
  aes(x = treatment, y = muts_per_mbp, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
      position = position_jitterdodge(
          jitter.width = 0.35, 
          dodge.width = 0.75), 
      size = 0.5, 
      alpha = 0.6) +
  labs(
      x = "", 
      y = "Mutations per Mbp") +
  my_theme +
  scale_fill_manual(values = col_palette$tissue) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1))
# tes_glm_tissues <- MASS::glm.nb(
#   mutations ~ tissue + offset(log_callable),
#   data = tes_mutations_count_stats %>%
#     filter(treatment == "DT"))
# tes_stats_tissues <- emm_pairs_tbl(
#   tes_glm_tissues, 
#   specs = ~ tissue,
#   adjust = "none")
tes_glm_tissues_time <- MASS::glm.nb(
  mutations ~ tissue + weeks + offset(log_callable),
  data = tes_mutations_count_stats %>% 
    filter(treatment == "DT"))
tes_stats_tissues_time <- emm_pairs_tbl(
  tes_glm_tissues_time, 
  specs = ~ tissue, 
  adjust = "BH")
tes_stats_tissues_time$pairs
export_plot_data(
  data = tes_mutations_count,
  file_name = file.path(output_dir, "Figure2_TES_muts_mbp"),
  cols = c(
    "sample_name",
    "tissue",
    "treatment",
    "muts_per_mbp"))
tes_muts_mbp

# WES DT
wes_muts_mbp <- ggplot(
  wes_mutations_count %>% filter(treatment == "DT"),
  aes(x = treatment, y = muts_per_mbp, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
      position = position_jitterdodge(
          jitter.width = 0.35, 
          dodge.width = 0.75), 
      size = 0.5, 
      alpha = 0.6) +
  labs(
      x = "", 
      y = "Mutations per Mbp") +
  my_theme +
  scale_fill_manual(values = col_palette$tissue) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1))
wes_glm_tissues <- MASS::glm.nb(
  mutations ~ tissue + offset(log_callable),
  data = wes_mutations_count %>% 
    filter(treatment == "DT"))
wes_stats_tissues <- emm_pairs_tbl(
  wes_glm_tissues, 
  specs = ~ tissue, 
  adjust = "BH")
wes_stats_tissues$pairs
export_plot_data(
  data = wes_mutations_count,
  file_name = file.path(output_dir, "Figure2_WES_muts_mbp"),
  cols = c(
    "sample_name",
    "tissue",
    "treatment",
    "muts_per_mbp"))
wes_muts_mbp

#-------------------------------------------------------------------------------
### Mutations over time ###
#-------------------------------------------------------------------------------
muts_over_time <- ggplot(
  tes_mutations_count %>% 
    filter(tissue == "Hair follicle", treatment == "DT", time != "Week 19"),
  aes(x = time, y = muts_per_mbp, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
      position = position_jitterdodge(
          jitter.width = 0.35, 
          dodge.width = 0.75), 
      size = 0.5, 
      alpha = 0.6) +
  labs(
      x = "", 
      y = "Mutations per Mbp") +
  my_theme +
  scale_fill_manual(values = col_palette$tissue) +
  theme(legend.position = "")
tes_glm_weeks <- MASS::glm.nb(
  mutations ~ weeks + offset(log_callable),
  data = tes_mutations_count_stats %>% 
    filter(treatment == "DT", tissue == "Hair follicle"))
tes_stats_weeks <- emtrends(
  tes_glm_weeks,
  specs = ~ 1,
  var = "weeks")
summary(tes_stats_weeks, infer = c(TRUE, TRUE), adjust = "BH")
export_plot_data(
  data = tes_mutations_count %>% filter(tissue == "Hair follicle", treatment == "DT"),
  file_name = file.path(output_dir, "Figure2_muts_over_time"),
  cols = c(
    "sample_name",
    "tissue",
    "time",
    "muts_per_mbp"))
muts_over_time

#-------------------------------------------------------------------------------
### Mutations by morphology and time ###
#-------------------------------------------------------------------------------
muts_morph_time <- ggplot(
  tes_mutations_count %>% 
    filter(tissue == "Hair follicle", treatment == "DT", time != "Week 19") %>%
    mutate(category = factor(
      category,
      levels = c("Visually normal", "Papilloma", "SCC"))),
  aes(x = time, y = muts_per_mbp, fill = category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75), 
    size = 0.5, 
    alpha = 0.6) +
  facet_wrap(~ category) +
  labs(x = "", y = "Mutations per Mbp") +
  my_theme +
  scale_x_discrete(labels = c(
    "Week 8" = "W 8",
    "Week 14" = "W 14",
    "Week 17" = "W 17")) +
  scale_fill_manual(
    values = col_palette$category,
    breaks = c("Visually normal", "Papilloma", "SCC"),
    labels = c(
      "Visually normal" = "MN",
      Papilloma = "PP",
      SCC = "SCC")) +
  theme(
    legend.title = element_blank(),
    legend.position = "")
tes_glm_weeks_category <- MASS::glm.nb(
  mutations ~ weeks * category + offset(log_callable),
  data = tes_mutations_count_stats %>%
    filter(treatment == "DT", tissue == "Hair follicle"))
tes_stats_category_within_weeks <- emtrends(
  tes_glm_weeks_category, 
  specs = ~ category, 
  var = "weeks")
summary(tes_stats_category_within_weeks, infer = c(TRUE, TRUE), adjust = "BH")
export_plot_data(
  data = tes_mutations_count %>% filter(tissue == "Hair follicle", treatment == "DT"),
  file_name = file.path(output_dir, "Figure2_muts_morph_time"),
  cols = c(
    "sample_name",
    "tissue",
    "time",
    "category",
    "muts_per_mbp"),
  panel_width = 13)
muts_morph_time

#-------------------------------------------------------------------------------
### Signature in HF vs skin ###
#-------------------------------------------------------------------------------

sig_hf_skin <- ggplot(
  tes_sig_df %>% filter(total >= 10),
  aes(x = mutation, y = ratio_diff)) +
  geom_boxplot(
    aes(
      fill = mutation,
      group = interaction(mutation, tissue)),
      position = position_dodge2(
          width = 0.75,
          preserve = "single"),
      outlier.shape = NA) +
  geom_jitter(
    aes(color = tissue),
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75), 
    size = 0.8, 
    alpha = 0.8) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_fill_manual(values = col_palette$mutation_signature) +
  scale_color_manual(values = col_palette$tissue) +
  labs(x = "", y = "Proportion (Sample - Acetone)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  my_theme +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1))
export_plot_data(
  data = tes_sig_df %>% filter(total >= 10),
  file_name = file.path(output_dir, "Figure2_signature_hf_skin"),
  cols = c(
    "sample_name",
    "tissue",
    "treatment",
    "mutation",
    "n",
    "total",
    "ratio",
    "acetone_ratio",
    "ratio_diff"),
  panel_width = 59)
sig_hf_skin

#-------------------------------------------------------------------------------
### Signature dynamics ###
#-------------------------------------------------------------------------------
build_signature_alluvial <- function(
  data,
  x_var,
  facet_var,
  x_levels = NULL,
  facet_levels = NULL) {
  plot_df <- data %>%
    transmute(
      sample_name = sample_name,
      mutation = mutation,
      ratio = ratio,
      plot_week = .data[[x_var]],
      plot_morphology = .data[[facet_var]]) %>%
    filter(
      !is.na(plot_week),
      !is.na(plot_morphology))

  if (!is.null(x_levels)) {
    plot_df <- plot_df %>%
      mutate(plot_week = factor(plot_week, levels = x_levels))
  }

  if (!is.null(facet_levels)) {
    plot_df <- plot_df %>%
      mutate(plot_morphology = factor(plot_morphology, levels = facet_levels))
  }

  plot_df %>%
    group_by(plot_morphology, plot_week, mutation) %>%
    summarise(
      mean_ratio = mean(ratio, na.rm = TRUE),
      n_samples = n_distinct(sample_name),
      .groups = "drop") %>%
    complete(
      plot_morphology,
      plot_week,
      mutation = c("C>A","C>G","C>T","T>A","T>C","T>G"),
      fill = list(mean_ratio = 0, n_samples = 0)) %>%
    group_by(plot_morphology, plot_week) %>%
    mutate(
      total_ratio = sum(mean_ratio),
      mean_ratio = if_else(
        total_ratio > 0,
        mean_ratio / total_ratio,
        0)) %>%
      ungroup() %>%
    dplyr::select(-total_ratio) %>%
    mutate(
      mutation = factor(
        mutation,
        levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")),
        flow_id = interaction(plot_morphology, mutation, drop = TRUE))
}

alluvial_df <- tes_sig_df %>%
  filter(
    total >= 10,
    treatment == "DT",
    category != "Acetone") %>%
  mutate(
    signature_group = case_when(
      tissue == "Hair follicle" & time == "Week 8" ~ "Week 8",
      tissue == "Hair follicle" & time == "Week 14" ~ "Week 14",
      tissue == "Hair follicle" & time == "Week 17" ~ "Week 17",
      # tissue == "Hair follicle" & time == "Week 19" ~ "Week 19",
      tissue == "Skin" ~ "Skin",
      TRUE ~ NA_character_)) %>%
  build_signature_alluvial(
    x_var = "signature_group",
    facet_var = "category",
    x_levels = c("Week 8", "Week 14", "Week 17", "Skin"),
    facet_levels = c("SCC", "Papilloma", "Visually normal")) %>%
  mutate(
    plot_week_num = case_when(
      plot_week == "Week 8" ~ 1,
      plot_week == "Week 14" ~ 2,
      plot_week == "Week 17" ~ 3,
      # plot_week == "Week 19" ~ 4,
      plot_week == "Skin" ~ 3.82))

flow_df <- alluvial_df %>% filter(plot_week %in% c("Week 8", "Week 14", "Week 17"))

sig_dynamics <- ggplot(
  alluvial_df,
  aes(
    x = plot_week_num,
    y = mean_ratio,
    alluvium = flow_id,
    stratum = mutation,
    fill = mutation)) +
  geom_flow(
    data = flow_df,
    width = 0.45,
    alpha = 0.35,
    curve_type = "sigmoid",
    knot.pos = 0.35,
    segments = 100) +
  geom_stratum(
    width = 0.45,
    color = "#7A5B50",
    linewidth = 0.45) +
  facet_wrap(~ plot_morphology, nrow = 1, axes = "all_y") +
  scale_fill_manual(values = col_palette$mutation_signature, drop = FALSE) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 3.82),
    labels = c("Week 8", "Week 14", "Week 17", "Skin"),
    expand = c(0.03, 0.03)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)) +
  labs(x = "", y = "Proportion") +
  my_theme +
  theme(legend.position = "")
TA_model <- glmmTMB(
  cbind(TA, Other) ~ category * time_num + log_callable,
  family = betabinomial(link="logit"),
  data = TA_counts)
TA_slopes <- emtrends(
  TA_model,
  ~ category,
  var = "time_num")
TA_slope_pairs <- pairs(TA_slopes, adjust = "BH")
TA_slope_pairs
# BH-adjusted P values for differences in time slopes:
# Papilloma - SCC: 0.09191; Papilloma - Visually normal: 0.09191;
# SCC - Visually normal: 0.97467.
summary(TA_slopes, infer = TRUE)
TA_categories_time <- emmeans(
    TA_model, ~ category | time_num, 
    at = list(time_num = c(8,14,17)))
TA_category_pairs <- pairs(TA_categories_time, adjust = "BH")
TA_category_pairs
# BH-adjusted P values for the Papilloma - SCC contrast:
# Week 8: 0.77782; Week 14: 0.05450; Week 17: 0.004662.
export_plot_data(
  data = alluvial_df,
  file_name = file.path(output_dir, "Figure2_signature_dynamics"),
  cols = c(
    "plot_morphology",
    "plot_week",
    "mutation",
    "mean_ratio",
    "n_samples"),
  width = 12)
sig_dynamics

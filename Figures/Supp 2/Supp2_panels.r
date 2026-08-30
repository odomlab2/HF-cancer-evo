# 27/04/2026
# Author: Yoav Avi-Guy
# Purpose: Supplementary Figure 2 panels - Mutation calling from TES and WES

#-------------------------------------------------------------------------------
# Libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(ggh4x)
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
tes_unique_all <- readRDS(file.path(repo_root, "TES/02_data/02_processed/mutations_unique.rds"))
wes_unique_all <- readRDS(file.path(repo_root, "WES/02_data/02_processed/mutations_unique.rds"))
tes_unique <- readRDS(file.path(repo_root, "TES/02_data/02_processed/mutations_unique_dp.rds"))
wes_unique <- readRDS(file.path(repo_root, "WES/02_data/02_processed/mutations_unique_dp.rds"))
tes_snvs <- readRDS(file.path(repo_root, "TES/02_data/02_processed/snvs_dp.rds"))
tes_outliers <- readRDS(file.path(repo_root, "TES/02_data/02_processed/technical_outliers.rds"))
wes_outliers <- readRDS(file.path(repo_root, "WES/02_data/02_processed/technical_outliers.rds"))
gene_panel <- readRDS(file.path(repo_root, "TES/02_data/02_processed/gene_panel.rds"))

supp2_dir <- Sys.getenv(
  "SUPP2_OUTPUT_DIR",
  unset = file.path(output_root, "Figures", "Supp 2"))
dir.create(supp2_dir, recursive = TRUE, showWarnings = FALSE)

tes_metadata <- tes_metadata %>% filter(!str_detect(sample_name, "ds"))

min_vaf_muts <- 9
panel_gene_n <- n_distinct(gene_panel$gene)
variant_levels <- c("SNV", "deletion", "insertion", "substitution")
skin_vaf_group_levels <- c("TES Acetone", "TES DT", "WES DT")
skin_vaf_group_palette <- c(
  "TES Acetone" = unname(col_palette$treatment["Acetone"]),
  "TES DT" = unname(col_palette$treatment["DT"]),
  "WES DT" = unname(col_palette$tissue["Skin"]))

sample_cols <- c(
  "sample_name",
  "tissue",
  "treatment",
  "time",
  "category",
  "last_cycle",
  "callable_mbp",
  "coverage")

mutation_cols <- c(
  "sample_name",
  "tissue",
  "treatment",
  "POS",
  "time",
  "category",
  "SYMBOL",
  "gt_AF",
  "gt_DP",
  "last_cycle")

filter_last_cycle <- function(df) {
  if ("last_cycle" %in% names(df)) {
    df <- df %>% filter(last_cycle == TRUE)
  }
  df
}

is_nonsynonymous_snv <- function(consequence, impact) {
  protein_altering <- paste(
    c(
      "missense_variant",
      "frameshift_variant",
      "protein_altering_variant",
      "inframe_insertion",
      "inframe_deletion",
      "splice_acceptor_variant",
      "splice_donor_variant",
      "splice_region_variant",
      "start_lost",
      "stop_gained",
      "stop_lost",
      "coding_sequence_variant"),
    collapse = "|")

  impact %in% c("HIGH", "MODERATE") |
    str_detect(coalesce(consequence, ""), protein_altering)
}

prepare_ta_signature_plot_data <- function(data) {
  data %>%
    filter(
      treatment == "DT",
      (tissue == "Hair follicle" &
        time %in% c("Week 8", "Week 14", "Week 17")) |
        tissue == "Skin") %>%
    mutate(
      plot_group = case_when(
        tissue == "Skin" ~ "Skin",
        TRUE ~ as.character(time)),
      plot_group = factor(
        plot_group,
        levels = c("Week 8", "Week 14", "Week 17", "Skin")),
      plot_x = c(
        "Week 8" = 1,
        "Week 14" = 2,
        "Week 17" = 3,
        "Skin" = 4.35)[as.character(plot_group)])
}

set_ta_signature_levels <- function(data) {
  data %>%
    plot_levels() %>%
    mutate(
      time = case_when(
        tissue == "Skin" ~ "Skin 19",
        TRUE ~ as.character(time)),
      time = factor(
        time,
        levels = c("Week 8", "Week 14", "Week 17", "Week 19", "Skin 19"),
        ordered = TRUE))
}

plot_ta_signature_over_time <- function(
  data,
  y_var,
  y_label,
  y_cap = NULL) {

  y_var <- enquo(y_var)
  point_data <- data
  if (!is.null(y_cap)) {
    point_data <- point_data %>%
      filter((!!y_var) <= y_cap)
  }

  plot <- ggplot(
    data,
    aes(
      x = plot_x,
      y = !!y_var,
      fill = tissue,
      group = plot_group)) +
    geom_boxplot(
      width = 0.65,
      outlier.shape = NA) +
    geom_jitter(
      data = point_data,
      width = 0.18,
      height = 0,
      size = 0.5,
      alpha = 0.6) +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4.35),
      labels = c("W 8", "W 14", "W 17", "Skin"),
      expand = expansion(add = 0.5)) +
    scale_fill_manual(values = col_palette$tissue) +
    labs(
      x = "",
      y = y_label) +
    my_theme +
    theme(legend.position = "none")

  if (is.null(y_cap)) {
    return(plot)
  }

  capped_data <- data %>%
    filter((!!y_var) > y_cap)

  plot +
    geom_jitter(
      data = capped_data,
      aes(
        x = plot_x,
        y = y_cap,
        group = plot_group),
      inherit.aes = FALSE,
      width = 0.18,
      height = 0,
      shape = 17,
      size = 1.2) +
    coord_cartesian(
      ylim = c(0, y_cap),
      clip = "off")
}

emm_pairs_tbl <- function(model, specs, by = NULL, adjust = "BH") {
  emm <- emmeans::emmeans(
    model,
    specs = specs,
    by = by,
    at = list(log_callable = 0),
    type = "response")

  list(
    emmeans = as_tibble(as.data.frame(summary(emm))),
    pairs = as_tibble(as.data.frame(summary(
      pairs(emm),
      adjust = adjust,
      infer = TRUE))))
}

count_sample_mutations <- function(metadata, mutations, outliers) {
  metadata %>%
    filter(!sample_name %in% outliers$sample_name) %>%
    select(any_of(sample_cols)) %>%
    distinct() %>%
    left_join(
      mutations %>%
        filter(
          !sample_name %in% outliers$sample_name,
          gt_AF >= 0.01) %>%
        count(sample_name, name = "mutations"),
      by = "sample_name") %>%
    mutate(
      mutations = coalesce(mutations, 0L),
      log_callable = log(callable_mbp),
      log_coverage = log(coverage)) %>%
    plot_levels()
}

count_ta_mutations <- function(metadata, snvs, outliers) {
  metadata %>%
    filter(!sample_name %in% outliers$sample_name) %>%
    left_join(
      snvs %>%
        filter(
          !sample_name %in% outliers$sample_name,
          mutation == "T>A",
          gt_AF >= 0.01) %>%
        count(sample_name, name = "mutations"),
      by = "sample_name") %>%
    mutate(
      mutations = coalesce(mutations, 0L),
      callable_mbp = coalesce(callable_mbp, 0),
      log_callable = log(callable_mbp),
      log_coverage = log(coverage),
      threshold = callable_mbp < 0.1,
      muts_per_mbp = if_else(
        callable_mbp > 0,
        mutations / callable_mbp,
        0)) %>%
    plot_levels()
}

prepare_ta_mutations_stats <- function(data) {
  data %>%
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
      category = droplevels(factor(category, ordered = FALSE))) %>%
    filter(is.finite(log_callable))
}

prepare_mutation_depth <- function(mutations, outliers) {
  mutations %>%
    filter(!sample_name %in% outliers$sample_name) %>%
    mutate(gt_DP = as.numeric(gt_DP)) %>%
    plot_levels()
}

bin_mutation_depth <- function(mutations) {
  depth_levels <- c(
    as.character(seq(1, 14)),
    "15-49",
    "50-99",
    "100-499",
    "500-1000",
    "1000+")

  mutations %>%
    mutate(
      gt_dp_bins = case_when(
        gt_DP > 1000 ~ "1000+",
        gt_DP > 499 ~ "500-1000",
        gt_DP > 99 ~ "100-499",
        gt_DP > 49 ~ "50-99",
        gt_DP > 14 ~ "15-49",
        gt_DP < 15 ~ as.character(gt_DP),
        TRUE ~ NA_character_),
      gt_dp_bins = factor(
        gt_dp_bins,
        levels = depth_levels,
        ordered = TRUE))
}

summarise_mutation_depth_bins <- function(mutations) {
  mutations %>%
    filter(!is.na(gt_dp_bins)) %>%
    count(gt_dp_bins, tissue, name = "n_mutations") %>%
    mutate(pct_mutations = n_mutations / sum(n_mutations) * 100)
}

summarise_mutation_depth_values <- function(mutations) {
  mutations %>%
    filter(!is.na(gt_dp_bins), !is.na(gt_DP)) %>%
    count(gt_DP, gt_dp_bins, tissue, name = "n_mutations") %>%
    mutate(pct_mutations = n_mutations / sum(n_mutations) * 100)
}

sample_vaf_density <- function(vaf, from = 0, to = 1, n = 512) {
  vaf <- as.numeric(vaf)
  vaf <- vaf[is.finite(vaf)]
  if (length(unique(vaf)) < 2) {
    return(tibble(x = seq(from, to, length.out = n), density = 0))
  }
  d <- stats::density(vaf, from = from, to = to, n = n)
  tibble(x = d$x, density = d$y)
}

vaf_density_by_group <- function(
  df,
  group_vars,
  min_mut = min_vaf_muts,
  from = 0,
  to = 1,
  n = 512,
  ci_level = 0.95) {

  df_filt <- df %>%
    mutate(gt_AF = as.numeric(gt_AF)) %>%
    filter(is.finite(gt_AF)) %>%
    group_by(across(all_of(group_vars)), sample_name) %>%
    filter(n() >= min_mut) %>%
    ungroup()

  sample_dens <- df_filt %>%
    group_by(across(all_of(group_vars)), sample_name) %>%
    summarise(
      n_mut = n(),
      dens = list(sample_vaf_density(gt_AF, from = from, to = to, n = n)),
      .groups = "drop") %>%
    tidyr::unnest(dens)

  alpha <- 1 - ci_level

  group_summary <- sample_dens %>%
    group_by(across(all_of(group_vars)), x) %>%
    summarise(
      n_samples = n_distinct(sample_name),
      mean = mean(density, na.rm = TRUE),
      sd = stats::sd(density, na.rm = TRUE),
      se = sd / sqrt(n_samples),
      tcrit = stats::qt(1 - alpha / 2, df = pmax(n_samples - 1, 1)),
      lwr = mean - tcrit * se,
      upr = mean + tcrit * se,
      .groups = "drop") %>%
    mutate(
      lwr = if_else(n_samples >= 2, pmax(lwr, 0), NA_real_),
      upr = if_else(n_samples >= 2, upr, NA_real_))

  mode_per_group <- group_summary %>%
    group_by(across(all_of(group_vars))) %>%
    slice_max(order_by = mean, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      across(all_of(group_vars)),
      mode_x = x,
      mode_density = mean)

  list(
    sample_densities = sample_dens,
    group_summary = group_summary,
    mode_per_group = mode_per_group)
}

build_skin_vaf <- function(mutations, outliers) {
  vaf_density_by_group(
    mutations %>%
      filter(
        !sample_name %in% outliers$sample_name,
        tissue == "Skin",
        treatment %in% c("Acetone", "DT")) %>%
      plot_levels(),
    group_vars = "treatment")
}

build_skin_vaf_combined <- function(
  tes_mutations,
  tes_outliers,
  wes_mutations,
  wes_outliers) {

  combined <- bind_rows(
    tes_mutations %>%
      filter(
        !sample_name %in% tes_outliers$sample_name,
        tissue == "Skin",
        treatment %in% c("Acetone", "DT")) %>%
      mutate(
        sequencing = "TES",
        vaf_group = paste(sequencing, treatment)),
    wes_mutations %>%
      filter(
        !sample_name %in% wes_outliers$sample_name,
        tissue == "Skin",
        treatment == "DT") %>%
      mutate(
        sequencing = "WES",
        vaf_group = paste(sequencing, treatment))) %>%
    mutate(
      vaf_group = factor(
        vaf_group,
        levels = skin_vaf_group_levels,
        ordered = TRUE))

  vaf_density_by_group(
    combined,
    group_vars = "vaf_group")
}

add_vaf_read_counts <- function(data) {
  data <- data %>%
    mutate(
      gt_AF = as.numeric(gt_AF),
      gt_DP = as.integer(gt_DP))

  if ("gt_AD" %in% names(data)) {
    data <- data %>%
      tidyr::separate(
        col = gt_AD,
        into = c("ref_reads", "alt_reads"),
        sep = ",",
        remove = FALSE,
        convert = TRUE)
  } else {
    data <- data %>%
      mutate(
        alt_reads = as.integer(round(gt_AF * gt_DP)),
        ref_reads = gt_DP - alt_reads)
  }

  data %>%
    mutate(
      alt_reads = as.integer(alt_reads),
      ref_reads = as.integer(ref_reads))
}

build_skin_vaf_model_df <- function(
  tes_mutations,
  tes_metadata,
  tes_outliers,
  wes_mutations,
  wes_metadata,
  wes_outliers,
  min_mut = min_vaf_muts) {

  coverage_df <- bind_rows(
    tes_metadata %>%
      transmute(
        sequencing = "TES",
        sample_name,
        coverage = as.numeric(coverage)),
    wes_metadata %>%
      transmute(
        sequencing = "WES",
        sample_name,
        coverage = as.numeric(coverage))) %>%
    distinct(sequencing, sample_name, .keep_all = TRUE)

  bind_rows(
    tes_mutations %>%
      filter(
        !sample_name %in% tes_outliers$sample_name,
        tissue == "Skin",
        treatment %in% c("Acetone", "DT")) %>%
      mutate(
        sequencing = "TES",
        vaf_group = paste(sequencing, treatment)),
    wes_mutations %>%
      filter(
        !sample_name %in% wes_outliers$sample_name,
        tissue == "Skin",
        treatment == "DT") %>%
      mutate(
        sequencing = "WES",
        vaf_group = paste(sequencing, treatment))) %>%
    select(-any_of("coverage")) %>%
    left_join(coverage_df, by = c("sequencing", "sample_name")) %>%
    add_vaf_read_counts() %>%
    filter(
      is.finite(gt_AF),
      is.finite(gt_DP),
      gt_DP > 0,
      alt_reads >= 0,
      ref_reads >= 0,
      alt_reads + ref_reads > 0,
      is.finite(coverage),
      coverage > 0,
      is.finite(callable_mbp),
      callable_mbp > 0) %>%
    group_by(vaf_group, sample_name) %>%
    filter(n() >= min_mut) %>%
    ungroup() %>%
    mutate(
      vaf_group = factor(
        vaf_group,
        levels = skin_vaf_group_levels,
        ordered = FALSE),
      log_coverage = log(coverage),
      log_callable = log(callable_mbp))
}

count_tes_panel_genes <- function(mutations, outliers) {
  mutations %>%
    filter(
      !sample_name %in% outliers$sample_name,
      gt_AF >= 0.01,
      !is.na(SYMBOL),
      SYMBOL %in% gene_panel$gene) %>%
    distinct(tissue, treatment, SYMBOL) %>%
    count(tissue, treatment, name = "n_genes_mutated") %>%
    mutate(pct_panel_genes = n_genes_mutated / panel_gene_n * 100) %>%
    plot_levels()
}

count_wes_genes <- function(mutations, outliers) {
  mutations %>%
    filter(
      !sample_name %in% outliers$sample_name,
      gt_AF >= 0.01,
      !is.na(SYMBOL)) %>%
    distinct(tissue, treatment, SYMBOL) %>%
    count(tissue, treatment, name = "n_genes_mutated") %>%
    plot_levels()
}

plot_sample_boxplot <- function(data, y_var, y_label) {
  ggplot(
    data,
    aes(
      x = treatment,
      y = {{ y_var }},
      fill = treatment)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.35,
        dodge.width = 0.75),
      size = 0.5,
      alpha = 0.6) +
    facet_wrap(~ tissue, scales = "free_y") +
    scale_fill_manual(values = col_palette$treatment) +
    scale_y_continuous(labels = my_label) +
    labs(
      x = "",
      y = y_label) +
    my_theme +
    theme(legend.position = "none") +
    expand_limits(y = 0)
}

plot_mutation_depth_bins <- function(data) {
  ggplot(
    data,
    aes(
      x = gt_dp_bins,
      y = pct_mutations,
      fill = tissue)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = col_palette$tissue) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "# of reads (bin)",
      y = "Mutations (%)") +
    my_theme +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5)) +
    expand_limits(y = 0)
}

summarise_mutation_types <- function(mutations, outliers) {
  if (!"time" %in% names(mutations)) {
    mutations <- mutations %>% mutate(time = NA_character_)
  }

  mutations %>%
    filter(
      !sample_name %in% outliers$sample_name,
      gt_AF >= 0.01,
      VARIANT_CLASS %in% variant_levels) %>%
    count(
      sample_name,
      tissue,
      treatment,
      time,
      category,
      VARIANT_CLASS,
      name = "n") %>%
    group_by(sample_name, tissue, treatment, time, category) %>%
    tidyr::complete(
      VARIANT_CLASS = variant_levels,
      fill = list(n = 0)) %>%
    mutate(
      total = sum(n),
      ratio = if_else(total > 0, n / total, 0)) %>%
    ungroup() %>%
    mutate(
      VARIANT_CLASS = factor(
        VARIANT_CLASS,
        levels = variant_levels,
        ordered = TRUE)) %>%
    plot_levels()
}

summarise_mutation_type_bars <- function(
  data,
  x_var,
  facet_var,
  x_levels = NULL,
  facet_levels = NULL) {
  plot_df <- data %>%
    transmute(
      sample_name = sample_name,
      VARIANT_CLASS = VARIANT_CLASS,
      ratio = ratio,
      plot_group = .data[[x_var]],
      plot_morphology = .data[[facet_var]]) %>%
    filter(
      !is.na(plot_group),
      !is.na(plot_morphology))

  if (!is.null(x_levels)) {
    plot_df <- plot_df %>%
      mutate(plot_group = factor(plot_group, levels = x_levels))
  }

  if (!is.null(facet_levels)) {
    plot_df <- plot_df %>%
      mutate(plot_morphology = factor(plot_morphology, levels = facet_levels))
  }

  out <- plot_df %>%
    group_by(plot_morphology, plot_group, VARIANT_CLASS) %>%
    summarise(
      mean_ratio = mean(ratio, na.rm = TRUE),
      sd_ratio = stats::sd(ratio, na.rm = TRUE),
      n_samples = n_distinct(sample_name),
      .groups = "drop") %>%
    mutate(
      sd_ratio = if_else(is.na(sd_ratio), 0, sd_ratio),
      se_ratio = sd_ratio / sqrt(n_samples),
      VARIANT_CLASS = factor(
        VARIANT_CLASS,
        levels = variant_levels,
        ordered = TRUE)) %>%
    group_by(plot_morphology, plot_group) %>%
    mutate(
      total_ratio = sum(mean_ratio),
      mean_ratio = if_else(total_ratio > 0, mean_ratio / total_ratio, 0)) %>%
    ungroup() %>%
    select(-total_ratio)

  if (!is.null(x_levels)) {
    out <- out %>% mutate(plot_group = factor(plot_group, levels = x_levels))
  }

  if (!is.null(facet_levels)) {
    out <- out %>%
      mutate(plot_morphology = factor(plot_morphology, levels = facet_levels))
  }

  out %>%
    arrange(plot_morphology, plot_group, VARIANT_CLASS) %>%
    group_by(plot_morphology, plot_group) %>%
    mutate(
      stack_top = cumsum(mean_ratio),
      stack_bottom = stack_top - mean_ratio,
      stack_mid = stack_bottom + mean_ratio / 2,
      err_ymin = pmax(stack_bottom, stack_mid - se_ratio),
      err_ymax = pmin(stack_top, stack_mid + se_ratio)) %>%
    ungroup()
}

wes_unique_last <- filter_last_cycle(wes_unique)
wes_unique_all_last <- filter_last_cycle(wes_unique_all)
wes_metadata_last <- filter_last_cycle(wes_metadata)

tes_muts <- count_sample_mutations(
  metadata = tes_metadata,
  mutations = tes_unique,
  outliers = tes_outliers)
wes_muts <- count_sample_mutations(
  metadata = wes_metadata_last,
  mutations = wes_unique_last,
  outliers = wes_outliers)
tes_ta_muts <- count_ta_mutations(
  metadata = tes_metadata,
  snvs = tes_snvs,
  outliers = tes_outliers)
tes_ta_muts_stats <- prepare_ta_mutations_stats(tes_ta_muts)
tes_ta_signature <- tes_snvs %>%
  filter(
    !sample_name %in% tes_outliers$sample_name,
    gt_AF >= 0.01) %>%
  count(sample_name, tissue, treatment, time, category, mutation) %>%
  group_by(sample_name, tissue, treatment, time, category) %>%
  tidyr::complete(
    mutation = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
    fill = list(n = 0)) %>%
  mutate(
    total_snvs = sum(n),
    pct_ta = if_else(total_snvs > 0, 100 * n / total_snvs, 0)) %>%
  ungroup() %>%
  filter(mutation == "T>A", total_snvs >= 10) %>%
  set_ta_signature_levels()
tes_nonsynonymous_snvs <- tes_snvs %>%
  filter(
    !sample_name %in% tes_outliers$sample_name,
    gt_AF >= 0.01,
    !(
      (treatment == "Acetone" |
        (treatment == "DT" & category == "Visually normal")) &
        IMPACT %in% c("LOW", "MODIFIER")),
    !is.na(SYMBOL),
    SYMBOL != "",
    Consequence != "intergenic_variant",
    Feature_type != "RegulatoryFeature",
    !str_detect(coalesce(Consequence, ""), "synonymous_variant"),
    is_nonsynonymous_snv(Consequence, IMPACT))
tes_nonsynonymous_burden <- tes_metadata %>%
  filter(!sample_name %in% tes_outliers$sample_name) %>%
  select(any_of(sample_cols)) %>%
  distinct() %>%
  left_join(
    tes_nonsynonymous_snvs %>%
      distinct(sample_name, CHROM, POS, REF, ALT) %>%
      count(sample_name, name = "n_nonsynonymous_mutations"),
    by = "sample_name") %>%
  mutate(
    n_nonsynonymous_mutations = coalesce(
      n_nonsynonymous_mutations,
      0L),
    nonsynonymous_mutations_per_mbp =
      n_nonsynonymous_mutations / callable_mbp) %>%
  set_ta_signature_levels()
tes_ta_signature_stats <- tes_ta_signature %>%
  filter(tissue == "Hair follicle", treatment == "DT", time != "Week 19") %>%
  left_join(
    tes_ta_muts_stats %>%
      select(sample_name, callable_mbp, log_callable, weeks),
    by = "sample_name") %>%
  filter(is.finite(log_callable)) %>%
  mutate(category = droplevels(factor(category, ordered = FALSE)))

tes_depth <- prepare_mutation_depth(
  mutations = tes_unique_all,
  outliers = tes_outliers)
wes_depth <- prepare_mutation_depth(
  mutations = wes_unique_all_last,
  outliers = wes_outliers)
tes_depth_binned <- bin_mutation_depth(tes_depth)
wes_depth_binned <- bin_mutation_depth(wes_depth)
tes_depth_bins <- summarise_mutation_depth_bins(tes_depth_binned)
wes_depth_bins <- summarise_mutation_depth_bins(wes_depth_binned)
tes_depth_values <- summarise_mutation_depth_values(tes_depth_binned)
wes_depth_values <- summarise_mutation_depth_values(wes_depth_binned)

skin_vaf <- build_skin_vaf_combined(
  tes_mutations = tes_unique,
  tes_outliers = tes_outliers,
  wes_mutations = wes_unique_last,
  wes_outliers = wes_outliers)
skin_vaf_model_df <- build_skin_vaf_model_df(
  tes_mutations = tes_unique,
  tes_metadata = tes_metadata,
  tes_outliers = tes_outliers,
  wes_mutations = wes_unique_last,
  wes_metadata = wes_metadata_last,
  wes_outliers = wes_outliers)

tes_panel_genes <- count_tes_panel_genes(
  mutations = tes_unique,
  outliers = tes_outliers)
wes_genes <- count_wes_genes(
  mutations = wes_unique_last,
  outliers = wes_outliers)

tes_mutation_types <- summarise_mutation_types(
  mutations = tes_unique,
  outliers = tes_outliers)
wes_mutation_types <- summarise_mutation_types(
  mutations = wes_unique_last,
  outliers = wes_outliers)

#-------------------------------------------------------------------------------
### Raw mutation counts ###
#-------------------------------------------------------------------------------

TES_raw_muts <- plot_sample_boxplot(
  data = tes_muts,
  y_var = mutations,
  y_label = "# of mutations")
export_plot_data(
  data = tes_muts,
  file_name = file.path(supp2_dir, "Supp2_TES_raw_mutation_count"),
  cols = c(sample_cols, "mutations"))
TES_raw_muts

WES_raw_muts <- plot_sample_boxplot(
  data = wes_muts,
  y_var = mutations,
  y_label = "# of mutations")
export_plot_data(
  data = wes_muts,
  file_name = file.path(supp2_dir, "Supp2_WES_raw_mutation_count"),
  cols = c(sample_cols, "mutations"))
wes_raw_glm_tissues <- MASS::glm.nb(
  mutations ~ tissue * treatment + offset(log_callable),
  data = wes_muts)
wes_raw_stats_tissues <- emm_pairs_tbl(
  wes_raw_glm_tissues,
  specs = ~ tissue,
  by = "treatment",
  adjust = "none")
WES_raw_muts

#-------------------------------------------------------------------------------
### T>A mutations by morphology and time ###
#-------------------------------------------------------------------------------

TA_muts_morph_time <- ggplot(
  tes_ta_muts %>%
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
  labs(x = "", y = "T>A mutations per Mbp") +
  my_theme +
  scale_x_discrete(labels = c(
    "Week 8" = "W 8",
    "Week 14" = "W 14",
    "Week 17" = "W 17")) +
  scale_fill_manual(
    values = col_palette$category,
    breaks = c("Visually normal", "Papilloma", "SCC")) +
  facet_wrap(~ category, nrow = 1) +
  theme(legend.position = "none")
tes_ta_glm_weeks_category <- MASS::glm.nb(
  mutations ~ weeks * category + offset(log_callable),
  data = tes_ta_muts_stats)
tes_ta_stats_category_within_weeks <- emmeans::emtrends(
  tes_ta_glm_weeks_category,
  specs = ~ category,
  var = "weeks")
summary(tes_ta_stats_category_within_weeks, infer = c(TRUE, TRUE), adjust = "holm")
export_plot_data(
  data = tes_ta_muts %>% filter(tissue == "Hair follicle", treatment == "DT"),
  file_name = file.path(supp2_dir, "Supp2_TA_muts_morph_time"),
  cols = c(
    "sample_name",
    "tissue",
    "time",
    "category",
    "muts_per_mbp"),
  panel_width = 13)
TA_muts_morph_time

#-------------------------------------------------------------------------------
### T>A signature proportion over time ###
#-------------------------------------------------------------------------------

tes_ta_signature_plot <- prepare_ta_signature_plot_data(tes_ta_signature)
TA_signature_over_time <- plot_ta_signature_over_time(
  data = tes_ta_signature_plot,
  y_var = pct_ta,
  y_label = "T>A mutations (%)")
export_plot_data(
  data = tes_ta_signature_plot,
  file_name = file.path(supp2_dir, "Supp2_TA_signature_over_time"),
  cols = c(
    "sample_name",
    "tissue",
    "time",
    "category",
    "n",
    "total_snvs",
    "pct_ta"))
tes_nonsynonymous_burden_plot <-
  prepare_ta_signature_plot_data(tes_nonsynonymous_burden)
TA_nonsynonymous_signature_over_time <- plot_ta_signature_over_time(
  data = tes_nonsynonymous_burden_plot,
  y_var = nonsynonymous_mutations_per_mbp,
  y_label = "Nonsynonymous mutations / Mb",
  y_cap = 20)
export_plot_data(
  data = tes_nonsynonymous_burden_plot,
  file_name = file.path(
    supp2_dir,
    "Supp2_TA_nonsynonymous_signature_over_time"),
  cols = c(
    "sample_name",
    "tissue",
    "time",
    "category",
    "n_nonsynonymous_mutations",
    "callable_mbp",
    "nonsynonymous_mutations_per_mbp"))
tes_ta_signature_glm_weeks_overall <- MASS::glm.nb(
  n ~ weeks + offset(log_callable),
  data = tes_ta_signature_stats)
tes_ta_signature_stats_overall_weeks <- emmeans::emtrends(
  tes_ta_signature_glm_weeks_overall,
  specs = ~ 1,
  var = "weeks")
summary(
  tes_ta_signature_stats_overall_weeks,
  infer = c(TRUE, TRUE))
TA_signature_over_time
TA_nonsynonymous_signature_over_time

#-------------------------------------------------------------------------------
### Informative reads per mutation ###
#-------------------------------------------------------------------------------

TES_informative_reads <- plot_mutation_depth_bins(tes_depth_bins)
export_plot_data(
  data = tes_depth_values,
  file_name = file.path(supp2_dir, "Supp2_TES_informative_reads"),
  cols = c("gt_DP", "gt_dp_bins", "tissue", "n_mutations", "pct_mutations"))
TES_informative_reads

WES_informative_reads <- plot_mutation_depth_bins(wes_depth_bins)
export_plot_data(
  data = wes_depth_values,
  file_name = file.path(supp2_dir, "Supp2_WES_informative_reads"),
  cols = c("gt_DP", "gt_dp_bins", "tissue", "n_mutations", "pct_mutations"))
WES_informative_reads

#-------------------------------------------------------------------------------
### Skin VAF distributions ###
#-------------------------------------------------------------------------------

skin_vaf_labels <- skin_vaf$group_summary %>%
  group_by(vaf_group) %>%
  summarise(n_samples = max(n_samples), .groups = "drop") %>%
  arrange(vaf_group) %>%
  mutate(
    x = 0.5,
    y = Inf,
    label = paste0(vaf_group, ": n = ", n_samples),
    vjust = row_number() * 1.2)

skin_vaf_plot <- ggplot(
  skin_vaf$group_summary,
  aes(
    x = x,
    y = mean,
    colour = vaf_group,
    fill = vaf_group)) +
  geom_ribbon(
    data = ~ filter(.x, !is.na(lwr), !is.na(upr)),
    aes(ymin = lwr, ymax = upr),
    alpha = 0.2,
    colour = NA) +
  geom_line(linewidth = 1) +
  geom_vline(
    data = skin_vaf$mode_per_group,
    aes(xintercept = mode_x, colour = vaf_group),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.7,
    alpha = 0.65) +
  geom_text(
    data = skin_vaf_labels,
    aes(
      x = x,
      y = y,
      label = label,
      color = vaf_group,
      vjust = vjust),
    inherit.aes = FALSE,
    hjust = 1.05,
    size = 2.6,
    show.legend = FALSE) +
  scale_colour_manual(values = skin_vaf_group_palette) +
  scale_fill_manual(values = skin_vaf_group_palette) +
  coord_cartesian(xlim = c(0, 0.5)) +
  labs(
    x = "VAF",
    y = "Density") +
  my_theme +
  theme(legend.title = element_blank())
export_plot_data(
  data = skin_vaf$group_summary,
  file_name = file.path(supp2_dir, "Supp2_skin_VAF_distribution"),
  cols = c("vaf_group", "x", "mean", "lwr", "upr", "n_samples"))
skin_vaf_fit <- glmmTMB(
  cbind(alt_reads, ref_reads) ~ vaf_group + log_callable,
  family = betabinomial(link = "logit"),
  data = skin_vaf_model_df)
skin_vaf_emmeans <- emmeans(
  skin_vaf_fit,
  specs = ~ vaf_group,
  type = "response")
skin_vaf_stats_pairs <- as_tibble(as.data.frame(summary(
  pairs(skin_vaf_emmeans),
  adjust = "BH",
  infer = TRUE)))
skin_vaf_plot

#-------------------------------------------------------------------------------
### Genes mutated ###
#-------------------------------------------------------------------------------

TES_panel_genes_plot <- ggplot(
  tes_panel_genes,
  aes(
    x = treatment,
    y = pct_panel_genes,
    fill = treatment)) +
  geom_col(width = 0.65) +
  facet_wrap(~ tissue) +
  scale_fill_manual(values = col_palette$treatment) +
  scale_y_continuous(
    name = "% of panel genes mutated",
    sec.axis = sec_axis(
      ~ . * panel_gene_n / 100,
      name = "# of panel genes mutated",
      breaks = c(0, panel_gene_n),
      labels = function(x) round(x))) +
  labs(x = "") +
  my_theme +
  theme(legend.position = "none") +
  expand_limits(y = 0)
export_plot_data(
  data = tes_panel_genes,
  file_name = file.path(supp2_dir, "Supp2_TES_panel_genes_mutated"),
  cols = c(
    "tissue",
    "treatment",
    "n_genes_mutated",
    "pct_panel_genes"))
TES_panel_genes_plot

WES_genes_plot <- ggplot(
  wes_genes,
  aes(
    x = treatment,
    y = n_genes_mutated,
    fill = treatment)) +
  geom_col(width = 0.65) +
  facet_wrap(~ tissue, scales = "free_y") +
  scale_fill_manual(values = col_palette$treatment) +
  scale_y_continuous(labels = my_label) +
  labs(
    x = "",
    y = "# of genes mutated") +
  my_theme +
  theme(legend.position = "none") +
  expand_limits(y = 0)
export_plot_data(
  data = wes_genes,
  file_name = file.path(supp2_dir, "Supp2_WES_genes_mutated"),
  cols = c(
    "tissue",
    "treatment",
    "n_genes_mutated"))
WES_genes_plot

#-------------------------------------------------------------------------------
### Mutation type distributions ###
#-------------------------------------------------------------------------------

tes_mutation_type_bars <- tes_mutation_types %>%
  mutate(
    plot_group = case_when(
      treatment == "Acetone" & tissue == "Hair follicle" ~ "Acetone HF",
      treatment == "Acetone" & tissue == "Skin" ~ "Acetone Skin",
      treatment == "DT" & tissue == "Hair follicle" ~ as.character(time),
      treatment == "DT" & tissue == "Skin" ~ "Skin DT",
      TRUE ~ NA_character_),
    plot_morphology = case_when(
      treatment == "Acetone" ~ "Acetone",
      TRUE ~ as.character(category))) %>%
  summarise_mutation_type_bars(
    x_var = "plot_group",
    facet_var = "plot_morphology",
    x_levels = c(
      "Acetone HF",
      "Week 8",
      "Week 14",
      "Week 17",
      "Week 19",
      "Skin DT",
      "Acetone Skin"),
    facet_levels = c("Acetone", "SCC", "Papilloma", "Visually normal"))

TES_mutation_type_dynamics <- ggplot(
  tes_mutation_type_bars,
  aes(
    x = plot_group,
    y = mean_ratio,
    fill = VARIANT_CLASS)) +
  geom_col(
    width = 0.65,
    position = position_stack(reverse = TRUE)) +
  geom_errorbar(
    aes(ymin = err_ymin, ymax = err_ymax),
    width = 0.18,
    linewidth = 0.25,
    color = "black") +
  facet_wrap(~ plot_morphology, nrow = 1, scales = "free_x", axes = "all_y") +
  scale_fill_manual(values = col_palette$mutation_type, drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)) +
  labs(
    x = "",
    y = "Proportion") +
  my_theme +
  theme(
    legend.title = element_blank())
export_plot_data(
  data = tes_mutation_type_bars,
  file_name = file.path(supp2_dir, "Supp2_TES_mutation_type_dynamics"),
  cols = c(
    "plot_morphology",
    "plot_group",
    "VARIANT_CLASS",
    "mean_ratio",
    "sd_ratio",
    "se_ratio",
    "n_samples"),
  width = 60)
TES_mutation_type_dynamics

wes_mutation_type_bars <- wes_mutation_types %>%
  mutate(
    plot_tissue = tissue,
    plot_morphology = category) %>%
  summarise_mutation_type_bars(
    x_var = "plot_tissue",
    facet_var = "plot_morphology",
    x_levels = c("Hair follicle", "Skin"),
    facet_levels = c("Acetone", "SCC", "Papilloma", "Visually normal"))

WES_mutation_type_tissue_category <- ggplot(
  wes_mutation_type_bars,
  aes(
    x = plot_group,
    y = mean_ratio,
    fill = VARIANT_CLASS)) +
  geom_col(
    width = 0.65,
    position = position_stack(reverse = TRUE)) +
  geom_errorbar(
    aes(ymin = err_ymin, ymax = err_ymax),
    width = 0.18,
    linewidth = 0.25,
    color = "black") +
  facet_wrap(~ plot_morphology, nrow = 1, scales = "free_x", axes = "all_y") +
  scale_fill_manual(values = col_palette$mutation_type, drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)) +
  labs(
    x = "",
    y = "Proportion") +
  my_theme +
  theme(
    legend.title = element_blank())
export_plot_data(
  data = wes_mutation_type_bars,
  file_name = file.path(supp2_dir, "Supp2_WES_mutation_type_tissue_category"),
  cols = c(
    "plot_morphology",
    "plot_group",
    "VARIANT_CLASS",
    "mean_ratio",
    "sd_ratio",
    "se_ratio",
    "n_samples"),
  width = 50)
WES_mutation_type_tissue_category

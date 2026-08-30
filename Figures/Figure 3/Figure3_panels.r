# 12/04/2026
# Author: Yoav Avi-Guy
# Purpose: Figure 3 panels - Hair follicles and skin share mutations that help 
# trace clonal evolution in tumor-prone regions
#-------------------------------------------------------------------------------
# Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(UpSetR)
library(emmeans)
library(glmmTMB)
library(ggbeeswarm)

# Colours and theme
repo_root <- Sys.getenv("HF_SCC_ROOT", unset = "")
if (!nzchar(repo_root)) stop("Set HF_SCC_ROOT to the repository root.")
repo_root <- normalizePath(repo_root, mustWork = TRUE)
figure3_source_dir <- file.path(repo_root, "Figures/Figure 3")
source(file.path(repo_root, "TES/01_source/01_plotting.r"))
source(file.path(repo_root, "scripts/external_data_mode.R"))
EXTERNAL_DATA_MODE <- publication_external_data_mode()
figure3_dir <- dirname(publication_output_path(
  file.path(figure3_source_dir, ".output-root")))
options(device = function(...) grDevices::pdf(
  file = publication_output_path(file.path(figure3_dir, "Rplots.pdf")), ...))

# Loading files
tes_metadata <- readRDS(file.path(
  repo_root, "TES/02_data/01_metadata/sample_metadata.rds"))
tes_mutations_df <- readRDS(file.path(
  repo_root, "TES/02_data/02_processed/mutations_unique_dp.rds"))
tes_outliers <- readRDS(file.path(
  repo_root, "TES/02_data/02_processed/technical_outliers.rds"))
wes_hf_skin <- readRDS(file.path(
  repo_root, "WES/02_data/02_processed/shared_muts_hf_skin.rds"))
refcds_path <- Sys.getenv("HF_SCC_REFCDS", unset = "")
if (!nzchar(refcds_path)) stop("Set HF_SCC_REFCDS to the reference CDS RDA file.")
load(refcds_path)

export_grid_plot_data <- function(
  plot_grob,
  data,
  file_name,
  width,
  height,
  cols = NULL) {

  out <- if (is.null(cols)) data else dplyr::select(data, any_of(cols))
  data_path <- publication_output_path(paste0(file_name, "_data.csv"))
  pdf_path <- publication_output_path(paste0(file_name, ".pdf"))
  readr::write_csv(out, data_path)

  grDevices::pdf(pdf_path, width = width, height = height)
  grid::grid.draw(plot_grob)
  dev.off()

  invisible(out)
}

export_pdf_plot_data <- function(
  data,
  file_name,
  draw,
  width,
  height,
  cols = NULL) {

  out <- if (is.null(cols)) data else dplyr::select(data, any_of(cols))
  data_path <- publication_output_path(paste0(file_name, "_data.csv"))
  pdf_path <- publication_output_path(paste0(file_name, ".pdf"))
  readr::write_csv(out, data_path)

  grDevices::pdf(pdf_path, width = width, height = height)
  draw()
  dev.off()

  invisible(out)
}

export_ggplot_data <- function(
  plot,
  file_name,
  data,
  cols = NULL,
  width = NULL,
  height = NULL,
  panel_width = 29,
  panel_height = 29,
  tick_size = 1) {

  out <- if (is.null(cols)) data else dplyr::select(data, any_of(cols))
  plot_layout <- ggplot_build(plot)$layout$layout
  n_panel_cols <- max(plot_layout$COL)
  n_panel_rows <- max(plot_layout$ROW)

  plot <- plot +
    ggh4x::force_panelsizes(
      cols = grid::unit(panel_width, "mm"),
      rows = grid::unit(panel_height, "mm"))

  if (is.null(width)) width <- n_panel_cols * (panel_width + tick_size) + 30
  if (is.null(height)) height <- n_panel_rows * (panel_height + tick_size) + 30

  data_path <- publication_output_path(paste0(file_name, "_data.csv"))
  pdf_path <- publication_output_path(paste0(file_name, ".pdf"))
  readr::write_csv(out, data_path)
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "mm")

  invisible(out)
}

#-------------------------------------------------------------------------------
### Panel 1: HF-SK shared mutations ###
#-------------------------------------------------------------------------------
# Mutation-level presence/absence matrix for epithelial tumour compartments.
# Acetone-treated and visually normal mutations are retained in the global
# mutation universe; only the displayed columns are restricted to HF/SK SCC and
# papilloma groups.

tissue_category_groups <- c(
  "HF SCC",
  "HF Papilloma",
  "SK SCC",
  "SK Papilloma")

mutation_base <- tes_mutations_df %>%
  plot_levels() %>%
  filter(
    # !is.na(sample_name), !is.na(SYMBOL),
    !sample_name %in% tes_outliers$sample_name,
    gt_AF >= 0.01,
    !(tissue == "Hair follicle" & as.character(time) == "Week 19"))

shared_mutations_tissue_category <- mutation_base %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("SCC", "Papilloma")) %>%
  mutate(
    tissue_label = recode(
      as.character(tissue),
      "Hair follicle" = "HF",
      "Skin" = "SK"),
    group = paste(tissue_label, category),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
  distinct(mutation_label, group) %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = group,
    values_from = present,
    values_fill = 0L)

shared_mutations_tissue_category[
  setdiff(tissue_category_groups, names(shared_mutations_tissue_category))
] <- 0L

shared_mutations_tissue_category_export <- shared_mutations_tissue_category %>%
  dplyr::select(mutation_label, all_of(tissue_category_groups)) %>%
  as.data.frame()
rownames(shared_mutations_tissue_category_export) <-
  shared_mutations_tissue_category_export$mutation_label

shared_mutations_tissue_category <- shared_mutations_tissue_category %>%
  dplyr::select(all_of(tissue_category_groups)) %>%
  as.data.frame()

HF_SK_shared_mutations <- UpSetR::upset(
  shared_mutations_tissue_category,
  sets = tissue_category_groups,
  nsets = length(tissue_category_groups),
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
  text.scale = c(1.6, 1.6, 1.4, 1.4, 1.6, 1.4),
  point.size = 3,
  line.size = 1)
export_pdf_plot_data(
  data = shared_mutations_tissue_category_export,
  file_name = file.path(figure3_dir, "Figure3_HF_SK_shared_muts"),
  draw = function() print(HF_SK_shared_mutations),
  width = 5,
  height = 5)
HF_SK_shared_mutations

#-------------------------------------------------------------------------------
### Panel 2: Adjacent-distal mutation sharing ###
#-------------------------------------------------------------------------------
# WES-derived mutation sharing is summarised at the strongest HF-SK pair per
# follicle and distance class, then plotted by anatomical adjacency.
wes_hf_skin_max <- wes_hf_skin %>%
  filter(
    hf_treatment == "DT",
    sk_treatment == "DT",
    hf_coverage >= 5) %>%
  group_by(hf_sample_name, distance) %>%
  slice_max(order_by = n_shared, with_ties = FALSE) %>%
  ungroup()

wes_hf_skin_plot_df <- wes_hf_skin_max %>%
  filter(!is.na(adjacency)) %>%
  group_by(hf_sample_name, adjacency) %>%
  slice_max(order_by = n_shared, with_ties = FALSE) %>%
  ungroup()

wes_hf_skin_shared <- ggplot(
  wes_hf_skin_plot_df,
  aes(x = adjacency, y = n_shared, fill = adjacency)) +
  geom_violin(
    colour = "black",
    linewidth = 0.4,
    trim = TRUE) +
  my_theme +
  labs(x = "", y = "# shared mutations") +
  scale_fill_manual(values = col_palette$adjacency_class) +
  theme(legend.position = "") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3)

export_ggplot_data(
  plot = wes_hf_skin_shared,
  data = wes_hf_skin_plot_df,
  file_name = file.path(figure3_dir, "Figure3_wes_adjacency"),
  cols = c(
    "n_shared",
    "distance",
    "adjacency",
    "hf_callable_mbp",
    "hf_coverage",
    "sk_callable_mbp",
    "sk_coverage",
    "hf_category",
    "sk_category"))
wes_hf_skin_shared

# Adjacency statistics from WES/04_results/00_scripts/06_grid_tracking.rmd.
hf_skin_adjacency <- wes_hf_skin %>%
  filter(
    hf_treatment == "DT",
    sk_treatment == "DT",
    hf_coverage >= 5,
    same_mouse,
    adjacency %in% c("Adjacent", "Distal")) %>%
  group_by(hf_sample_name, adjacency) %>%
  slice_max(order_by = n_shared, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    adjacency = factor(adjacency, levels = c("Distal", "Adjacent")),
    pair_callable_mbp = pmin(hf_callable_mbp, sk_callable_mbp),
    sk_cov_z = as.numeric(scale(sk_coverage)))

hf_skin_poisson <- glmmTMB(
  n_shared ~ adjacency +
    (1 | hf_sample_name) +
    offset(log(hf_coverage)) +
    offset(log(sk_coverage)) +
    offset(log(pair_callable_mbp)),
  family = poisson,
  data = hf_skin_adjacency)
summary(hf_skin_poisson)

hf_skin_poisson_null <- glmmTMB(
  n_shared ~ 1 + (1 | hf_sample_name) +
    offset(log(pair_callable_mbp)) +
    offset(log(sk_coverage)) +
    offset(log(hf_coverage)),
  family = poisson,
  data = hf_skin_adjacency)
hf_skin_poisson_lrt <- anova(hf_skin_poisson_null, hf_skin_poisson)
hf_skin_poisson_lrt

hf_skin_poisson_fixed_effects <- as.data.frame(
  summary(hf_skin_poisson)$coefficients$cond)
readr::write_csv(
  data.frame(
    term = rownames(hf_skin_poisson_fixed_effects),
    hf_skin_poisson_fixed_effects,
    row.names = NULL,
    check.names = FALSE),
  file.path(figure3_dir, "Figure3_wes_adjacency_poisson_fixed_effects.csv"))
readr::write_csv(
  data.frame(
    model = rownames(as.data.frame(hf_skin_poisson_lrt)),
    as.data.frame(hf_skin_poisson_lrt),
    row.names = NULL,
    check.names = FALSE),
  file.path(figure3_dir, "Figure3_wes_adjacency_poisson_lrt.csv"))

set.seed(20260504)
n_distal_resample <- 10000L

hf_perm_input <- wes_hf_skin %>%
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
    adjacent_shared = list(n_shared[adjacency == "Adjacent"]),
    distal_shared = list(n_shared[adjacency == "Distal"]),
    .groups = "drop")

hf_distal_wr_input <- hf_perm_input %>%
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

hf_distal_wr_null_median <- apply(distal_wr_max_by_hf, 1, stats::median)
hf_distal_wr_null_q025 <- apply(
  distal_wr_max_by_hf, 1, function(x) as.numeric(stats::quantile(x, 0.025)))
hf_distal_wr_null_q975 <- apply(
  distal_wr_max_by_hf, 1, function(x) as.numeric(stats::quantile(x, 0.975)))
hf_distal_wr_p <- vapply(seq_len(nrow(hf_distal_wr_test_input)), function(i) {
  (sum(distal_wr_max_by_hf[i, ] >=
    hf_distal_wr_test_input$observed_adj_max[[i]]) + 1) /
    (n_distal_resample + 1)
}, numeric(1))

hf_distal_wr_results <- hf_distal_wr_test_input %>%
  mutate(
    distal_wr_null_median = hf_distal_wr_null_median,
    distal_wr_null_q025 = hf_distal_wr_null_q025,
    distal_wr_null_q975 = hf_distal_wr_null_q975,
    effect_vs_null_median = observed_adj_max - distal_wr_null_median,
    p_one_sided = hf_distal_wr_p) %>%
  dplyr::select(
    hf_sample_name,
    n_adj,
    n_distal,
    observed_adj_max,
    distal_wr_null_median,
    distal_wr_null_q025,
    distal_wr_null_q975,
    effect_vs_null_median,
    p_one_sided)

global_distal_wr_summary
hf_distal_wr_results

readr::write_csv(
  global_distal_wr_summary,
  file.path(figure3_dir,
    "Figure3_wes_adjacency_replacement_distal_null_global.csv"))
readr::write_csv(
  hf_distal_wr_results,
  file.path(figure3_dir,
    "Figure3_wes_adjacency_replacement_distal_null_by_hf.csv"))

#-------------------------------------------------------------------------------
### Panel 3: Mutations over time and space ###
#-------------------------------------------------------------------------------

timepoint_levels <- c("Week 8", "Week 14", "Week 17", "Skin")
time_space_plot_timepoints <- c("Week 8", "Week 14", "Week 17")
sharing_class_levels <- c(
  "Shared with SK",
  "Shared with later HF",
  "Shared with earlier HF")

mutation_timepoint_presence <- mutation_base %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "SCC")) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
  filter(!is.na(plot_timepoint)) %>%
  filter(!(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
  distinct(mutation_label, plot_timepoint, tissue, category) %>%
  mutate(
    plot_timepoint = factor(plot_timepoint, levels = timepoint_levels),
    timepoint_rank = as.integer(plot_timepoint))

weighted_se <- function(x, w) {
  valid <- !is.na(x) & !is.na(w) & w > 0
  x <- x[valid]
  w <- w[valid]

  if (length(x) < 2) {
    return(NA_real_)
  }

  w_sum <- sum(w)
  w_mean <- weighted.mean(x, w)
  w_var <- sum(w * (x - w_mean)^2) / w_sum
  n_eff <- w_sum^2 / sum(w^2)

  if (n_eff <= 1) {
    return(NA_real_)
  }

  sqrt(w_var / n_eff)
}

weighted_sd <- function(x, w) {
  valid <- !is.na(x) & !is.na(w) & w > 0
  x <- x[valid]
  w <- w[valid]

  if (length(x) < 2) {
    return(NA_real_)
  }

  w_mean <- weighted.mean(x, w)
  sqrt(sum(w * (x - w_mean)^2) / sum(w))
}

mutation_sample_timepoint_flags <- mutation_base %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("Papilloma", "SCC")) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
  filter(!is.na(plot_timepoint)) %>%
  filter(!(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
  distinct(sample_name, mutation_label, plot_timepoint) %>%
  mutate(
    plot_timepoint = factor(plot_timepoint, levels = timepoint_levels),
    timepoint_rank = as.integer(plot_timepoint)) %>%
  filter(plot_timepoint %in% time_space_plot_timepoints) %>%
  mutate(
    shared_skin = mutation_label %in%
      mutation_timepoint_presence$mutation_label[
        mutation_timepoint_presence$plot_timepoint == "Skin" &
          mutation_timepoint_presence$category %in% c("SCC", "Papilloma")]) %>%
  rowwise() %>%
  mutate(
    shared_later_hf = any(
      mutation_timepoint_presence$mutation_label == mutation_label &
        mutation_timepoint_presence$tissue == "Hair follicle" &
        mutation_timepoint_presence$timepoint_rank > timepoint_rank),
    shared_earlier_hf = any(
      mutation_timepoint_presence$mutation_label == mutation_label &
        mutation_timepoint_presence$tissue == "Hair follicle" &
        mutation_timepoint_presence$timepoint_rank < timepoint_rank)) %>%
  ungroup()

mutation_sample_timepoint_totals <- mutation_sample_timepoint_flags %>%
  count(plot_timepoint, sample_name, name = "total_mutations")

mutation_sample_timepoint_long <- mutation_sample_timepoint_flags %>%
  pivot_longer(
    cols = c(shared_skin, shared_later_hf, shared_earlier_hf),
    names_to = "sharing_class",
    values_to = "included") %>%
  filter(included) %>%
  mutate(
    sharing_class = recode(
      sharing_class,
      shared_skin = "Shared with SK",
      shared_later_hf = "Shared with later HF",
      shared_earlier_hf = "Shared with earlier HF"),
    sharing_class = factor(sharing_class, levels = sharing_class_levels))

mutation_sample_timepoint_pct <- mutation_sample_timepoint_totals %>%
  dplyr::select(plot_timepoint, sample_name, total_mutations) %>%
  mutate(dummy_join = 1L) %>%
  left_join(
    data.frame(
      sharing_class = factor(sharing_class_levels, levels = sharing_class_levels),
      dummy_join = 1L),
    by = "dummy_join",
    relationship = "many-to-many") %>%
  dplyr::select(-dummy_join) %>%
  left_join(
    mutation_sample_timepoint_long %>%
      count(
        plot_timepoint,
        sample_name,
        sharing_class,
        name = "n_mutations") %>%
      mutate(sharing_class = factor(sharing_class, levels = sharing_class_levels)),
    by = c("plot_timepoint", "sample_name", "sharing_class")) %>%
  mutate(
    n_mutations = coalesce(n_mutations, 0L),
    pct_mutations = n_mutations / total_mutations)

mutations_time_space_bb_counts <- mutation_sample_timepoint_pct %>%
  filter(sharing_class == "Shared with SK") %>%
  left_join(
    tes_metadata %>%
      dplyr::select(sample_name, callable_mbp, coverage),
    by = "sample_name") %>%
  mutate(
    time_num = case_when(
      plot_timepoint == "Week 8" ~ 8,
      plot_timepoint == "Week 14" ~ 14,
      plot_timepoint == "Week 17" ~ 17),
    n_shared_skin = n_mutations,
    other_shared_mutations = total_mutations - n_shared_skin,
    log_coverage = log(coverage),
    shared_skin_pct = n_shared_skin / total_mutations) %>%
  filter(
    !is.na(time_num),
    total_mutations > 0,
    is.finite(log_coverage),
    is.finite(callable_mbp))

mutations_time_space_bb_model <- glmmTMB(
  cbind(n_shared_skin, other_shared_mutations) ~
    time_num + log_coverage,
  family = binomial(link = "logit"),
  data = mutations_time_space_bb_counts)
summary(mutations_time_space_bb_model)
mutations_time_space_bb_emmeans <- emmeans(
  mutations_time_space_bb_model,
  ~ time_num,
  at = list(time_num = c(8, 14, 17)),
  type = "response")
summary(mutations_time_space_bb_emmeans)
mutations_time_space_bb_trend <- emtrends(
  mutations_time_space_bb_model,
  specs = ~ 1,
  var = "time_num")
summary(mutations_time_space_bb_trend)

mutations_time_space_weighted_df <- mutation_sample_timepoint_pct %>%
  group_by(plot_timepoint, sharing_class) %>%
  summarise(
    n_samples = n_distinct(sample_name),
    total_mutations_weight = sum(total_mutations),
    n_mutations = sum(n_mutations),
    weighted_mean_pct = weighted.mean(pct_mutations, total_mutations),
    weighted_sd_pct = weighted_sd(pct_mutations, total_mutations),
    weighted_se_pct = weighted_se(pct_mutations, total_mutations),
    ymin = pmax(weighted_mean_pct - weighted_se_pct, 0),
    ymax = pmin(weighted_mean_pct + weighted_se_pct, 1),
    .groups = "drop")

mutations_over_time_space_weighted <- ggplot(
  mutations_time_space_weighted_df,
  aes(x = plot_timepoint, y = weighted_mean_pct, fill = sharing_class)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "black") +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    position = position_dodge(width = 0.8),
    width = 0.2) +
  scale_fill_manual(
    values = c(
      "Shared with SK" = "#f0b981",
      "Shared with later HF" = "#300358",
      "Shared with earlier HF" = "#9780aa"),
    drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "", y = "Weighted mean % of mutations", fill = "") +
  my_theme +
  theme(legend.position = "bottom")

readr::write_csv(
  mutations_time_space_bb_counts %>%
    dplyr::select(any_of(c(
    "sample_name",
    "plot_timepoint",
    "time_num",
    "n_shared_skin",
    "other_shared_mutations",
    "total_mutations",
    "shared_skin_pct",
    "callable_mbp",
    "coverage",
    "log_coverage"))),
  file.path(figure3_dir, "Figure3_mutations_over_time_space_bb_counts.csv"))
readr::write_csv(
  as_tibble(as.data.frame(summary(mutations_time_space_bb_emmeans))),
  file.path(figure3_dir, "Figure3_mutations_over_time_space_bb_emmeans.csv"))
readr::write_csv(
  as_tibble(as.data.frame(summary(mutations_time_space_bb_trend, infer = TRUE))),
  file.path(figure3_dir, "Figure3_mutations_over_time_space_bb_trend.csv"))
export_ggplot_data(
  plot = mutations_over_time_space_weighted,
  data = mutations_time_space_weighted_df,
  file_name = file.path(figure3_dir,
    "Figure3_mutations_over_time_space_weighted"),
  cols = c(
    "plot_timepoint",
    "sharing_class",
    "n_samples",
    "total_mutations_weight",
    "n_mutations",
    "weighted_mean_pct",
    "weighted_sd_pct",
    "weighted_se_pct",
    "ymin",
    "ymax"))
mutations_over_time_space_weighted

# Restore the submitted Figure 4F calculation without changing the current
# non-exclusive weighted summary used by supplementary analyses.
source(file.path(figure3_source_dir, "Figure3_mutations_over_time_space.R"))

#-------------------------------------------------------------------------------
### Panel 3b: Early mutations retained in SCC-associated regions ###
#-------------------------------------------------------------------------------
# A region is defined by mouse and grid and is considered SCC-associated when
# any DT-treated HF or skin sample from that region is classified as SCC. This
# retrospective definition retains regions that were morphologically normal at
# Week 8/14 but subsequently developed SCC. Exact alleles are tracked only
# within the same mouse and grid, preventing cross-animal or cross-region
# matches. Week 19 HF samples are excluded, consistently with the other panels;
# the Week 19 skin sample is retained.

early_scc_timepoint_levels <- c("Week 8", "Week 14")
early_scc_detection_levels <- c(
  "Later HF and Skin",
  "Later HF only",
  "Skin only",
  "Not detected later")

scc_region_definitions <- tes_metadata %>%
  plot_levels() %>%
  filter(
    treatment == "DT",
    category == "SCC",
    tissue %in% c("Hair follicle", "Skin"),
    !is.na(mouse),
    !is.na(grid)) %>%
  group_by(mouse, grid) %>%
  summarise(
    scc_sample_names = paste(sort(unique(sample_name)), collapse = "; "),
    scc_tissues = paste(sort(unique(as.character(tissue))), collapse = "; "),
    scc_timepoints = paste(sort(unique(as.character(time))), collapse = "; "),
    .groups = "drop")

scc_region_tracking_base <- tes_mutations_df %>%
  plot_levels() %>%
  mutate(
    gt_AF = suppressWarnings(as.numeric(gt_AF)),
    time_character = as.character(time),
    tracking_timepoint = if_else(
      tissue == "Skin",
      "Skin",
      time_character),
    time_num = case_when(
      tissue == "Skin" ~ 19L,
      tissue == "Hair follicle" ~ suppressWarnings(
        as.integer(str_extract(time_character, "[0-9]+")))),
    compartment = recode(
      as.character(tissue),
      "Hair follicle" = "HF",
      "Skin" = "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    !sample_name %in% tes_outliers$sample_name,
    !is.na(sample_name),
    !is.na(mouse),
    !is.na(grid),
    !is.na(time_num),
    gt_AF >= 0.01,
    !(tissue == "Hair follicle" & time_character == "Week 19")) %>%
  semi_join(
    scc_region_definitions,
    by = c("mouse", "grid")) %>%
  distinct(sample_name, mutation_id, .keep_all = TRUE)

early_scc_region_mutations <- scc_region_tracking_base %>%
  filter(
    tissue == "Hair follicle",
    tracking_timepoint %in% early_scc_timepoint_levels) %>%
  transmute(
    mouse,
    grid,
    grid_x,
    grid_y,
    early_sample_name = sample_name,
    early_timepoint = factor(
      tracking_timepoint,
      levels = early_scc_timepoint_levels),
    early_time_num = time_num,
    early_category = as.character(category),
    mutation_id,
    CHROM,
    POS,
    REF,
    ALT,
    mutation,
    SYMBOL,
    Consequence,
    VARIANT_CLASS,
    early_gt_AF = gt_AF,
    early_gt_DP = gt_DP,
    early_coverage = coverage,
    early_callable_mbp = callable_mbp)

early_scc_region_later_hits <- early_scc_region_mutations %>%
  dplyr::select(
    mouse,
    grid,
    early_sample_name,
    early_timepoint,
    early_time_num,
    mutation_id) %>%
  inner_join(
    scc_region_tracking_base %>%
      transmute(
        mouse,
        grid,
        mutation_id,
        later_sample_name = sample_name,
        later_timepoint = tracking_timepoint,
        later_time_num = time_num,
        later_compartment = compartment,
        later_category = as.character(category),
        later_gt_AF = gt_AF,
        later_gt_DP = gt_DP,
        later_coverage = coverage,
        later_callable_mbp = callable_mbp),
    by = c("mouse", "grid", "mutation_id"),
    relationship = "many-to-many") %>%
  filter(later_time_num > early_time_num) %>%
  arrange(
    early_time_num,
    mouse,
    grid,
    mutation_id,
    later_time_num,
    later_compartment,
    later_sample_name)

early_scc_region_later_flags <- early_scc_region_later_hits %>%
  group_by(
    mouse,
    grid,
    early_sample_name,
    early_timepoint,
    early_time_num,
    mutation_id) %>%
  summarise(
    detected_later_hf = any(later_compartment == "HF"),
    detected_later_skin = any(later_compartment == "Skin"),
    n_later_samples = n_distinct(later_sample_name),
    later_sample_names = paste(
      sort(unique(later_sample_name)),
      collapse = "; "),
    later_timepoints = paste(
      sort(unique(later_timepoint)),
      collapse = "; "),
    later_categories = paste(
      sort(unique(later_category)),
      collapse = "; "),
    max_later_gt_AF = max(later_gt_AF, na.rm = TRUE),
    .groups = "drop")

early_scc_region_tracking <- early_scc_region_mutations %>%
  left_join(
    early_scc_region_later_flags,
    by = c(
      "mouse",
      "grid",
      "early_sample_name",
      "early_timepoint",
      "early_time_num",
      "mutation_id")) %>%
  mutate(
    detected_later_hf = coalesce(detected_later_hf, FALSE),
    detected_later_skin = coalesce(detected_later_skin, FALSE),
    detected_later = detected_later_hf | detected_later_skin,
    n_later_samples = coalesce(n_later_samples, 0L),
    detection_outcome = case_when(
      detected_later_hf & detected_later_skin ~ "Later HF and Skin",
      detected_later_hf ~ "Later HF only",
      detected_later_skin ~ "Skin only",
      TRUE ~ "Not detected later"),
    detection_outcome = factor(
      detection_outcome,
      levels = early_scc_detection_levels)) %>%
  arrange(early_time_num, mouse, grid, mutation_id)

early_scc_region_tracking_summary <- early_scc_region_tracking %>%
  count(early_timepoint, detection_outcome, name = "n_mutations") %>%
  complete(
    early_timepoint = factor(
      early_scc_timepoint_levels,
      levels = early_scc_timepoint_levels),
    detection_outcome = factor(
      early_scc_detection_levels,
      levels = early_scc_detection_levels),
    fill = list(n_mutations = 0L)) %>%
  group_by(early_timepoint) %>%
  mutate(
    total_early_mutations = sum(n_mutations),
    pct_early_mutations = n_mutations / total_early_mutations) %>%
  ungroup()

early_scc_region_sample_totals <- early_scc_region_tracking %>%
  count(
    mouse,
    grid,
    early_sample_name,
    early_timepoint,
    name = "total_early_mutations")

early_scc_region_sample_summary <- early_scc_region_tracking %>%
  dplyr::select(
    mouse,
    grid,
    early_sample_name,
    early_timepoint,
    mutation_id,
    detected_later_hf,
    detected_later_skin) %>%
  pivot_longer(
    cols = c(detected_later_hf, detected_later_skin),
    names_to = "later_compartment",
    values_to = "detected") %>%
  mutate(
    later_compartment = recode(
      later_compartment,
      detected_later_hf = "Later HF",
      detected_later_skin = "Skin"),
    later_compartment = factor(
      later_compartment,
      levels = c("Later HF", "Skin"))) %>%
  group_by(
    mouse,
    grid,
    early_sample_name,
    early_timepoint,
    later_compartment) %>%
  summarise(
    n_detected_later = sum(detected),
    .groups = "drop") %>%
  left_join(
    early_scc_region_sample_totals,
    by = c("mouse", "grid", "early_sample_name", "early_timepoint")) %>%
  mutate(pct_detected_later = n_detected_later / total_early_mutations)

early_scc_region_tracked_genes <- early_scc_region_tracking %>%
  filter(
    detected_later,
    !is.na(SYMBOL),
    SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(
    n_early_observations = n(),
    n_mutations = n_distinct(mutation_id),
    n_regions = n_distinct(paste(mouse, grid, sep = ":")),
    n_early_samples = n_distinct(early_sample_name),
    later_hf_only = sum(detected_later_hf & !detected_later_skin),
    skin_only = sum(!detected_later_hf & detected_later_skin),
    later_hf_and_skin = sum(detected_later_hf & detected_later_skin),
    n_later_hf_observations = sum(detected_later_hf),
    n_skin_observations = sum(detected_later_skin),
    detected_later_hf = any(detected_later_hf),
    detected_later_skin = any(detected_later_skin),
    .groups = "drop") %>%
  arrange(
    dplyr::desc(n_early_observations),
    dplyr::desc(n_mutations),
    SYMBOL)

early_scc_gene_order <- early_scc_region_tracked_genes$SYMBOL

early_scc_region_tracked_gene_plot_df <-
  early_scc_region_tracked_genes %>%
  dplyr::select(
    SYMBOL,
    later_hf_only,
    skin_only,
    later_hf_and_skin) %>%
  pivot_longer(
    cols = c(later_hf_only, skin_only, later_hf_and_skin),
    names_to = "detection_class",
    values_to = "n") %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = rev(early_scc_gene_order)),
    plot_panel = factor("Tracked early observations"),
    detection_class = factor(
      recode(
        detection_class,
        later_hf_only = "Later HF only",
        skin_only = "Skin only",
        later_hf_and_skin = "Later HF and Skin"),
      levels = c(
        "Later HF only",
        "Skin only",
        "Later HF and Skin")))

stopifnot(
  nrow(scc_region_definitions) > 0,
  nrow(early_scc_region_mutations) > 0,
  all(
    early_scc_region_later_hits$later_time_num >
      early_scc_region_later_hits$early_time_num),
  all(
    (early_scc_region_tracking$n_later_samples == 0L) ==
      !early_scc_region_tracking$detected_later),
  sum(early_scc_region_tracked_genes$n_early_observations) ==
    sum(
      early_scc_region_tracked_genes$later_hf_only +
        early_scc_region_tracked_genes$skin_only +
        early_scc_region_tracked_genes$later_hf_and_skin))

early_scc_mutations_later_detection <- ggplot(
  early_scc_region_sample_summary,
  aes(x = early_timepoint, y = pct_detected_later)) +
  ggbeeswarm::geom_quasirandom(
    width = 0.18,
    size = 1.8,
    alpha = 0.8) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.45,
    colour = "#D55E00") +
  facet_wrap(~ later_compartment, nrow = 1) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Early HF timepoint",
    y = "% early mutations detected later") +
  my_theme +
  theme(legend.position = "none")

early_scc_tracked_gene_summary_plot <- ggplot(
  early_scc_region_tracked_gene_plot_df,
  aes(x = n, y = SYMBOL, fill = detection_class)) +
  geom_col(
    width = 0.75,
    colour = "black",
    linewidth = 0.25) +
  facet_grid(. ~ plot_panel) +
  scale_fill_manual(
    values = c(
      "Later HF only" = "#300358",
      "Skin only" = "#f0b981",
      "Later HF and Skin" = "#008B8B"),
    drop = FALSE) +
  scale_x_continuous(
    breaks = scales::breaks_width(5),
    expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "# mutations",
    y = "",
    fill = "") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  my_theme +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8))

readr::write_csv(
  scc_region_definitions,
  file.path(figure3_dir, "Figure3_early_SCC_region_definitions.csv"))
readr::write_csv(
  early_scc_region_later_hits,
  file.path(figure3_dir, "Figure3_early_SCC_region_tracking_hits.csv"))
readr::write_csv(
  early_scc_region_tracking,
  file.path(figure3_dir, "Figure3_early_SCC_region_tracking_mutations.csv"))
readr::write_csv(
  early_scc_region_tracking_summary,
  file.path(figure3_dir, "Figure3_early_SCC_region_tracking_summary.csv"))
readr::write_csv(
  early_scc_region_tracked_genes,
  file.path(figure3_dir, "Figure3_early_SCC_region_tracking_genes.csv"))
export_ggplot_data(
  plot = early_scc_mutations_later_detection,
  data = early_scc_region_sample_summary,
  file_name = file.path(figure3_dir,
    "Figure3_early_SCC_region_mutations_later_detection"),
  cols = c(
    "mouse",
    "grid",
    "early_sample_name",
    "early_timepoint",
    "later_compartment",
    "n_detected_later",
    "total_early_mutations",
    "pct_detected_later"))
early_scc_mutations_later_detection

export_ggplot_data(
  plot = early_scc_tracked_gene_summary_plot,
  data = early_scc_region_tracked_gene_plot_df,
  file_name = file.path(figure3_dir,
    "Figure3_early_SCC_region_tracked_gene_summary"),
  cols = c(
    "SYMBOL",
    "plot_panel",
    "detection_class",
    "n"),
  width = 100,
  height = 185,
  panel_width = 55,
  panel_height = 145)
early_scc_tracked_gene_summary_plot

#-------------------------------------------------------------------------------
### Panel 4: Recurrent mutation oncoplot ###
#-------------------------------------------------------------------------------
# Recurrent mutations are defined as variants observed in at least two HF/SK
# timepoints. Genes are ordered by recurrence, broad functional grouping, and
# first detection, with canonical cSCC genes highlighted in red.

morphology_levels <- c("SCC", "Papilloma")
mutation_nature_levels <- c(
  "Missense",
  "Insertion",
  "Deletion",
  "DBS")
mutation_effect_levels <- c("Nonsynonymous", "Synonymous")
plot_text_size <- 8
plot_text_size_mm <- plot_text_size / ggplot2::.pt
gene_family_levels <- c(
  "Keratin-associated proteins",
  "Genome integrity",
  "Tumour suppressor genes",
  "Epigenetic regulators",
  "Cell adhesion and cytoskeleton",
  "Ras",
  "Metabolism",
  "Other")
filter_panel4_human_orthologues <- TRUE
epigenetic_plot_order <- c(
  "Kmt2d",
  "Kdm6a",
  "Kmt2c",
  "Atrx",
  "Kdm5b",
  "Ash1l",
  "Cmtr1",
  "Smarca4",
  "Setd2")
KNOWN_CSCC_GENES <- c(
  "Ajuba", "Arid2", "Asxl1", "Casp8", "Card11", "Ccnd1", "Cdkn2a", "Chuk",
  "Crebbp", "Egfr", "Ep300", "Erbb2", "Erbb3", "Erbb4", "Ezh2", "Fat1",
  "Hras", "Iqgap1", "Irf6", "Kmt2a", "Kmt2c", "Kmt2d", "Knstrn", "Kras",
  "Mtor", "Myh9", "Ncor1", "Nfe2l2", "Notch1", "Notch2", "Notch3", "Nrp1",
  "Pbrm1", "Pik3ca", "Ptch1", "Pten", "Rhbdf2", "Ros1", "Rras2", "Runx1",
  "Setd2", "Tert", "Tp53", "Tp63", "Tsc1", "Usp28", "Vegfa")

plot_class <- function(consequence, variant_class, ref, alt) {
  consequence <- coalesce(as.character(consequence), "")
  variant_class <- coalesce(as.character(variant_class), "")
  ref <- coalesce(as.character(ref), "")
  alt <- coalesce(as.character(alt), "")

  case_when(
    str_detect(consequence, "stop_gained|stop_lost") ~ "Stop gain/loss",
    variant_class == "insertion" | nchar(ref) < nchar(alt) ~ "Insertion",
    variant_class == "deletion" | nchar(ref) > nchar(alt) ~ "Deletion",
    variant_class %in% c("substitution", "MNV", "DNV") |
      (nchar(ref) == 2 & nchar(alt) == 2) ~ "DBS",
    str_detect(consequence, "missense_variant") ~ "Missense",
    str_detect(consequence, "synonymous_variant") ~ "Synonymous",
    TRUE ~ "Other")
}

plot_class_with_intronic_utr <- function(consequence, variant_class, ref, alt) {
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
    TRUE ~ "Other")
}

plot_effect <- function(consequence) {
  if_else(
    str_detect(coalesce(as.character(consequence), ""), "synonymous_variant"),
    "Synonymous",
    "Nonsynonymous")
}

plot_effect_with_intronic_utr <- function(consequence) {
  consequence <- coalesce(as.character(consequence), "")

  case_when(
    str_detect(consequence, "intron_variant|UTR_variant") ~ "Gene-flanking",
    str_detect(consequence, "synonymous_variant") ~ "Synonymous",
    TRUE ~ "Nonsynonymous")
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

interval_union_width <- function(start, end) {
  intervals <- tibble(start = start, end = end) %>%
    filter(!is.na(start), !is.na(end)) %>%
    arrange(start, end)

  if (nrow(intervals) == 0) {
    return(NA_real_)
  }

  total_width <- 0
  current_start <- intervals$start[[1]]
  current_end <- intervals$end[[1]]

  for (i in seq_len(nrow(intervals))[-1]) {
    if (intervals$start[[i]] <= current_end + 1) {
      current_end <- max(current_end, intervals$end[[i]])
    } else {
      total_width <- total_width + current_end - current_start + 1
      current_start <- intervals$start[[i]]
      current_end <- intervals$end[[i]]
    }
  }

  total_width + current_end - current_start + 1
}

biopsy_group_label <- function(plot_timepoint, morphology) {
  morphology <- if_else(
    morphology == "Morphologically normal",
    "Morphologically\nnormal",
    morphology)

  if_else(
    plot_timepoint == "Skin",
    paste("Skin", morphology, sep = "\n"),
    paste(plot_timepoint, morphology, sep = "\n"))
}

biopsy_group_levels <- unlist(lapply(
  morphology_levels,
  function(morphology) {
    biopsy_group_label(timepoint_levels, morphology)
  }))

biopsy_order <- tes_metadata %>%
  plot_levels() %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("SCC", "Papilloma"),
    !sample_name %in% tes_outliers$sample_name) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    morphology = recode(
      as.character(category),
      "Visually normal" = "Morphologically normal"),
    biopsy_group = biopsy_group_label(plot_timepoint, morphology),
    plot_timepoint = factor(plot_timepoint, levels = timepoint_levels),
    morphology = factor(morphology, levels = morphology_levels),
    biopsy_group = factor(biopsy_group, levels = biopsy_group_levels)) %>%
  filter(!is.na(plot_timepoint), !is.na(morphology)) %>%
  arrange(morphology, plot_timepoint, sample_name) %>%
  mutate(biopsy_rank = row_number()) %>%
  dplyr::select(
    sample_name,
    biopsy_rank,
    plot_timepoint,
    morphology,
    biopsy_group)
biopsy_order_all <- biopsy_order

if (!exists("gr_genes")) {
  load(refcds_path)
}

gene_lengths_kb <- as.data.frame(gr_genes) %>%
  as_tibble() %>%
  group_by(gene_symbol = names) %>%
  summarise(
    gene_length_kb = interval_union_width(start, end) / 1e3,
    .groups = "drop")

if (identical(EXTERNAL_DATA_MODE, "frozen")) {
  frozen_orthologue_sources <- c(
    publication_repo_input(
      "Figures/Figure 3/Figure3_oncoplot_msigdb_annotations_data.csv",
      "a564167836a9228b6f7d7795e62326e1dbcdc566d24777b1ec32137783c4df7b"),
    publication_repo_input(
      paste0("Figures/Figure 3/",
        "Figure3_oncoplot_msigdb_annotations_synonymous_data.csv"),
      "56b657035d11674014e6114bfdf58e046056e89039f8516492799e624b7ef669"))
  mouse_human_orthologues <- frozen_orthologue_sources %>%
    lapply(readr::read_csv, show_col_types = FALSE) %>%
    bind_rows() %>%
    select(symbol = gene_symbol, human_symbols) %>%
    tidyr::separate_rows(human_symbols, sep = ";") %>%
    transmute(symbol, human_symbol = trimws(human_symbols)) %>%
    filter(!is.na(human_symbol), human_symbol != "") %>%
    distinct()
} else {
  publication_require_live_refresh()
  mouse_human_orthologues <- babelgene::orthologs(
    genes = sort(unique(na.omit(tes_mutations_df$SYMBOL))),
    species = "mouse",
    human = FALSE) %>%
    filter(!is.na(human_symbol), human_symbol != "") %>%
    distinct(symbol, human_symbol)
}

gene_family_assignments <- mouse_human_orthologues %>%
  group_by(symbol) %>%
  summarise(
    human_symbols = paste(sort(unique(human_symbol)), collapse = ";"),
    .groups = "drop") %>%
  mutate(
    gene_symbol = symbol,
    gene_family = gene_family_label(gene_symbol, human_symbols)) %>%
  dplyr::select(gene_symbol, human_symbols, gene_family)

filter_human_orthologues <- function(data) {
  if (!filter_panel4_human_orthologues) {
    return(data)
  }

  semi_join(data, mouse_human_orthologues, by = c("SYMBOL" = "symbol"))
}

mutation_biopsy_base <- tes_mutations_df %>%
  plot_levels() %>%
  mutate(gt_AF = suppressWarnings(as.numeric(gt_AF))) %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("SCC", "Papilloma"),
    gt_AF >= 0.01,
    Consequence != "intergenic_variant",
    Feature_type != "RegulatoryFeature",
    !is.na(sample_name),
    !is.na(SYMBOL),
    SYMBOL != "",
    !sample_name %in% tes_outliers$sample_name) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_nature = plot_class(Consequence, VARIANT_CLASS, REF, ALT),
    mutation_effect = plot_effect(Consequence)) %>%
  filter(!is.na(plot_timepoint)) %>%
  semi_join(biopsy_order, by = "sample_name") %>%
  filter_human_orthologues() %>%
  distinct(
    sample_name,
    plot_timepoint,
    mutation_id,
    SYMBOL,
    mutation_nature,
    mutation_effect,
    .keep_all = TRUE)

recurrent_mutation_ids <- mutation_biopsy_base %>%
  distinct(mutation_id, plot_timepoint) %>%
  group_by(mutation_id) %>%
  summarise(n_timepoints = n_distinct(plot_timepoint), .groups = "drop") %>%
  filter(n_timepoints >= 2)

recurrent_mutations <- mutation_biopsy_base %>%
  semi_join(recurrent_mutation_ids, by = "mutation_id")

biopsy_order <- biopsy_order %>%
  semi_join(
    recurrent_mutations %>%
      filter(mutation_nature %in% mutation_nature_levels) %>%
      distinct(sample_name),
    by = "sample_name") %>%
  arrange(morphology, plot_timepoint, sample_name) %>%
  mutate(biopsy_rank = row_number())

gene_order <- recurrent_mutations %>%
  filter(mutation_nature %in% mutation_nature_levels) %>%
  group_by(gene_symbol = SYMBOL) %>%
  summarise(
    n_mutated_biopsies = n_distinct(sample_name),
    n_recurrent_mutations = n_distinct(mutation_id),
    n_nonsynonymous_mutations = n_distinct(
      mutation_id[mutation_effect == "Nonsynonymous"]),
    first_biopsy_rank = min(
      biopsy_order$biopsy_rank[
        match(sample_name, biopsy_order$sample_name)
      ],
      na.rm = TRUE),
    .groups = "drop") %>%
  left_join(gene_family_assignments, by = "gene_symbol") %>%
  mutate(
    gene_family = coalesce(gene_family, "Other"),
    gene_family = factor(gene_family, levels = gene_family_levels)) %>%
  left_join(gene_lengths_kb, by = "gene_symbol") %>%
  mutate(
    epigenetic_plot_rank = coalesce(
      match(gene_symbol, epigenetic_plot_order),
      Inf)) %>%
  group_by(gene_family) %>%
  mutate(
    family_mutated_biopsies = sum(n_mutated_biopsies),
    family_recurrent_mutations = sum(n_recurrent_mutations),
    family_n_genes = n()) %>%
  ungroup() %>%
  arrange(
    gene_family == "Other",
    dplyr::desc(family_mutated_biopsies),
    dplyr::desc(family_recurrent_mutations),
    gene_family,
    epigenetic_plot_rank,
    dplyr::desc(n_mutated_biopsies),
    dplyr::desc(n_recurrent_mutations),
    dplyr::desc(n_nonsynonymous_mutations),
    first_biopsy_rank,
    gene_symbol) %>%
  mutate(
    gene_rank = row_number(),
    family_rank = as.integer(factor(gene_family, levels = unique(gene_family))),
    n_gene_families = n_distinct(gene_family),
    plot_y = n() - gene_rank + 1 + (n_gene_families - family_rank) * 0.35,
    known_cscc_gene = gene_symbol %in% KNOWN_CSCC_GENES)

# MSigDB annotations for manual review of oncoplot functional grouping.
# The table keeps the oncoplot mouse symbols and preserves the original
# MSigDB/source symbols. Set MSIGDB_DB_SPECIES to "HS" to use human MSigDB with
# mouse ortholog mapping, or "MM" to use mouse-native MSigDB gene sets.
MSIGDB_DB_SPECIES <- "MM"

collapse_unique <- function(x) {
  x <- sort(unique(na.omit(as.character(x))))
  x <- x[x != ""]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

msigdb_collection_code <- function(collection_name, db_species = MSIGDB_DB_SPECIES) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "The msigdbr package is required to build the oncoplot gene annotation table. ",
      "Install it with install.packages('msigdbr').",
      call. = FALSE)
  }

  collection_map <- list(
    HS = c(hallmark = "H", pathways = "C2", go = "C5"),
    MM = c(hallmark = "MH", pathways = "M2", go = "M5"))
  db_species <- toupper(db_species)

  if (!db_species %in% names(collection_map)) {
    stop("Unsupported MSigDB db_species: ", db_species, call. = FALSE)
  }

  collection <- unname(collection_map[[db_species]][[collection_name]])
  available_collections <- msigdbr::msigdbr_collections(
    db_species = db_species)$gs_collection

  if (is.na(collection) || !collection %in% available_collections) {
    stop(
      "MSigDB collection '",
      collection_name,
      "' is unavailable for db_species = '",
      db_species,
      "'. Available collections: ",
      paste(sort(unique(available_collections)), collapse = ", "),
      call. = FALSE)
  }

  collection
}

load_msigdbr_collection <- function(collection_name, db_species = MSIGDB_DB_SPECIES, species = "Mus musculus") {
  msigdbr_args <- list(species = species)
  msigdbr_formals <- names(formals(msigdbr::msigdbr))
  collection <- msigdb_collection_code(collection_name, db_species)

  if ("db_species" %in% msigdbr_formals) {
    msigdbr_args$db_species <- db_species
  }
  if ("collection" %in% msigdbr_formals) {
    msigdbr_args$collection <- collection
  } else {
    msigdbr_args$category <- collection
  }

  gene_sets <- do.call(msigdbr::msigdbr, msigdbr_args)

  if (!"gs_collection" %in% names(gene_sets) && "gs_cat" %in% names(gene_sets)) {
    gene_sets <- gene_sets %>% dplyr::rename(gs_collection = gs_cat)
  }
  if (!"gs_subcollection" %in% names(gene_sets) && "gs_subcat" %in% names(gene_sets)) {
    gene_sets <- gene_sets %>% dplyr::rename(gs_subcollection = gs_subcat)
  }
  if (!"db_version" %in% names(gene_sets) && "msigdb_version" %in% names(gene_sets)) {
    gene_sets <- gene_sets %>% dplyr::rename(db_version = msigdb_version)
  }
  if (!"gs_collection" %in% names(gene_sets)) {
    gene_sets$gs_collection <- collection
  }
  if (!"gs_subcollection" %in% names(gene_sets)) {
    gene_sets$gs_subcollection <- NA_character_
  }

  gene_sets
}

standardise_msigdb_annotation_cols <- function(data) {
  required_cols <- c(
    "gene_symbol",
    "ncbi_gene",
    "ensembl_gene",
    "db_gene_symbol",
    "db_ncbi_gene",
    "db_ensembl_gene",
    "source_gene",
    "gs_id",
    "gs_name",
    "gs_collection",
    "gs_subcollection",
    "gs_collection_name",
    "gs_description",
    "gs_source_species",
    "gs_pmid",
    "gs_geoid",
    "gs_exact_source",
    "gs_url",
    "db_version",
    "db_target_species",
    "ortholog_sources",
    "num_ortholog_sources")

  missing_cols <- setdiff(required_cols, names(data))
  data[missing_cols] <- NA
  data %>% dplyr::select(all_of(required_cols))
}

build_oncoplot_msigdb_annotations <- function(gene_tbl) {
  gene_tbl_base <- gene_tbl %>%
    dplyr::select(
      gene_symbol,
      gene_rank,
      gene_family,
      human_symbols,
      gene_length_kb,
      n_mutated_biopsies,
      n_recurrent_mutations,
      n_nonsynonymous_mutations,
      known_cscc_gene)

  hallmark_sets <- load_msigdbr_collection("hallmark") %>%
    mutate(
      annotation_database = "Hallmark",
      annotation_source = "Hallmark")

  pathway_sets <- load_msigdbr_collection("pathways") %>%
    mutate(
      annotation_database = case_when(
        str_detect(gs_subcollection, "REACTOME") ~ "Reactome",
        str_detect(gs_subcollection, "KEGG") ~ "KEGG",
        str_detect(gs_subcollection, "WIKIPATHWAYS") ~ "WikiPathways",
        TRUE ~ NA_character_),
      annotation_source = case_when(
        str_detect(gs_subcollection, "REACTOME") ~ "Reactome",
        str_detect(gs_subcollection, "KEGG") ~ "KEGG",
        str_detect(gs_subcollection, "WIKIPATHWAYS") ~ "WikiPathways",
        TRUE ~ NA_character_)) %>%
    filter(!is.na(annotation_database))

  go_sets <- load_msigdbr_collection("go") %>%
    mutate(
      annotation_database = "GO",
      annotation_source = case_when(
        str_detect(gs_subcollection, "GO:BP|(^|:)BP$") ~ "GO:BP",
        str_detect(gs_subcollection, "GO:CC|(^|:)CC$") ~ "GO:CC",
        str_detect(gs_subcollection, "GO:MF|(^|:)MF$") ~ "GO:MF",
        TRUE ~ NA_character_)) %>%
    filter(!is.na(annotation_source))

  msigdb_annotations_raw <- bind_rows(
    hallmark_sets,
    pathway_sets,
    go_sets) %>%
    filter(gene_symbol %in% gene_tbl_base$gene_symbol)

  msigdb_annotations <- msigdb_annotations_raw %>%
    standardise_msigdb_annotation_cols() %>%
    bind_cols(
      msigdb_annotations_raw %>%
        dplyr::select(annotation_database, annotation_source)) %>%
    distinct(
      gene_symbol,
      annotation_database,
      annotation_source,
      gs_name,
      db_gene_symbol,
      .keep_all = TRUE)

  annotation_long <- gene_tbl_base %>%
    left_join(msigdb_annotations, by = "gene_symbol") %>%
    mutate(has_msigdb_annotation = !is.na(gs_name)) %>%
    arrange(
      gene_rank,
      annotation_database,
      annotation_source,
      gs_name,
      db_gene_symbol)

  annotation_wide <- annotation_long %>%
    group_by(gene_symbol) %>%
    summarise(
      n_msigdb_annotations = sum(has_msigdb_annotation),
      msigdb_versions = collapse_unique(db_version),
      hallmark_annotations = collapse_unique(
        gs_name[annotation_source == "Hallmark"]),
      reactome_annotations = collapse_unique(
        gs_name[annotation_source == "Reactome"]),
      kegg_annotations = collapse_unique(
        gs_name[annotation_source == "KEGG"]),
      wikipathways_annotations = collapse_unique(
        gs_name[annotation_source == "WikiPathways"]),
      go_bp_annotations = collapse_unique(
        gs_name[annotation_source == "GO:BP"]),
      go_cc_annotations = collapse_unique(
        gs_name[annotation_source == "GO:CC"]),
      go_mf_annotations = collapse_unique(
        gs_name[annotation_source == "GO:MF"]),
      go_annotations = collapse_unique(
        gs_name[annotation_database == "GO"]),
      .groups = "drop") %>%
    right_join(gene_tbl_base, by = "gene_symbol") %>%
    arrange(gene_rank) %>%
    dplyr::select(
      gene_symbol,
      gene_rank,
      gene_family,
      human_symbols,
      known_cscc_gene,
      gene_length_kb,
      n_mutated_biopsies,
      n_recurrent_mutations,
      n_nonsynonymous_mutations,
      n_msigdb_annotations,
      msigdb_versions,
      hallmark_annotations,
      reactome_annotations,
      kegg_annotations,
      wikipathways_annotations,
      go_bp_annotations,
      go_cc_annotations,
      go_mf_annotations,
      go_annotations)

  list(long = annotation_long, wide = annotation_wide)
}

load_submission_msigdb_annotations <- function(synonymous = FALSE) {
  suffix <- if (synonymous) "_synonymous" else ""
  long_hash <- if (synonymous) {
    "2a01bf56915da115aca19c6d2b531203e3189e183419f3e800dd7649dd211b7f"
  } else {
    "1f12fed501726233339c5f97b8854a8157b07947e2486f283f77ff76920ca09b"
  }
  wide_hash <- if (synonymous) {
    "56b657035d11674014e6114bfdf58e046056e89039f8516492799e624b7ef669"
  } else {
    "a564167836a9228b6f7d7795e62326e1dbcdc566d24777b1ec32137783c4df7b"
  }
  long_name <- paste0("Figure3_oncoplot_msigdb_annotations", suffix,
    "_long_data.csv")
  wide_name <- paste0("Figure3_oncoplot_msigdb_annotations", suffix,
    "_data.csv")
  long_path <- publication_repo_input(
    file.path("Figures/Figure 3", long_name), long_hash)
  wide_path <- publication_repo_input(
    file.path("Figures/Figure 3", wide_name), wide_hash)
  list(
    long = readr::read_csv(long_path, show_col_types = FALSE),
    wide = readr::read_csv(wide_path, show_col_types = FALSE),
    source = c(long = long_path, wide = wide_path),
    hash = c(long = long_hash, wide = wide_hash))
}

if (identical(EXTERNAL_DATA_MODE, "frozen")) {
  oncoplot_msigdb_annotations <- load_submission_msigdb_annotations()
  publication_copy_frozen(
    oncoplot_msigdb_annotations$source[["long"]],
    file.path(figure3_dir,
      "Figure3_oncoplot_msigdb_annotations_long_data.csv"),
    oncoplot_msigdb_annotations$hash[["long"]])
  publication_copy_frozen(
    oncoplot_msigdb_annotations$source[["wide"]],
    file.path(figure3_dir, "Figure3_oncoplot_msigdb_annotations_data.csv"),
    oncoplot_msigdb_annotations$hash[["wide"]])
} else {
  publication_require_live_refresh()
  oncoplot_msigdb_annotations <- build_oncoplot_msigdb_annotations(gene_order)
  readr::write_csv(
    oncoplot_msigdb_annotations$long,
    file.path(figure3_dir,
      "Figure3_oncoplot_msigdb_annotations_long_data.csv"))
  readr::write_csv(
    oncoplot_msigdb_annotations$wide,
    file.path(figure3_dir, "Figure3_oncoplot_msigdb_annotations_data.csv"))
}

plot_y_limits <- range(gene_order$plot_y) + c(-0.5, 0.5)

gene_group_label_df <- gene_order %>%
  group_by(gene_family) %>%
  summarise(
    plot_y = mean(range(plot_y)),
    family_mutated_biopsies = dplyr::first(family_mutated_biopsies),
    family_recurrent_mutations = dplyr::first(family_recurrent_mutations),
    .groups = "drop") %>%
  arrange(
    gene_family == "Other",
    dplyr::desc(family_mutated_biopsies),
    dplyr::desc(family_recurrent_mutations),
    gene_family)

gene_label_df <- gene_order %>%
  mutate(gene_label_strip = factor("Week 8\nMorphologically normal"))

gene_mutation_counts <- recurrent_mutations %>%
  filter(SYMBOL %in% gene_order$gene_symbol) %>%
  distinct(SYMBOL, mutation_id, mutation_effect) %>%
  group_by(gene_symbol = SYMBOL, mutation_id) %>%
  summarise(
    mutation_effect = if_else(
      any(mutation_effect == "Nonsynonymous"),
      "Nonsynonymous",
      "Synonymous"),
    .groups = "drop") %>%
  count(gene_symbol, mutation_effect, name = "n_mutations")

gene_mutation_bar_df <- expand_grid(
  gene_symbol = gene_order$gene_symbol,
  mutation_effect = mutation_effect_levels) %>%
  left_join(
    gene_mutation_counts,
    by = c("gene_symbol", "mutation_effect")) %>%
  left_join(
    gene_order %>% dplyr::select(gene_symbol, plot_y),
    by = "gene_symbol") %>%
  mutate(
    n_mutations = coalesce(n_mutations, 0L),
    bar_strip = factor("Week 8\nMorphologically normal"),
    mutation_effect = factor(
      mutation_effect,
      levels = mutation_effect_levels))

gene_mutation_totals <- gene_mutation_bar_df %>%
  pivot_wider(
    names_from = mutation_effect,
    values_from = n_mutations,
    values_fill = 0L) %>%
  mutate(total_mutations = Nonsynonymous + Synonymous) %>%
  dplyr::select(gene_symbol, Nonsynonymous, Synonymous, total_mutations)

gene_mutation_bar_df <- gene_mutation_bar_df %>%
  left_join(gene_mutation_totals, by = "gene_symbol") %>%
  mutate(pct_mutations = if_else(total_mutations > 0, n_mutations / total_mutations, 0))

gene_biopsy_hits <- recurrent_mutations %>%
  filter(
    SYMBOL %in% gene_order$gene_symbol,
    mutation_nature %in% mutation_nature_levels) %>%
  mutate(
    mutation_nature = factor(
      mutation_nature,
      levels = mutation_nature_levels)) %>%
  group_by(gene_symbol = SYMBOL, sample_name) %>%
  summarise(
    mutation_classes = paste(
      sort(unique(as.character(mutation_nature))),
      collapse = "; "),
    mutation_ids = paste(sort(unique(mutation_id)), collapse = "; "),
    n_cell_mutations = n_distinct(mutation_id),
    mutation_nature = as.character(
      mutation_nature[which.min(as.integer(mutation_nature))]),
    .groups = "drop")

muts_shared_timepoints_biopsy_df <- expand_grid(
  gene_symbol = gene_order$gene_symbol,
  sample_name = biopsy_order$sample_name) %>%
  left_join(biopsy_order, by = "sample_name") %>%
  left_join(gene_biopsy_hits, by = c("gene_symbol", "sample_name")) %>%
  left_join(gene_mutation_totals, by = "gene_symbol") %>%
  left_join(gene_order, by = "gene_symbol") %>%
  mutate(
    sample_name = factor(sample_name, levels = biopsy_order$sample_name),
    mutation_nature = factor(
      mutation_nature,
      levels = mutation_nature_levels),
    mutation_nature_plot = factor(
      coalesce(as.character(mutation_nature), "Not detected"),
      levels = c(mutation_nature_levels, "Not detected")))

muts_shared_timepoints_biopsy_bar <- ggplot(
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
    breaks = gene_order$plot_y,
    labels = NULL,
    expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(
    values = c(
      Nonsynonymous = "#496B80",
      Synonymous = "#B2585E"),
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

muts_shared_timepoints_biopsy_gene_labels <- ggplot(
  gene_label_df,
  aes(x = 1, y = plot_y, label = gene_symbol, colour = known_cscc_gene)) +
  geom_text(hjust = 1, size = plot_text_size_mm) +
  facet_grid(. ~ gene_label_strip) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    limits = plot_y_limits,
    breaks = gene_order$plot_y,
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

muts_shared_timepoints_biopsy_heatmap <- ggplot(
  muts_shared_timepoints_biopsy_df,
  aes(x = sample_name, y = plot_y, fill = mutation_nature_plot)) +
  geom_tile(color = "white", linewidth = 0.2, height = 0.85) +
  facet_grid(
    . ~ biopsy_group,
    scales = "free_x",
    space = "free_x") +
  scale_y_continuous(
    limits = plot_y_limits,
    breaks = gene_order$plot_y,
    labels = NULL,
    expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(
    values = c(
      Missense = "#3D1B1B",
      Insertion = "#7EA3CF",
      Deletion = "#436503",
      DBS = "#2F9A8B",
      "Not detected" = "#D9D9D9"),
    limits = c(mutation_nature_levels, "Not detected"),
    breaks = mutation_nature_levels,
    drop = FALSE) +
  labs(x = "", y = "", fill = "") +
  my_theme +
  theme(
    text = element_text(size = plot_text_size),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = plot_text_size,
      angle = 0,
      hjust = 0.5,
      lineheight = 0.85),
    legend.position = "bottom",
    legend.text = element_text(size = plot_text_size),
    legend.key.size = grid::unit(0.35, "cm"),
    panel.spacing.x = grid::unit(0.08, "lines"),
    plot.margin = margin(5.5, 5.5, 5.5, 0)) +
  guides(fill = guide_legend(nrow = 1))

gene_length_df <- gene_order %>%
  mutate(
    length_strip = factor("Week 8\nMorphologically normal"),
    gene_length_label = if_else(
      is.na(gene_length_kb),
      "NA",
      formatC(gene_length_kb, format = "f", digits = 2)))

muts_shared_timepoints_biopsy_gene_length <- ggplot(
  gene_length_df,
  aes(x = 0, y = plot_y, label = gene_length_label)) +
  geom_text(hjust = 0, size = plot_text_size_mm) +
  geom_text(
    data = gene_group_label_df,
    aes(x = 0.35, y = plot_y, label = gene_family),
    inherit.aes = FALSE,
    hjust = 0,
    size = plot_text_size_mm) +
  facet_grid(. ~ length_strip) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    limits = plot_y_limits,
    breaks = gene_order$plot_y,
    labels = NULL,
    expand = expansion(mult = c(0, 0))) +
  labs(x = "Gene length\n(kb)", y = "") +
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

muts_shared_timepoints_biopsy <- grid::grid.grabExpr(
  {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(
        nrow = 1,
        ncol = 4,
        widths = grid::unit(c(1.4, 0.9, 6, 2.2), "null"))))
    print(
      muts_shared_timepoints_biopsy_bar,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(
      muts_shared_timepoints_biopsy_gene_labels,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(
      muts_shared_timepoints_biopsy_heatmap,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
    print(
      muts_shared_timepoints_biopsy_gene_length,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 4))
  },
  wrap.grobs = TRUE)

readr::write_csv(
  muts_shared_timepoints_biopsy_df,
  file.path(figure3_dir,
    "Figure3_muts_shared_timepoints_biopsy_data.csv"))
grid::grid.draw(muts_shared_timepoints_biopsy)

mutation_nature_levels_synonymous <- c(
  mutation_nature_levels,
  "Synonymous",
  "Gene-flanking")
mutation_effect_levels_synonymous <- c(
  "Nonsynonymous",
  "Synonymous",
  "Gene-flanking")

mutation_biopsy_base_synonymous <- tes_mutations_df %>%
  plot_levels() %>%
  mutate(gt_AF = suppressWarnings(as.numeric(gt_AF))) %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category %in% c("SCC", "Papilloma"),
    gt_AF >= 0.01,
    !is.na(sample_name),
    !is.na(SYMBOL),
    SYMBOL != "",
    !sample_name %in% tes_outliers$sample_name) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_nature = plot_class_with_intronic_utr(
      Consequence,
      VARIANT_CLASS,
      REF,
      ALT),
    mutation_effect = plot_effect_with_intronic_utr(Consequence)) %>%
  filter(!is.na(plot_timepoint)) %>%
  filter(!(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
  semi_join(biopsy_order_all, by = "sample_name") %>%
  filter_human_orthologues() %>%
  distinct(
    sample_name,
    plot_timepoint,
    mutation_id,
    SYMBOL,
    mutation_nature,
    mutation_effect,
    .keep_all = TRUE)

# Visually normal regions are not displayed in the tumour-prone-region
# oncoplot, but are retained here so the companion table can report whether
# the exact shared tumour-region alleles were also detected in normal regions.
normal_mutation_biopsy_base_synonymous <- tes_mutations_df %>%
  plot_levels() %>%
  mutate(gt_AF = suppressWarnings(as.numeric(gt_AF))) %>%
  filter(
    treatment == "DT",
    tissue %in% c("Hair follicle", "Skin"),
    category == "Visually normal",
    gt_AF >= 0.01,
    !is.na(sample_name),
    !is.na(SYMBOL),
    SYMBOL != "",
    !sample_name %in% tes_outliers$sample_name) %>%
  mutate(
    plot_timepoint = case_when(
      tissue == "Hair follicle" ~ as.character(time),
      tissue == "Skin" ~ "Skin"),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
    mutation_nature = plot_class_with_intronic_utr(
      Consequence,
      VARIANT_CLASS,
      REF,
      ALT),
    mutation_effect = plot_effect_with_intronic_utr(Consequence)) %>%
  filter(!is.na(plot_timepoint)) %>%
  filter(!(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
  filter_human_orthologues() %>%
  distinct(
    sample_name,
    plot_timepoint,
    mutation_id,
    SYMBOL,
    mutation_nature,
    mutation_effect,
    .keep_all = TRUE)

recurrent_mutation_ids_synonymous <- mutation_biopsy_base_synonymous %>%
  distinct(mutation_id, plot_timepoint) %>%
  group_by(mutation_id) %>%
  summarise(n_timepoints = n_distinct(plot_timepoint), .groups = "drop") %>%
  filter(n_timepoints >= 2)

recurrent_mutations_synonymous <- mutation_biopsy_base_synonymous %>%
  semi_join(recurrent_mutation_ids_synonymous, by = "mutation_id")

biopsy_order_synonymous <- biopsy_order_all %>%
  semi_join(
    recurrent_mutations_synonymous %>%
      filter(mutation_nature %in% mutation_nature_levels_synonymous) %>%
      distinct(sample_name),
    by = "sample_name") %>%
  arrange(morphology, plot_timepoint, sample_name) %>%
  mutate(biopsy_rank = row_number())

gene_order_synonymous <- recurrent_mutations_synonymous %>%
  filter(mutation_nature %in% mutation_nature_levels_synonymous) %>%
  group_by(gene_symbol = SYMBOL) %>%
  summarise(
    n_mutated_biopsies = n_distinct(sample_name),
    n_recurrent_mutations = n_distinct(mutation_id),
    n_nonsynonymous_mutations = n_distinct(
      mutation_id[mutation_effect == "Nonsynonymous"]),
    first_biopsy_rank = min(
      biopsy_order_synonymous$biopsy_rank[
        match(sample_name, biopsy_order_synonymous$sample_name)
      ],
      na.rm = TRUE),
    .groups = "drop") %>%
  left_join(gene_family_assignments, by = "gene_symbol") %>%
  mutate(
    gene_family = coalesce(gene_family, "Other"),
    gene_family = factor(gene_family, levels = gene_family_levels)) %>%
  left_join(gene_lengths_kb, by = "gene_symbol") %>%
  mutate(
    epigenetic_plot_rank = coalesce(
      match(gene_symbol, epigenetic_plot_order),
      Inf)) %>%
  group_by(gene_family) %>%
  mutate(
    family_mutated_biopsies = sum(n_mutated_biopsies),
    family_recurrent_mutations = sum(n_recurrent_mutations),
    family_n_genes = n()) %>%
  ungroup() %>%
  arrange(
    gene_family == "Other",
    dplyr::desc(family_mutated_biopsies),
    dplyr::desc(family_recurrent_mutations),
    gene_family,
    epigenetic_plot_rank,
    dplyr::desc(n_mutated_biopsies),
    dplyr::desc(n_recurrent_mutations),
    dplyr::desc(n_nonsynonymous_mutations),
    first_biopsy_rank,
    gene_symbol) %>%
  mutate(
    gene_rank = row_number(),
    family_rank = as.integer(factor(gene_family, levels = unique(gene_family))),
    n_gene_families = n_distinct(gene_family),
    plot_y = n() - gene_rank + 1 + (n_gene_families - family_rank) * 0.35,
    known_cscc_gene = gene_symbol %in% KNOWN_CSCC_GENES)

if (identical(EXTERNAL_DATA_MODE, "frozen")) {
  oncoplot_msigdb_annotations_synonymous <-
    load_submission_msigdb_annotations(synonymous = TRUE)
  publication_copy_frozen(
    oncoplot_msigdb_annotations_synonymous$source[["long"]],
    file.path(figure3_dir,
      "Figure3_oncoplot_msigdb_annotations_synonymous_long_data.csv"),
    oncoplot_msigdb_annotations_synonymous$hash[["long"]])
  publication_copy_frozen(
    oncoplot_msigdb_annotations_synonymous$source[["wide"]],
    file.path(figure3_dir,
      "Figure3_oncoplot_msigdb_annotations_synonymous_data.csv"),
    oncoplot_msigdb_annotations_synonymous$hash[["wide"]])
} else {
  publication_require_live_refresh()
  oncoplot_msigdb_annotations_synonymous <-
    build_oncoplot_msigdb_annotations(gene_order_synonymous)
  readr::write_csv(
    oncoplot_msigdb_annotations_synonymous$long,
    file.path(figure3_dir,
      "Figure3_oncoplot_msigdb_annotations_synonymous_long_data.csv"))
  readr::write_csv(
    oncoplot_msigdb_annotations_synonymous$wide,
    file.path(figure3_dir,
      "Figure3_oncoplot_msigdb_annotations_synonymous_data.csv"))
}

functional_group_rule_definitions <- tibble::tribble(
  ~functional_group, ~priority, ~rule_term, ~match_scope, ~regex,
  "Keratin-associated proteins", 1L, "symbol_starts_Krtap", "symbol", "^KRTAP",
  "Keratin-associated proteins", 1L, "symbol_Gm115_cluster", "symbol", "^GM115",
  "DNA damage and repair", 2L, "DNA_DAMAGE_RESPONSE", "msigdb", "DNA_DAMAGE_RESPONSE|DNA DAMAGE RESPONSE",
  "DNA damage and repair", 2L, "DNA_REPAIR", "msigdb", "DNA_REPAIR|DNA REPAIR",
  "DNA damage and repair", 2L, "P53_PATHWAY", "msigdb", "P53_PATHWAY|P53 PATHWAY|TP53",
  "DNA damage and repair", 2L, "CHECKPOINT", "msigdb", "CHECKPOINT",
  "DNA damage and repair", 2L, "UV_RESPONSE", "msigdb", "UV_RESPONSE|UV RESPONSE",
  "DNA damage and repair", 2L, "CHROMOSOME_MAINTENANCE", "msigdb", "CHROMOSOME_MAINTENANCE|CHROMOSOME MAINTENANCE",
  "DNA damage and repair", 2L, "TELOMERE", "msigdb", "TELOMERE",
  "MAPK / Ras", 3L, "RAS", "both", "(^|[^A-Z0-9])(RAS|KRAS|HRAS|RRAS|RRAS2)([^A-Z0-9]|$)",
  "MAPK / Ras", 3L, "MAPK", "msigdb", "MAPK|MITOGEN_ACTIVATED_PROTEIN_KINASE",
  "MAPK / Ras", 3L, "ERK", "msigdb", "(^|[^A-Z0-9])ERK([^A-Z0-9]|$)|EXTRACELLULAR_SIGNAL_REGULATED_KINASE",
  "MAPK / Ras", 3L, "EGFR", "both", "EGFR|EPIDERMAL_GROWTH_FACTOR_RECEPTOR",
  "MAPK / Ras", 3L, "GROWTH_FACTOR", "msigdb", "GROWTH_FACTOR|GROWTH FACTOR",
  "MAPK / Ras", 3L, "RTK", "msigdb", "(^|[^A-Z0-9])RTK([^A-Z0-9]|$)|RECEPTOR_TYROSINE_KINASE",
  "MAPK / Ras", 3L, "ACTIVIN", "msigdb", "ACTIVIN",
  "MAPK / Ras", 3L, "PKC", "msigdb", "(^|[^A-Z0-9])PKC([^A-Z0-9]|$)|PROTEIN_KINASE_C",
  "Epigenetic regulators", 4L, "CHROMATIN", "msigdb", "CHROMATIN",
  "Epigenetic regulators", 4L, "HISTONE", "msigdb", "HISTONE",
  "Epigenetic regulators", 4L, "EPIGENETIC", "msigdb", "EPIGENETIC",
  "Epigenetic regulators", 4L, "METHYL", "msigdb", "METHYL",
  "Epigenetic regulators", 4L, "DEMETHYL", "msigdb", "DEMETHYL",
  "Epigenetic regulators", 4L, "ATP_DEPENDENT_CHROMATIN", "msigdb", "ATP_DEPENDENT_CHROMATIN|ATP DEPENDENT CHROMATIN",
  "Epigenetic regulators", 4L, "CHROMATIN_REMODELING", "msigdb", "CHROMATIN_REMODEL|CHROMATIN REMODEL",
  "Epigenetic regulators", 4L, "HAT_HDAC", "msigdb", "(^|[^A-Z0-9])HATS?([^A-Z0-9]|$)|HDAC",
  "Gene expression regulation", 5L, "TRANSCRIPTION", "msigdb", "TRANSCRIPTION|TRANSCRIPTIONAL",
  "Gene expression regulation", 5L, "RNA_POLYMERASE", "msigdb", "RNA_POLYMERASE|RNA POLYMERASE|POL_II|POL II",
  "Gene expression regulation", 5L, "TFIID_MEDIATOR", "msigdb", "TFIID|TFII_D|MEDIATOR_COMPLEX|MEDIATOR COMPLEX",
  "Gene expression regulation", 5L, "PROMOTER_ENHANCER", "msigdb", "PROMOTER|ENHANCER",
  "Gene expression regulation", 5L, "DNA_BINDING_TRANSCRIPTION_FACTOR", "msigdb", "DNA_BINDING_TRANSCRIPTION_FACTOR|DNA BINDING TRANSCRIPTION FACTOR|SEQUENCE_SPECIFIC_DNA_BINDING|SEQUENCE SPECIFIC DNA BINDING|REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING|REGULATORY REGION SEQUENCE SPECIFIC DNA BINDING",
  "Gene expression regulation", 5L, "RNA_PROCESSING", "msigdb", "RNA_PROCESSING|RNA PROCESSING|MRNA_PROCESSING|MRNA PROCESSING|RNA_SPLICING|RNA SPLICING|SPLICEOSOME",
  "Gene expression regulation", 5L, "RNA_SURVEILLANCE_DECAY", "msigdb", "NONSENSE_MEDIATED_DECAY|NONSENSE MEDIATED DECAY|MRNA_SURVEILLANCE|MRNA SURVEILLANCE|DEADENYLATION|DECAPPING|RNA_STABILITY|RNA STABILITY",
  "Gene expression regulation", 5L, "RNA_LOCALIZATION_EXPORT", "msigdb", "RNA_LOCALIZATION|RNA LOCALIZATION|RNA_EXPORT|RNA EXPORT|MRNA_EXPORT|MRNA EXPORT",
  "Gene expression regulation", 5L, "TRANSLATION", "msigdb", "TRANSLATION|TRANSLATIONAL",
  "Gene expression regulation", 5L, "RIBOSOME", "msigdb", "RIBOSOME|RIBOSOMAL",
  "Gene expression regulation", 5L, "EIF_TRANSLATION_INITIATION", "both", "(^|[^A-Z0-9])EIF[0-9A-Z]*([^A-Z0-9]|$)|EUKARYOTIC_TRANSLATION_INITIATION|EUKARYOTIC TRANSLATION INITIATION|TRANSLATION_INITIATION|TRANSLATION INITIATION",
  "Ion transport and metabolism", 6L, "METABOLISM", "msigdb", "METABOLISM",
  "Ion transport and metabolism", 6L, "XENOBIOTIC", "msigdb", "XENOBIOTIC",
  "Ion transport and metabolism", 6L, "OXIDATIVE", "msigdb", "OXIDATIVE",
  "Ion transport and metabolism", 6L, "GLYCOLYSIS", "msigdb", "GLYCOLYSIS",
  "Ion transport and metabolism", 6L, "LIPID", "msigdb", "LIPID",
  "Ion transport and metabolism", 6L, "TRANSPORT", "msigdb", "TRANSPORT",
  "Ion transport and metabolism", 6L, "SLC", "both", "(^|[^A-Z0-9])SLC([^A-Z0-9]|$)|^SLC",
  "Ion transport and metabolism", 6L, "CALCIUM", "msigdb", "CALCIUM",
  "Cell adhesion and cytoskeleton", 7L, "CELL_ADHESION", "msigdb", "CELL_ADHESION|CELL ADHESION",
  "Cell adhesion and cytoskeleton", 7L, "APICAL_JUNCTION", "msigdb", "APICAL_JUNCTION|APICAL JUNCTION",
  "Cell adhesion and cytoskeleton", 7L, "CYTOSKELETON", "msigdb", "CYTOSKELETON",
  "Cell adhesion and cytoskeleton", 7L, "ACTIN", "msigdb", "ACTIN",
  "Cell adhesion and cytoskeleton", 7L, "ECM_EXTRACELLULAR_MATRIX", "msigdb", "(^|[^A-Z0-9])ECM([^A-Z0-9]|$)|EXTRACELLULAR_MATRIX|EXTRACELLULAR MATRIX",
  "Cell adhesion and cytoskeleton", 7L, "INTEGRIN", "msigdb", "INTEGRIN",
  "Cell adhesion and cytoskeleton", 7L, "RHO_GTPASE", "msigdb", "RHO_GTPASE|RHO GTPASE",
  "Cell adhesion and cytoskeleton", 7L, "CELL_MIGRATION", "msigdb", "CELL_MIGRATION|CELL MIGRATION",
  "Cell adhesion and cytoskeleton", 7L, "EPIDERMIS_KERATINIZATION_CORNIFIED", "msigdb", "EPIDERMIS|EPIDERMAL|KERATINIZATION|CORNIFIED",
  "Immune", 8L, "INFLAMMATORY", "msigdb", "INFLAMMATORY",
  "Immune", 8L, "INTERFERON", "msigdb", "INTERFERON",
  "Immune", 8L, "NF_KAPPAB", "both", "NF_KAPPAB|NF KAPPAB|NF_KB|NFKB",
  "Immune", 8L, "TNFA", "msigdb", "TNFA|TNFALPHA|TNF_ALPHA|TNF ALPHA",
  "Immune", 8L, "CYTOKINE", "msigdb", "CYTOKINE",
  "Immune", 8L, "IMMUNE", "msigdb", "IMMUNE")

build_functional_group_rule_table <- function(gene_tbl, annotation_long, rule_tbl) {
  gene_tbl_base <- gene_tbl %>%
    dplyr::select(
      gene_symbol,
      gene_rank,
      human_symbols,
      known_cscc_gene,
      gene_length_kb,
      n_mutated_biopsies,
      n_recurrent_mutations,
      n_nonsynonymous_mutations)

  symbol_matches <- tidyr::crossing(
    gene_symbol = gene_tbl_base$gene_symbol,
    rule_tbl %>% filter(match_scope %in% c("symbol", "both"))) %>%
    mutate(
      symbol_rule_text = str_to_upper(gene_symbol),
      rule_matched = str_detect(symbol_rule_text, regex)) %>%
    filter(rule_matched)

  keratin_override <- symbol_matches %>%
    filter(
      functional_group == "Keratin-associated proteins",
      rule_term %in% c("symbol_starts_Krtap", "symbol_Gm115_cluster")) %>%
    group_by(gene_symbol) %>%
    summarise(
      override_group = "Keratin-associated proteins",
      override_rule_terms = collapse_unique(rule_term),
      .groups = "drop")

  annotation_rule_input <- annotation_long %>%
    filter(has_msigdb_annotation) %>%
    mutate(
      symbol_rule_text = str_to_upper(gene_symbol),
      msigdb_rule_text = str_to_upper(str_squish(paste(
        gs_name,
        gs_description,
        gs_collection_name,
        gs_exact_source,
        annotation_source))),
      both_rule_text = str_squish(paste(symbol_rule_text, msigdb_rule_text))) %>%
    dplyr::select(
      gene_symbol,
      gs_name,
      symbol_rule_text,
      msigdb_rule_text,
      both_rule_text)

  msigdb_rule_matches <- annotation_rule_input %>%
    tidyr::crossing(rule_tbl %>% filter(match_scope %in% c("msigdb", "both"))) %>%
    mutate(
      rule_text = case_when(
        match_scope == "msigdb" ~ msigdb_rule_text,
        match_scope == "both" ~ both_rule_text,
        TRUE ~ msigdb_rule_text),
      rule_matched = str_detect(rule_text, regex)) %>%
    filter(rule_matched)

  group_scores <- msigdb_rule_matches %>%
    group_by(gene_symbol, functional_group, priority) %>%
    summarise(
      matching_gs_names = n_distinct(gs_name),
      rule_term_score = n_distinct(rule_term),
      matched_rule_terms = collapse_unique(rule_term),
      matched_gs_names = collapse_unique(gs_name),
      .groups = "drop")

  top_scores <- group_scores %>%
    group_by(gene_symbol) %>%
    arrange(
      desc(matching_gs_names),
      priority,
      functional_group,
      .by_group = TRUE) %>%
    mutate(
      top_group = dplyr::first(functional_group),
      top_group_score = dplyr::first(matching_gs_names),
      second_group_score = if_else(
        n() >= 2,
        dplyr::nth(matching_gs_names, 2L),
        0L),
      score_margin_to_second = top_group_score - second_group_score,
      top_tie_count = sum(matching_gs_names == top_group_score),
      top_tie_groups = collapse_unique(
        functional_group[matching_gs_names == top_group_score])) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::select(
      gene_symbol,
      top_group,
      top_group_score,
      second_group_score,
      score_margin_to_second,
      top_tie_count,
      top_tie_groups)

  assigned_groups <- top_scores %>%
    left_join(keratin_override, by = "gene_symbol") %>%
    mutate(
      has_keratin_override = !is.na(override_group),
      functional_group_assignment = case_when(
        has_keratin_override ~ override_group,
        top_group_score > 0 &
          top_tie_count == 1 &
          score_margin_to_second >= 2 ~ top_group,
        top_group_score > 0 ~ "For manual review",
        TRUE ~ "Other"),
      assignment_rule = case_when(
        has_keratin_override ~ "Krtap/Gm115 symbol override",
        functional_group_assignment == "For manual review" ~
          "Top group margin < 2 or tie",
        functional_group_assignment == "Other" ~
          "No matching MSigDB gene sets or symbol override",
        TRUE ~ "Top group has >=2 more matching unique gs_names than second group")) %>%
    dplyr::select(
      gene_symbol,
      functional_group_assignment,
      assignment_rule,
      has_keratin_override,
      override_rule_terms,
      top_group,
      top_group_score,
      second_group_score,
      score_margin_to_second,
      top_tie_count,
      top_tie_groups)

  matched_group_summary <- group_scores %>%
    group_by(gene_symbol) %>%
    summarise(
      matched_group_count = n_distinct(functional_group),
      matched_multiple_groups = matched_group_count > 1,
      matched_groups = paste(
        functional_group,
        matching_gs_names,
        sep = " (gs_names=",
        collapse = "; "),
      matched_groups = str_replace_all(matched_groups, ";", ");"),
      matched_groups = paste0(matched_groups, ")"),
      .groups = "drop")

  score_wide <- group_scores %>%
    dplyr::select(gene_symbol, functional_group, matching_gs_names) %>%
    tidyr::pivot_wider(
      names_from = functional_group,
      values_from = matching_gs_names,
      values_fill = 0,
      names_prefix = "n_gs_names_") %>%
    rename_with(
      ~ str_replace_all(
        str_to_lower(.x),
        "[^a-z0-9]+",
        "_"),
      starts_with("n_gs_names_"))

  rule_term_score_wide <- group_scores %>%
    dplyr::select(gene_symbol, functional_group, rule_term_score) %>%
    tidyr::pivot_wider(
      names_from = functional_group,
      values_from = rule_term_score,
      values_fill = 0,
      names_prefix = "n_rule_terms_") %>%
    rename_with(
      ~ str_replace_all(
        str_to_lower(.x),
        "[^a-z0-9]+",
        "_"),
      starts_with("n_rule_terms_"))

  term_wide <- group_scores %>%
    dplyr::select(gene_symbol, functional_group, matched_rule_terms) %>%
    tidyr::pivot_wider(
      names_from = functional_group,
      values_from = matched_rule_terms,
      names_prefix = "terms_") %>%
    rename_with(
      ~ str_replace_all(
        str_to_lower(.x),
        "[^a-z0-9]+",
        "_"),
      starts_with("terms_"))

  gs_name_wide <- group_scores %>%
    dplyr::select(gene_symbol, functional_group, matched_gs_names) %>%
    tidyr::pivot_wider(
      names_from = functional_group,
      values_from = matched_gs_names,
      names_prefix = "gs_names_") %>%
    rename_with(
      ~ str_replace_all(
        str_to_lower(.x),
        "[^a-z0-9]+",
        "_"),
      starts_with("gs_names_"))

  gene_tbl_base %>%
    left_join(assigned_groups, by = "gene_symbol") %>%
    left_join(keratin_override, by = "gene_symbol") %>%
    mutate(
      functional_group_assignment = coalesce(
        functional_group_assignment,
        if_else(
          !is.na(override_group),
          override_group,
          "Other")),
      assignment_rule = coalesce(
        assignment_rule,
        if_else(
          !is.na(override_group),
          "Krtap/Gm115 symbol override",
          "No matching MSigDB gene sets or symbol override")),
      has_keratin_override = coalesce(
        has_keratin_override,
        !is.na(override_group)),
      override_rule_terms = coalesce(
        override_rule_terms.x,
        override_rule_terms.y)) %>%
    dplyr::select(-override_group, -override_rule_terms.x, -override_rule_terms.y) %>%
    left_join(matched_group_summary, by = "gene_symbol") %>%
    left_join(score_wide, by = "gene_symbol") %>%
    left_join(rule_term_score_wide, by = "gene_symbol") %>%
    left_join(term_wide, by = "gene_symbol") %>%
    left_join(gs_name_wide, by = "gene_symbol") %>%
    mutate(
      top_group = coalesce(top_group, "None"),
      top_group_score = coalesce(top_group_score, 0L),
      second_group_score = coalesce(second_group_score, 0L),
      score_margin_to_second = coalesce(score_margin_to_second, 0L),
      top_tie_count = coalesce(top_tie_count, 0L),
      top_tie_groups = coalesce(top_tie_groups, "None"),
      matched_group_count = coalesce(matched_group_count, 0L),
      matched_multiple_groups = coalesce(matched_multiple_groups, FALSE),
      matched_groups = coalesce(matched_groups, "None")) %>%
    arrange(gene_rank)
}

oncoplot_functional_group_rules_synonymous <-
  build_functional_group_rule_table(
    gene_tbl = gene_order_synonymous,
    annotation_long = oncoplot_msigdb_annotations_synonymous$long,
    rule_tbl = functional_group_rule_definitions)

manual_functional_group_annotations_synonymous <- tibble::tribble(
  ~gene_symbol, ~manual_functional_group_assignment,
  "Dnm2", "Cell adhesion and cytoskeleton",
  "Syne1", "Cell adhesion and cytoskeleton",
  "Rcsd1", "Cell adhesion and cytoskeleton",
  "Smg1", "DNA damage and repair",
  "Pten", "Other",
  "Huwe1", "DNA damage and repair",
  "Prkce", "Cell adhesion and cytoskeleton",
  "Bbs9", "Other",
  "Tbc1d5", "Other",
  "Exoc6b", "Other",
  "Unc13a", "Other",
  "Notch4", "Other",
  "Sufu", "Other",
  "Smarca4", "Epigenetic regulators",
  "Setd2", "DNA damage and repair",
  "Phlpp1", "MAPK / Ras",
  "Spon1", "Cell adhesion and cytoskeleton",
  "Snx25", "Other",
  "Slit3", "Other",
  "Nrsn1", "Other",
  "H2bc14", "Epigenetic regulators",
  "Fat4", "Cell adhesion and cytoskeleton",
  "Psme4", "DNA damage and repair")

oncoplot_functional_group_rules_synonymous <-
  oncoplot_functional_group_rules_synonymous %>%
  mutate(
    rule_functional_group_assignment = functional_group_assignment,
    rule_assignment_rule = assignment_rule) %>%
  left_join(
    manual_functional_group_annotations_synonymous,
    by = "gene_symbol") %>%
  mutate(
    functional_group_assignment = coalesce(
      manual_functional_group_assignment,
      rule_functional_group_assignment),
    assignment_rule = if_else(
      is.na(manual_functional_group_assignment),
      rule_assignment_rule,
      "Manual review annotation"),
    manual_annotation_applied = !is.na(manual_functional_group_assignment)) %>%
  relocate(
    rule_functional_group_assignment,
    rule_assignment_rule,
    manual_functional_group_assignment,
    manual_annotation_applied,
    .after = assignment_rule)

readr::write_csv(
  oncoplot_functional_group_rules_synonymous,
  file.path(figure3_dir,
    "Figure3_oncoplot_functional_group_rules_synonymous_data.csv"))

functional_group_levels <- c(
  "Keratin-associated proteins",
  "DNA damage and repair",
  "MAPK / Ras",
  "Epigenetic regulators",
  "Gene expression regulation",
  "Ion transport and metabolism",
  "Cell adhesion and cytoskeleton",
  "Immune",
  "Other")

apply_functional_group_order <- function(gene_tbl, functional_group_tbl) {
  gene_tbl %>%
    dplyr::select(
      -any_of(c(
        "original_gene_family",
        "family_mutated_biopsies",
        "family_recurrent_mutations",
        "family_n_genes",
        "functional_group_rank",
        "gene_rank",
        "family_rank",
        "n_gene_families",
        "plot_y"))) %>%
    mutate(original_gene_family = as.character(gene_family)) %>%
    left_join(
      functional_group_tbl %>%
        dplyr::select(gene_symbol, functional_group_assignment),
      by = "gene_symbol") %>%
    mutate(
      gene_family = coalesce(functional_group_assignment, "Other"),
      gene_family = factor(gene_family, levels = functional_group_levels),
      functional_group_rank = match(as.character(gene_family), functional_group_levels)) %>%
    group_by(gene_family) %>%
    mutate(
      family_mutated_biopsies = sum(n_mutated_biopsies),
      family_recurrent_mutations = sum(n_recurrent_mutations),
      family_n_genes = n()) %>%
    ungroup() %>%
    arrange(
      functional_group_rank,
      epigenetic_plot_rank,
      dplyr::desc(n_mutated_biopsies),
      dplyr::desc(n_recurrent_mutations),
      dplyr::desc(n_nonsynonymous_mutations),
      first_biopsy_rank,
      gene_symbol) %>%
    mutate(
      gene_rank = row_number(),
      family_rank = as.integer(factor(gene_family, levels = unique(gene_family))),
      n_gene_families = n_distinct(gene_family),
      plot_y = n() - gene_rank + 1 + (n_gene_families - family_rank) * 0.35)
}

save_recurrent_mutation_oncoplot <- function(
  mutation_base_input,
  normal_mutation_input,
  file_stub) {
  recurrent_ids <- mutation_base_input %>%
    distinct(mutation_id, plot_timepoint) %>%
    group_by(mutation_id) %>%
    summarise(n_timepoints = n_distinct(plot_timepoint), .groups = "drop") %>%
    filter(n_timepoints >= 2)

  recurrent_mutations <- mutation_base_input %>%
    semi_join(recurrent_ids, by = "mutation_id")

  biopsy_order_plot <- biopsy_order_all %>%
    semi_join(
      recurrent_mutations %>%
        filter(mutation_nature %in% mutation_nature_levels_synonymous) %>%
        distinct(sample_name),
      by = "sample_name") %>%
    arrange(morphology, plot_timepoint, sample_name) %>%
    mutate(biopsy_rank = row_number())

  gene_order_plot_base <- recurrent_mutations %>%
    filter(mutation_nature %in% mutation_nature_levels_synonymous) %>%
    group_by(gene_symbol = SYMBOL) %>%
    summarise(
      n_mutated_biopsies = n_distinct(sample_name),
      n_recurrent_mutations = n_distinct(mutation_id),
      n_nonsynonymous_mutations = n_distinct(
        mutation_id[mutation_effect == "Nonsynonymous"]),
      first_biopsy_rank = min(
        biopsy_order_plot$biopsy_rank[
          match(sample_name, biopsy_order_plot$sample_name)
        ],
        na.rm = TRUE),
      .groups = "drop") %>%
    left_join(gene_family_assignments, by = "gene_symbol") %>%
    mutate(
      gene_family = coalesce(gene_family, "Other"),
      gene_family = factor(gene_family, levels = gene_family_levels)) %>%
    left_join(gene_lengths_kb, by = "gene_symbol") %>%
    mutate(
      epigenetic_plot_rank = coalesce(
        match(gene_symbol, epigenetic_plot_order),
        Inf),
      known_cscc_gene = gene_symbol %in% KNOWN_CSCC_GENES)

  gene_order_plot <- apply_functional_group_order(
    gene_tbl = gene_order_plot_base,
    functional_group_tbl = oncoplot_functional_group_rules_synonymous)

  gene_timepoint_summary <- bind_rows(
    recurrent_mutations,
    normal_mutation_input %>%
      semi_join(recurrent_ids, by = "mutation_id")) %>%
    filter(
      SYMBOL %in% gene_order_plot$gene_symbol,
      mutation_nature %in% mutation_nature_levels_synonymous) %>%
    transmute(
      Gene = SYMBOL,
      morphology = recode(
        as.character(category),
        "Visually normal" = "Normal"),
      timepoint = recode(
        plot_timepoint,
        "Week 8" = "w8",
        "Week 14" = "w14",
        "Week 17" = "w17",
        "Skin" = "SK"),
      timepoint_rank = match(
        plot_timepoint,
        c("Week 8", "Week 14", "Week 17", "Skin"))) %>%
    filter(
      morphology %in% c("Normal", "Papilloma", "SCC"),
      !is.na(timepoint_rank)) %>%
    distinct(Gene, morphology, timepoint, timepoint_rank) %>%
    arrange(Gene, morphology, timepoint_rank) %>%
    group_by(Gene, morphology) %>%
    summarise(
      detected_timepoints = paste(timepoint, collapse = "; "),
      .groups = "drop") %>%
    right_join(
      expand_grid(
        Gene = gene_order_plot$gene_symbol,
        morphology = c("Normal", "Papilloma", "SCC")),
      by = c("Gene", "morphology")) %>%
    mutate(detected_timepoints = coalesce(detected_timepoints, "")) %>%
    pivot_wider(
      names_from = morphology,
      values_from = detected_timepoints) %>%
    arrange(match(Gene, gene_order_plot$gene_symbol)) %>%
    dplyr::select(Gene, Normal, Papilloma, SCC)

  readr::write_csv(
    gene_timepoint_summary,
    paste0(file_stub, "_gene_timepoint_summary.csv"))

  plot_y_limits_plot <- range(gene_order_plot$plot_y) + c(-0.5, 0.5)

  gene_group_label_df_plot <- gene_order_plot %>%
    group_by(gene_family) %>%
    summarise(
      plot_y = mean(range(plot_y)),
      family_mutated_biopsies = dplyr::first(family_mutated_biopsies),
      family_recurrent_mutations = dplyr::first(family_recurrent_mutations),
      .groups = "drop") %>%
    arrange(
      gene_family == "Other",
      dplyr::desc(family_mutated_biopsies),
      dplyr::desc(family_recurrent_mutations),
      gene_family)

  gene_label_df_plot <- gene_order_plot %>%
    mutate(gene_label_strip = factor("Week 8\nMorphologically normal"))

  gene_mutation_counts_plot <- recurrent_mutations %>%
    filter(SYMBOL %in% gene_order_plot$gene_symbol) %>%
    distinct(SYMBOL, mutation_id, mutation_effect) %>%
    group_by(gene_symbol = SYMBOL, mutation_id) %>%
    summarise(
      mutation_effect = case_when(
        any(mutation_effect == "Nonsynonymous") ~ "Nonsynonymous",
        any(mutation_effect == "Gene-flanking") ~ "Gene-flanking",
        TRUE ~ "Synonymous"),
      .groups = "drop") %>%
    count(gene_symbol, mutation_effect, name = "n_mutations")

  gene_mutation_bar_df_plot <- expand_grid(
    gene_symbol = gene_order_plot$gene_symbol,
    mutation_effect = mutation_effect_levels_synonymous) %>%
    left_join(gene_mutation_counts_plot, by = c("gene_symbol", "mutation_effect")) %>%
    left_join(gene_order_plot %>% dplyr::select(gene_symbol, plot_y), by = "gene_symbol") %>%
    mutate(
      n_mutations = coalesce(n_mutations, 0L),
      bar_strip = factor("Week 8\nMorphologically normal"),
      mutation_effect = factor(mutation_effect, levels = mutation_effect_levels_synonymous))

  gene_mutation_totals_plot <- gene_mutation_bar_df_plot %>%
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

  gene_mutation_bar_df_plot <- gene_mutation_bar_df_plot %>%
    left_join(gene_mutation_totals_plot, by = "gene_symbol") %>%
    mutate(pct_mutations = if_else(total_mutations > 0, n_mutations / total_mutations, 0))

  gene_biopsy_hits_plot <- recurrent_mutations %>%
    filter(
      SYMBOL %in% gene_order_plot$gene_symbol,
      mutation_nature %in% mutation_nature_levels_synonymous) %>%
    mutate(
      mutation_nature = factor(mutation_nature, levels = mutation_nature_levels_synonymous)) %>%
    group_by(gene_symbol = SYMBOL, sample_name) %>%
    summarise(
      mutation_classes = paste(sort(unique(as.character(mutation_nature))), collapse = "; "),
      mutation_ids = paste(sort(unique(mutation_id)), collapse = "; "),
      n_cell_mutations = n_distinct(mutation_id),
      mutation_nature = as.character(
        mutation_nature[which.min(as.integer(mutation_nature))]),
      .groups = "drop")

  plot_df <- expand_grid(
    gene_symbol = gene_order_plot$gene_symbol,
    sample_name = biopsy_order_plot$sample_name) %>%
    left_join(biopsy_order_plot, by = "sample_name") %>%
    left_join(gene_biopsy_hits_plot, by = c("gene_symbol", "sample_name")) %>%
    left_join(gene_mutation_totals_plot, by = "gene_symbol") %>%
    left_join(gene_order_plot, by = "gene_symbol") %>%
    mutate(
      sample_name = factor(sample_name, levels = biopsy_order_plot$sample_name),
      mutation_nature = factor(mutation_nature, levels = mutation_nature_levels_synonymous),
      mutation_nature_plot = factor(
        coalesce(as.character(mutation_nature), "Not detected"),
        levels = c(mutation_nature_levels_synonymous, "Not detected")))

  bar_plot <- ggplot(
    gene_mutation_bar_df_plot,
    aes(x = pct_mutations, y = plot_y, fill = mutation_effect)) +
    geom_col(width = 0.85, orientation = "y") +
    facet_grid(. ~ bar_strip) +
    scale_x_reverse(
      limits = c(1, 0),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = plot_y_limits_plot,
      breaks = gene_order_plot$plot_y,
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
      strip.text.x = element_text(size = plot_text_size, angle = 90, hjust = 0, colour = "white"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)) +
    guides(fill = guide_legend(ncol = 1))

  gene_labels_plot <- ggplot(
    gene_label_df_plot,
    aes(x = 1, y = plot_y, label = gene_symbol, colour = known_cscc_gene)) +
    geom_text(hjust = 1, size = plot_text_size_mm) +
    facet_grid(. ~ gene_label_strip) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      limits = plot_y_limits_plot,
      breaks = gene_order_plot$plot_y,
      labels = NULL,
      expand = expansion(mult = c(0, 0))) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "#D63B2A"), guide = "none") +
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
      strip.text.x = element_text(size = plot_text_size, angle = 90, hjust = 0, colour = "white"),
      plot.margin = margin(5.5, 2, 5.5, 0))

  heatmap_plot <- ggplot(
    plot_df,
    aes(x = sample_name, y = plot_y, fill = mutation_nature_plot)) +
    geom_tile(color = "white", linewidth = 0.2, height = 0.85) +
    facet_grid(. ~ biopsy_group, scales = "free_x", space = "free_x") +
    scale_y_continuous(
      limits = plot_y_limits_plot,
      breaks = gene_order_plot$plot_y,
      labels = NULL,
      expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(
      values = c(
        Missense = "#3D1B1B",
        Insertion = "#7EA3CF",
        Deletion = "#436503",
        DBS = "#2F9A8B",
        Synonymous = "#B2585E",
        `Gene-flanking` = "#C7A439",
        "Not detected" = "#D9D9D9"),
      limits = c(mutation_nature_levels_synonymous, "Not detected"),
      breaks = mutation_nature_levels_synonymous,
      drop = FALSE) +
    labs(x = "", y = "", fill = "") +
    my_theme +
    theme(
      text = element_text(size = plot_text_size),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(size = plot_text_size, angle = 0, hjust = 0.5, lineheight = 0.85),
      legend.position = "bottom",
      legend.text = element_text(size = plot_text_size),
      legend.key.size = grid::unit(0.35, "cm"),
      panel.spacing.x = grid::unit(0.08, "lines"),
      plot.margin = margin(5.5, 5.5, 5.5, 0)) +
    guides(fill = guide_legend(nrow = 1))

  gene_length_plot <- gene_order_plot %>%
    mutate(
      length_strip = factor("Week 8\nMorphologically normal"),
      gene_length_label = if_else(
        is.na(gene_length_kb),
        "NA",
        formatC(gene_length_kb, format = "f", digits = 2))) %>%
    ggplot(aes(x = 0, y = plot_y, label = gene_length_label)) +
    geom_text(hjust = 0, size = plot_text_size_mm) +
    geom_text(
      data = gene_group_label_df_plot,
      aes(x = 0.35, y = plot_y, label = gene_family),
      inherit.aes = FALSE,
      hjust = 0,
      size = plot_text_size_mm) +
    facet_grid(. ~ length_strip) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      limits = plot_y_limits_plot,
      breaks = gene_order_plot$plot_y,
      labels = NULL,
      expand = expansion(mult = c(0, 0))) +
    labs(x = "Gene length\n(kb)", y = "") +
    my_theme +
    theme(
      text = element_text(size = plot_text_size),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(size = plot_text_size, angle = 90, hjust = 0, colour = "white"),
      plot.margin = margin(5.5, 5.5, 5.5, 2))

  plot_grob <- grid::grid.grabExpr(
    {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(
        layout = grid::grid.layout(
          nrow = 1,
          ncol = 4,
          widths = grid::unit(c(1.4, 0.9, 6, 2.2), "null"))))
      print(bar_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
      print(gene_labels_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
      print(heatmap_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
      print(gene_length_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 4))
    },
    wrap.grobs = TRUE)

  export_grid_plot_data(
    plot_grob = plot_grob,
    data = plot_df,
    file_name = file_stub,
    width = max(14, min(24, 7 + 0.16 * nrow(biopsy_order_plot))),
    height = max(5, min(12, 2 + 0.22 * nrow(gene_order_plot))))

  invisible(plot_grob)
}

muts_shared_timepoints_biopsy_with_synonymous_vaf <-
  save_recurrent_mutation_oncoplot(
    mutation_base_input = mutation_biopsy_base_synonymous,
    normal_mutation_input = normal_mutation_biopsy_base_synonymous,
    file_stub = file.path(figure3_dir,
      "Figure3_muts_shared_timepoints_biopsy_with_synonymous_vaf"))
grid::grid.draw(muts_shared_timepoints_biopsy_with_synonymous_vaf)

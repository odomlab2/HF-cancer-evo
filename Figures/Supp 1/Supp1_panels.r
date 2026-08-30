# 09/04/2026
# Author: Yoav Avi-Guy
# Purpose: Supplementary Figure 1 panels - Sequencing QC for TES and WES

#-------------------------------------------------------------------------------
# Libraries
library(ggplot2)
library(dplyr)
library(ggh4x)
repo_root <- normalizePath(Sys.getenv("HF_SCC_ROOT", unset = "."), mustWork = TRUE)
output_root <- Sys.getenv("HF_SCC_OUTPUT_ROOT")
if (!nzchar(output_root)) stop("Set HF_SCC_OUTPUT_ROOT to an external output directory.")

# Colours
source(file.path(repo_root, "TES", "01_source", "01_plotting.r"))

# Loading files
TES_seq_df <- readRDS(file.path(repo_root, "TES/02_data/02_processed/seq_qc_summary.rds"))
WES_seq_df <- readRDS(file.path(repo_root, "WES/02_data/02_processed/seq_qc_summary.rds"))
WES_seq_df <- WES_seq_df %>% filter(last_cycle == TRUE)

supp1_dir <- file.path(output_root, "Figures", "Supp 1")
dir.create(supp1_dir, recursive = TRUE, showWarnings = FALSE)

seq_qc_cols <- c(
  "sample_name",
  "tissue",
  "treatment",
  "unique_pct",
  "non_duplicate_mapped_reads",
  "on_target_ratio",
  "coverage",
  "qcBasesMapped")

TES_plot_cols <- c(
  seq_qc_cols,
  "y_plot_reads",
  "above_cap_reads",
  "y_plot_cov",
  "above_cap_cov")

#-------------------------------------------------------------------------------
### Mapping rate for TES ###
#-------------------------------------------------------------------------------

TES_map_rate <- ggplot(
  TES_seq_df,
  aes(
    x = treatment,
    y = unique_pct,
    fill = treatment)) +
  scale_fill_manual(values = col_palette$treatment) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ tissue, scales = "free_y") +
  my_theme +
  labs(
    x = "", 
    y = "Mapped reads (%)") +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75),
    size = 0.5, 
    alpha = 0.6) +
  theme(legend.position = "") + 
  expand_limits(y = 0)
export_plot_data(
  data = TES_seq_df,
  file_name = file.path(supp1_dir, "Supp1_TES_map_rate"),
  cols = seq_qc_cols)
TES_map_rate

WES_map_rate <- ggplot(
  WES_seq_df,
  aes(
    x = treatment,
    y = unique_pct,
    fill = treatment)) +
  scale_fill_manual(values = col_palette$treatment) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ tissue, scales = "free_y") +
  my_theme +
  labs(
    x = "", 
    y = "Mapped reads (%)") +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75),
    size = 0.5, 
    alpha = 0.6) +
  theme(legend.position = "") + 
  expand_limits(y = 0)
export_plot_data(
  data = WES_seq_df,
  file_name = file.path(supp1_dir, "Supp1_WES_map_rate"),
  cols = c(seq_qc_cols, "last_cycle"))
WES_map_rate

#-------------------------------------------------------------------------------
### Unique reads ###
#-------------------------------------------------------------------------------
TES_plot_df <- TES_seq_df %>%
  mutate(
    y_plot_reads = if_else(
      tissue == "Hair follicle",
      pmin(non_duplicate_mapped_reads, 2.5e7),
      non_duplicate_mapped_reads),
    above_cap_reads = tissue == 
      "Hair follicle" & non_duplicate_mapped_reads > 2.5e7,
    y_plot_cov = if_else(
      tissue == "Hair follicle",
      pmin(coverage, 400),
      coverage),
    above_cap_cov = tissue ==
      "Hair follicle" & coverage > 400)

TES_unique_reads <- ggplot(
  TES_plot_df,
  aes(
    x = treatment,
    y = y_plot_reads,
    fill = treatment)) +
  facet_wrap(~ tissue, scales = "free_y") +
  ggh4x::facetted_pos_scales(
    y = list(
      tissue == "Hair follicle" ~ scale_y_continuous(
        limits = c(0, 2.5e7),
        labels = my_label
      ),
      tissue != "Hair follicle" ~ scale_y_continuous(
        labels = my_label))) +
  scale_fill_manual(values = col_palette$treatment) +
  geom_boxplot(outlier.shape = NA) +
  my_theme +
  labs(
    x = "", 
    y = "# mapped reads") +
  geom_jitter(
    data = ~ dplyr::filter(.x, !above_cap_reads),
    position = position_jitterdodge(
      jitter.width = 0.35,
      dodge.width = 0.75),
    size = 0.5,
    alpha = 0.6,
    shape = 16) +
  geom_jitter(
    data = ~ dplyr::filter(.x, above_cap_reads),
    position = position_jitterdodge(
      jitter.width = 0.35,
      dodge.width = 0.75),
    size = 1.8,
    alpha = 0.9,
    shape = 17) +
  theme(legend.position = "none") + 
  expand_limits(y = 0)
export_plot_data(
  data = TES_plot_df,
  file_name = file.path(supp1_dir, "Supp1_TES_unique_reads"),
  cols = TES_plot_cols)
TES_unique_reads

WES_unique_reads <- ggplot(
  WES_seq_df,
  aes(
    x = treatment,
    y = non_duplicate_mapped_reads,
    fill = treatment)) +
  facet_wrap(~ tissue, scales = "free_y") +
  scale_fill_manual(values = col_palette$treatment) +
  geom_boxplot(outlier.shape = NA) +
  my_theme +
  labs(
    x = "", 
    y = "# mapped reads") +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75),
    size = 0.5, 
    alpha = 0.6) +
  theme(legend.position = "none") + 
  scale_y_continuous(labels = my_label) +
  expand_limits(y = 0)
export_plot_data(
  data = WES_seq_df,
  file_name = file.path(supp1_dir, "Supp1_WES_unique_reads"),
  cols = c(seq_qc_cols, "last_cycle"))
WES_unique_reads

#-------------------------------------------------------------------------------
### On-target ratio ###
#-------------------------------------------------------------------------------

TES_OTR <- ggplot(
  TES_seq_df,
  aes(
    x = treatment,
    y = on_target_ratio,
    fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75),
    size = 0.5, 
    alpha = 0.6) +
  facet_wrap(~ tissue, axes = "all_y") +
  scale_fill_manual(values = col_palette$treatment) +
  labs(
    x = "", 
    y = "On-target ratio") +
  ylim(0, 1) +
  my_theme +
  theme(legend.position = "") 
export_plot_data(
  data = TES_seq_df,
  file_name = file.path(supp1_dir, "Supp1_TES_OTR"),
  cols = seq_qc_cols)
TES_OTR

WES_OTR <- ggplot(
  WES_seq_df,
  aes(
    x = treatment,
    y = on_target_ratio,
    fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75),
    size = 0.5, 
    alpha = 0.6) +
  facet_wrap(~ tissue, axes = "all_y") +
  scale_fill_manual(values = col_palette$treatment) +
  labs(
    x = "", 
    y = "On-target ratio") +
  ylim(0, 1) +
  my_theme +
  theme(legend.position = "") 
export_plot_data(
  data = WES_seq_df,
  file_name = file.path(supp1_dir, "Supp1_WES_OTR"),
  cols = c(seq_qc_cols, "last_cycle"))
WES_OTR

#-------------------------------------------------------------------------------
### Coverage ###
#-------------------------------------------------------------------------------

TES_coverage <- ggplot(
  TES_plot_df,
  aes(
    x = treatment,
    y = y_plot_cov,
    fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    data = ~ dplyr::filter(.x, !above_cap_cov),
    position = position_jitterdodge(
      jitter.width = 0.35,
      dodge.width = 0.75),
    size = 0.5,
    alpha = 0.6,
    shape = 16) +
  geom_jitter(
    data = ~ dplyr::filter(.x, above_cap_cov),
    position = position_jitterdodge(
      jitter.width = 0.35,
      dodge.width = 0.75),
    size = 1.8,
    alpha = 0.9,
    shape = 17) +
  facet_wrap(~ tissue, scales = "free_y") +
  scale_fill_manual(values = col_palette$treatment) +
  labs(
    x = "", 
    y = "Coverage (X)") +
  my_theme +
  theme(legend.position = "")
export_plot_data(
  data = TES_plot_df,
  file_name = file.path(supp1_dir, "Supp1_TES_coverage"),
  cols = TES_plot_cols)
TES_coverage

WES_coverage <- ggplot(
  WES_seq_df,
  aes(
    x = treatment,
    y = coverage,
    fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75),
    size = 0.5, 
    alpha = 0.6) +
  facet_wrap(~ tissue, axes = "all_y") +
  scale_fill_manual(values = col_palette$treatment) +
  labs(
    x = "", 
    y = "Coverage (X)") +
  my_theme +
  theme(legend.position = "")
export_plot_data(
  data = WES_seq_df,
  file_name = file.path(supp1_dir, "Supp1_WES_coverage"),
  cols = c(seq_qc_cols, "last_cycle"))
WES_coverage

#-------------------------------------------------------------------------------
### Bases covered ###
#-------------------------------------------------------------------------------

TES_bases_covered <- ggplot(
  TES_seq_df,
  aes(
    x = treatment,
    y = qcBasesMapped,
    fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75),
    size = 0.5, 
    alpha = 0.6) +
  facet_wrap(~ tissue, axes = "all_y", scales = "free_y") +
  scale_fill_manual(values = col_palette$treatment) +
  labs(
    x = "", 
    y = "# of bases covered") +
  my_theme +
  scale_y_continuous(labels = my_label) +
  theme(legend.position = "")
export_plot_data(
  data = TES_seq_df,
  file_name = file.path(supp1_dir, "Supp1_TES_bases_covered"),
  cols = seq_qc_cols)
TES_bases_covered

WES_bases_covered <- ggplot(
  WES_seq_df,
  aes(
    x = treatment,
    y = qcBasesMapped,
    fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.35, 
      dodge.width = 0.75),
    size = 0.5, 
    alpha = 0.6) +
  facet_wrap(~ tissue, axes = "all_y", scales = "free_y") +
  scale_fill_manual(values = col_palette$treatment) +
  labs(
    x = "", 
    y = "# of bases covered") +
  my_theme +
  scale_y_continuous(labels = my_label) +
  theme(legend.position = "")
export_plot_data(
  data = WES_seq_df,
  file_name = file.path(supp1_dir, "Supp1_WES_bases_covered"),
  cols = c(seq_qc_cols, "last_cycle"))
WES_bases_covered

#-------------------------------------------------------------------------------

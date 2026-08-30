# Reproduces the submitted Figure 4F calculation.
# Set HF_SCC_ROOT and FIGURE3_OUTPUT_DIR to isolate inputs and outputs.
local({
  library(dplyr)
  library(tidyr)

  repo_root <- normalizePath(Sys.getenv(
    "HF_SCC_ROOT",
    unset = "."))
  output_dir <- Sys.getenv(
    "FIGURE3_OUTPUT_DIR",
    unset = "")
  if (!nzchar(output_dir)) {
    stop("Set FIGURE3_OUTPUT_DIR to an external output directory.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  source(file.path(repo_root, "TES", "01_source", "01_plotting.r"))

  mutations <- readRDS(file.path(
    repo_root, "TES", "02_data", "02_processed", "mutations_unique_dp.rds"))
  outliers <- readRDS(file.path(
    repo_root, "TES", "02_data", "02_processed", "technical_outliers.rds"))
  timepoints <- c("Week 8", "Week 14", "Week 17", "Skin")
  plot_timepoints <- timepoints[1:3]
  classes <- c(
    "Shared with SK", "Shared with later HF", "Shared with earlier HF")

  base <- mutations %>%
    plot_levels() %>%
    filter(
      !sample_name %in% outliers$sample_name,
      gt_AF >= 0.01,
      !(tissue == "Hair follicle" & as.character(time) == "Week 19"))
  presence <- base %>%
    filter(
      treatment == "DT",
      tissue %in% c("Hair follicle", "Skin"),
      category %in% c("Papilloma", "SCC")) %>%
    mutate(
      plot_timepoint = if_else(tissue == "Hair follicle", as.character(time), "Skin"),
      mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
      mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
    filter(!is.na(plot_timepoint),
           !(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
    distinct(mutation_label, plot_timepoint, tissue, category) %>%
    mutate(
      plot_timepoint = factor(plot_timepoint, levels = timepoints),
      timepoint_rank = as.integer(plot_timepoint))
  flags <- base %>%
    filter(
      treatment == "DT",
      tissue %in% c("Hair follicle", "Skin"),
      category %in% c("Papilloma", "SCC")) %>%
    mutate(
      plot_timepoint = if_else(tissue == "Hair follicle", as.character(time), "Skin"),
      mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
      mutation_label = paste(mutation_id, SYMBOL, sep = "_")) %>%
    filter(!is.na(plot_timepoint),
           !(tissue == "Hair follicle" & plot_timepoint == "Week 19")) %>%
    distinct(sample_name, mutation_label, plot_timepoint) %>%
    mutate(
      plot_timepoint = factor(plot_timepoint, levels = timepoints),
      timepoint_rank = as.integer(plot_timepoint)) %>%
    filter(plot_timepoint %in% plot_timepoints) %>%
    mutate(shared_skin = mutation_label %in% presence$mutation_label[
      presence$plot_timepoint == "Skin" & presence$category %in% c("SCC", "Papilloma")]) %>%
    rowwise() %>%
    mutate(
      shared_later_hf = any(presence$mutation_label == mutation_label &
        presence$tissue == "Hair follicle" & presence$timepoint_rank > timepoint_rank),
      shared_earlier_hf = any(presence$mutation_label == mutation_label &
        presence$tissue == "Hair follicle" & presence$timepoint_rank < timepoint_rank)) %>%
    ungroup() %>%
    mutate(sharing_class = case_when(
      shared_skin ~ classes[1],
      shared_later_hf ~ classes[2],
      shared_earlier_hf ~ classes[3],
      TRUE ~ NA_character_)) %>%
    filter(!is.na(sharing_class))

  totals <- count(flags, plot_timepoint, sample_name, name = "total_mutations")
  sample_pct <- totals %>%
    select(plot_timepoint, sample_name, total_mutations) %>%
    tidyr::crossing(sharing_class = factor(classes, levels = classes)) %>%
    left_join(
      flags %>% count(plot_timepoint, sample_name, sharing_class, name = "n_mutations") %>%
        mutate(sharing_class = factor(sharing_class, levels = classes)),
      by = c("plot_timepoint", "sample_name", "sharing_class")) %>%
    mutate(
      n_mutations = coalesce(n_mutations, 0L),
      pct_mutations = n_mutations / total_mutations)
  weighted_sd <- function(x, w) {
    valid <- !is.na(x) & !is.na(w) & w > 0
    sqrt(sum(w[valid] * (x[valid] - weighted.mean(x[valid], w[valid]))^2) /
      sum(w[valid]))
  }
  stacked <- sample_pct %>%
    group_by(plot_timepoint, sharing_class) %>%
    summarise(
      n_samples = n_distinct(sample_name),
      weighted_mean_pct = weighted.mean(pct_mutations, total_mutations),
      weighted_sd_pct = weighted_sd(pct_mutations, total_mutations),
      .groups = "drop") %>%
    mutate(
      pct_mutations = weighted_mean_pct,
      sd_pct_mutations = coalesce(weighted_sd_pct, 0),
      plot_timepoint_num = as.integer(plot_timepoint),
      xmin = plot_timepoint_num - 0.375,
      xmax = plot_timepoint_num + 0.375) %>%
    arrange(plot_timepoint, sharing_class) %>%
    group_by(plot_timepoint) %>%
    mutate(
      segment_ymax = cumsum(pct_mutations),
      segment_ymin = lag(segment_ymax, default = 0),
      segment_mid = segment_ymin + pct_mutations / 2,
      errorbar_ymin = pmax(segment_mid - sd_pct_mutations / 2, 0),
      errorbar_ymax = pmin(segment_mid + sd_pct_mutations / 2, 1)) %>%
    ungroup()
  plot <- ggplot(stacked, aes(fill = sharing_class)) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = segment_ymin, ymax = segment_ymax),
              colour = "black") +
    geom_errorbar(
      data = filter(stacked, sd_pct_mutations > 0),
      aes(x = plot_timepoint_num, ymin = errorbar_ymin, ymax = errorbar_ymax,
          group = sharing_class),
      inherit.aes = FALSE, width = 0.18, colour = "black") +
    scale_fill_manual(values = c(
      "Shared with SK" = "#f0b981", "Shared with later HF" = "#300358",
      "Shared with earlier HF" = "#9780aa"), drop = FALSE) +
    scale_x_continuous(breaks = seq_along(plot_timepoints), labels = plot_timepoints,
                       expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "", y = "% of shared mutations", fill = "") + my_theme +
    theme(legend.position = "bottom") +
    ggh4x::force_panelsizes(cols = grid::unit(29, "mm"), rows = grid::unit(29, "mm"))
  export <- select(stacked, plot_timepoint, sharing_class, n_samples, pct_mutations,
                   sd_pct_mutations, segment_ymin, segment_ymax, segment_mid,
                   errorbar_ymin, errorbar_ymax)
  readr::write_csv(export, file.path(output_dir,
    "Figure3_mutations_over_time_space_data.csv"))
  ggsave(file.path(output_dir, "Figure3_mutations_over_time_space.pdf"), plot,
         width = 60, height = 60, units = "mm")
})

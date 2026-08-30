suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(readr)})
root <- normalizePath(Sys.getenv("HF_SCC_ROOT", unset = "."), mustWork = TRUE)
fig_dir <- file.path(root, "Figures", "Figure 4")
out_dir <- Sys.getenv("FIGURE4_OUTPUT_DIR", unset = fig_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "TES", "01_source", "01_plotting.r"))
source(file.path(fig_dir, "Figure4_TES_WES_candidate_burden_helpers.r"))
source(file.path(fig_dir, "Figure4_TES_WES_candidate_burden_data.r"))

nonsynonymous_test <- run_candidate_burden_spearman(
  plot_data,
  y = "nonsynonymous_mutations_per_mbp",
  analysis_id = "TES_WES_candidate_count_nonsynonymous_burden")
all_mutation_test <- run_candidate_burden_spearman(
  plot_data,
  y = "all_mutations_per_mbp",
  analysis_id = "TES_WES_candidate_count_all_mutation_burden")
nonsynonymous_plot <- plot_candidate_burden_spearman(
  plot_data,
  y = "nonsynonymous_mutations_per_mbp",
  y_label = "Nonsynonymous mutations / Mb",
  test_data = nonsynonymous_test)
all_mutation_plot <- plot_candidate_burden_spearman(
  plot_data,
  y = "all_mutations_per_mbp",
  y_label = "Mutations / Mb",
  test_data = all_mutation_test)

write_csv(
  nonsynonymous_test,
  file.path(
    out_dir,
    "Figure4_TES_WES_candidate_gene_count_TMB_spearman_tests_data.csv"))
write_csv(
  all_mutation_test,
  file.path(
    out_dir,
    paste0(
      "Figure4_TES_WES_candidate_gene_count_all_mutations",
      "_spearman_tests_data.csv")))
ggplot2::set_last_plot(nonsynonymous_plot)
export_plot_data(
  data = plot_data,
  file_name = file.path(
    out_dir,
    paste0(
      "Figure4_TES_WES_candidate_gene_count_vs_",
      "nonsynonymous_TMB_spearman")))
ggplot2::set_last_plot(all_mutation_plot)
export_plot_data(
  data = plot_data,
  file_name = file.path(
    out_dir,
    paste0(
      "Figure4_TES_WES_candidate_gene_count_vs_",
      "all_mutations_per_Mb_spearman")))

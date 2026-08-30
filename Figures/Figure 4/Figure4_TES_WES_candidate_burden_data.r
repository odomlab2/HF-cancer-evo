suppressPackageStartupMessages({library(dplyr); library(readr); library(stringr)})
root <- normalizePath(Sys.getenv("HF_SCC_ROOT", unset = "."), mustWork = TRUE)
fig_dir <- file.path(root, "Figures", "Figure 4")
path <- function(assay, section, file) {
  file.path(root, assay, "02_data", section, file)
}
mutations <- bind_rows(
  TES = readRDS(path("TES", "02_processed", "mutations_unique_dp.rds")),
  WES = readRDS(path("WES", "02_processed", "mutations_unique_dp.rds")),
  .id = "dataset")
metadata <- bind_rows(
  TES = readRDS(path("TES", "01_metadata", "sample_metadata.rds")),
  WES = readRDS(path("WES", "01_metadata", "sample_metadata.rds")),
  .id = "dataset")
outliers <- bind_rows(
  TES = readRDS(path("TES", "02_processed", "technical_outliers.rds")),
  WES = readRDS(path("WES", "02_processed", "technical_outliers.rds")),
  .id = "dataset")
candidates <- read_csv(
  file.path(fig_dir, "Figure4_candidate_genes_data.csv"),
  show_col_types = FALSE) %>%
  distinct(gene_symbol)
tes_panel <- readRDS(path("TES", "02_processed", "gene_panel.rds"))
tes_targeted_candidates <- semi_join(
  candidates, tes_panel, by = c("gene_symbol" = "gene"))
stopifnot(nrow(candidates) == 30L, nrow(tes_targeted_candidates) == 8L)
protein_altering <- paste(c(
  "missense_variant", "frameshift_variant", "protein_altering_variant",
  "inframe_insertion", "inframe_deletion", "splice_acceptor_variant",
  "splice_donor_variant", "splice_region_variant", "start_lost",
  "stop_gained", "stop_lost", "coding_sequence_variant"), collapse = "|")
analysis_mutations <- mutations %>%
  anti_join(outliers, by = c("dataset", "sample_name")) %>%
  mutate(gt_AF = suppressWarnings(as.numeric(gt_AF)),
    mutation_id = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
  filter(dataset == "TES" | last_cycle %in% TRUE,
    !is.na(sample_name), !is.na(gt_AF), gt_AF >= 0.01)
nonsynonymous <- analysis_mutations %>%
  filter(!((category %in% "Visually normal" | treatment %in% "Acetone") &
      IMPACT %in% c("LOW", "MODIFIER")),
    !is.na(SYMBOL), SYMBOL != "", Consequence != "intergenic_variant",
    Feature_type != "RegulatoryFeature",
    !str_detect(coalesce(Consequence, ""), "synonymous_variant"),
    IMPACT %in% c("HIGH", "MODERATE") |
      str_detect(coalesce(Consequence, ""), protein_altering))
candidate_counts <- nonsynonymous %>%
  semi_join(candidates, by = c("SYMBOL" = "gene_symbol")) %>%
  distinct(dataset, sample_name, SYMBOL) %>%
  summarise(n_mutated_candidate_genes = n(),
    mutated_candidate_genes = paste(sort(SYMBOL), collapse = ";"),
    .by = c(dataset, sample_name))
nonsynonymous_counts <- nonsynonymous %>%
  distinct(dataset, sample_name, mutation_id) %>%
  count(dataset, sample_name, name = "n_nonsynonymous_mutations")
all_mutation_counts <- analysis_mutations %>%
  distinct(dataset, sample_name, mutation_id) %>%
  count(dataset, sample_name, name = "n_all_mutations")
plot_data <- metadata %>%
  anti_join(outliers, by = c("dataset", "sample_name")) %>%
  filter(dataset == "TES" | last_cycle %in% TRUE) %>%
  left_join(nonsynonymous_counts, by = c("dataset", "sample_name")) %>%
  left_join(all_mutation_counts, by = c("dataset", "sample_name")) %>%
  left_join(candidate_counts, by = c("dataset", "sample_name")) %>%
  mutate(sample_uid = paste(dataset, sample_name, sep = "::"),
    n_nonsynonymous_mutations = coalesce(n_nonsynonymous_mutations, 0L),
    n_all_mutations = coalesce(n_all_mutations, 0L),
    n_mutated_candidate_genes = coalesce(n_mutated_candidate_genes, 0L),
    mutated_candidate_genes = coalesce(mutated_candidate_genes, ""),
    n_candidate_genes_considered = nrow(candidates),
    candidate_genes_considered = paste(
      sort(candidates$gene_symbol), collapse = ";"),
    n_tes_targeted_candidate_genes = nrow(tes_targeted_candidates),
    tes_targeted_candidate_genes = paste(
      sort(tes_targeted_candidates$gene_symbol), collapse = ";"),
    nonsynonymous_mutations_per_mbp =
      n_nonsynonymous_mutations / callable_mbp,
    all_mutations_per_mbp = n_all_mutations / callable_mbp,
    candidate_gene_count_definition = paste(
      "Distinct genes among all 30 final candidates with at least one",
      "configured nonsynonymous mutation observed in the mutation table"),
    candidate_profiling_caveat = paste(
      "TES uniformly targets only 8/30 candidates; TES zero counts do not",
      "establish wild-type status across all 30 genes"),
    burden_vaf_definition =
      "TES and WES mutation burdens require gt_AF >= 0.01") %>%
  arrange(dataset, sample_name)
stopifnot(nrow(plot_data) == 171L,
  nrow(plot_data) == n_distinct(plot_data$sample_uid),
  all(is.finite(plot_data$callable_mbp) & plot_data$callable_mbp > 0))

################################################################################
# @Project - TES Hair Follicle
# @Date - 05/01/2026
# @Author - Yoav Avi-Guy
# @Description - This file generates the working dataframes used in most of the
# following analysis
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "TES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# Loading work files
mutations_annotated <- readRDS(
  file.path(PROCESSED_PATH, "mutations_annotated.rds"))
gene_panel <- readRDS(
  file.path(PROCESSED_PATH, "gene_panel.rds"))

#-------------------------------------------------------------------------------
# Constants
panel_genes <- gene_panel %>% distinct(gene) %>% pull()
impact_levels <- c("HIGH","MODERATE","LOW","MODIFIER")
depth_cutoff <- 10

#-------------------------------------------------------------------------------
# Generating a dataframe with the following rules:
# 1) Single row per position, keeping preferably the row with the assigned gene 
# that matches the targeted gene panel, and then genes with the highest impact
# 2) Removing mitochondrial-assigned reads
# 3) Passing all MutectFilterCalls filters ("PASS")
# 4) Converting gt_DP to numeric
mutations_unique <- mutations_annotated %>%
  mutate(
    in_panel = SYMBOL %in% panel_genes,
    IMPACT   = factor(IMPACT, levels = impact_levels, ordered = TRUE),
    gt_DP = as.numeric(gt_DP),
    gt_AF = as.numeric(gt_AF)
  ) %>%
  group_by(sample_name, CHROM, POS, mutation) %>%
  mutate(any_in_panel = any(in_panel)) %>%
  filter(
    if_else(any_in_panel, in_panel, TRUE),
    CHROM %in% c(1:19, "X"),
    FILTER == "PASS",
    !grepl("ds", sample_name) # remove if you want to include downsampled samps
  ) %>%
  slice_min(order_by = IMPACT, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-in_panel, -any_in_panel)

#-------------------------------------------------------------------------------
# Generating a dataframe with mutations above the depth cutoff
mutations_unique_dp <- mutations_unique %>%
  filter(gt_DP >= depth_cutoff)

#-------------------------------------------------------------------------------
# Generating a dataframe with mutations within the targeted region
mutations_target <- inner_join(
  mutations_unique_dp, 
  gene_panel,
  by = join_by(
    CHROM, 
    between(POS, start, end))
)

#-------------------------------------------------------------------------------
# Generating a dataframes with SNVs according to the rules above
snvs <- mutations_unique %>%
  filter(VARIANT_CLASS == "SNV")

snvs_dp <- mutations_unique_dp %>%
  filter(VARIANT_CLASS == "SNV")

snvs_target <- mutations_target %>%
  filter(VARIANT_CLASS == "SNV")

#-------------------------------------------------------------------------------
# saving working dataframes
saveRDS(mutations_unique,
file.path(PROCESSED_PATH, "/mutations_unique.rds"))

saveRDS(mutations_unique_dp,
file.path(PROCESSED_PATH, "/mutations_unique_dp.rds"))

saveRDS(mutations_target,
file.path(PROCESSED_PATH, "/mutations_target.rds"))

saveRDS(snvs,
file.path(PROCESSED_PATH, "/snvs.rds"))

saveRDS(snvs_dp,
file.path(PROCESSED_PATH, "/snvs_dp.rds"))

saveRDS(snvs_target,
file.path(PROCESSED_PATH, "/snvs_target.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "07_working_dataframes_session.txt"))

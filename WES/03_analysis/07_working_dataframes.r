################################################################################
# @Project - WES Hair Follicle
# @Date - 02/03/2026
# @Author - Yoav Avi-Guy
# @Description - This file generates the working dataframes used in most of the
# following analysis
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "WES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# Libraries
library(stringr)

#-------------------------------------------------------------------------------
# Loading work files
mutations_annotated <- readRDS(
  file.path(PROCESSED_PATH, "mutations_annotated.rds"))

#-------------------------------------------------------------------------------
# Constants
impact_levels <- c("HIGH","MODERATE","LOW","MODIFIER")
depth_cutoff <- 8

#-------------------------------------------------------------------------------
# Generating a dataframe with the following rules:
# 1) Single row per position, keeping preferably genes with the highest impact
# 2) Removing mitochondrial-assigned reads
# 3) Passing all MutectFilterCalls filters ("PASS")
# 4) Converting gt_DP to numeric
mutations_unique <- mutations_annotated %>%
  mutate(
    IMPACT   = factor(IMPACT, levels = impact_levels, ordered = TRUE),
    gt_DP = as.numeric(gt_DP),
    gt_AF = as.numeric(gt_AF)
  ) %>%
  group_by(sample_name, CHROM, POS, mutation) %>%
  filter(
    CHROM %in% c(seq(1:19), "X"),
    FILTER == "PASS"
  ) %>%
  slice_min(order_by = IMPACT, with_ties = FALSE) %>%
  ungroup()

#-------------------------------------------------------------------------------
# Filtering for the depth cutoff
mutations_unique_dp <- mutations_unique %>%
  filter(gt_DP >= depth_cutoff)

#-------------------------------------------------------------------------------
# Generating a dataframes with SNVs according to the rules above
snvs <- mutations_unique %>%
  filter(VARIANT_CLASS == "SNV")

snvs_dp <- mutations_unique_dp %>%
  filter(VARIANT_CLASS == "SNV")

#-------------------------------------------------------------------------------
# saving working dataframes
saveRDS(mutations_unique,
file.path(PROCESSED_PATH, "/mutations_unique.rds"))

saveRDS(mutations_unique_dp,
file.path(PROCESSED_PATH, "/mutations_unique_dp.rds"))

saveRDS(snvs,
file.path(PROCESSED_PATH, "/snvs.rds"))

saveRDS(snvs_dp,
file.path(PROCESSED_PATH, "/snvs_dp.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "07_working_dataframes_session.txt"))

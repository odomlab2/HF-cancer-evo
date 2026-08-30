################################################################################
# @Project - WES Hair Follicle
# @Date - 04/03/2026
# @Author - Yoav Avi-Guy
# @Description - This file generates dndscv objects for gene selection analysis
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "WES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# Libraries
library(purrr)
library(stringr)
library(dndscv)

#-------------------------------------------------------------------------------
# Paths
REF_PATH <- Sys.getenv("HF_SCC_REFCDS")
if (!nzchar(REF_PATH)) stop("Set HF_SCC_REFCDS to the reference CDS RDA file.")

#-------------------------------------------------------------------------------
# Loading work files
mutations_unique_dp <- readRDS(file.path(PROCESSED_PATH, "mutations_unique_dp.rds"))
outliers <- readRDS(file.path(PROCESSED_PATH, "technical_outliers.rds"))

#-------------------------------------------------------------------------------
# Building the input dataframe
dnds_df <- mutations_unique_dp %>%
  filter(
    last_cycle == TRUE,
    !sample_name %in% outliers$sample_name,
    gt_AF >= 0.01,
    tissue == "Skin",
    treatment == "DT") %>%
  dplyr::select(sample_name, CHROM, POS, REF, ALT, tissue, treatment, category)

#-------------------------------------------------------------------------------
# Generating a dndscv object for tissue
dnds_results <- dnds_df %>%
  group_split(.keep = TRUE) %>%
  map(~ dndscv(
    .x,
    refdb = REF_PATH,
    max_muts_per_gene_per_sample = Inf,
    cv = NULL))

#-------------------------------------------------------------------------------
# Generating a dndscv object for category comparison in each tissue
dnds_category_results <- dnds_df %>%
  filter(!category %in% c("Acetone", "SCC")) %>%
  split(.$category) %>%
  map(~ dndscv(
    .x,
    refdb = REF_PATH,
    max_muts_per_gene_per_sample = Inf,
    cv = NULL))

#-------------------------------------------------------------------------------
# Saving dndscv objects
saveRDS(dnds_results,
  file.path(PROCESSED_PATH, "/dnds_results.rds"))

saveRDS(dnds_category_results,
  file.path(PROCESSED_PATH, "/dnds_category_results.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "09_dnds_analysis_session.txt"))

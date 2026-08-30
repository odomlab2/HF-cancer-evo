################################################################################
# @Project - TES Hair Follicle
# @Date - 11/02/2026
# @Author - Yoav Avi-Guy
# @Description - This file generates dndscv objects for gene selection analysis
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "TES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# Libraries
  library(purrr)
  library(dndscv)

#-------------------------------------------------------------------------------
# Paths
REF_PATH <- Sys.getenv("HF_SCC_REFCDS")
if (!nzchar(REF_PATH)) stop("Set HF_SCC_REFCDS to the reference CDS RDA file.")

#-------------------------------------------------------------------------------
# Loading work files
mutations_target <- readRDS(file.path(PROCESSED_PATH, "mutations_target.rds"))
mutations_unique_dp <- readRDS(file.path(PROCESSED_PATH, "mutations_unique_dp.rds"))
outliers <- readRDS(file.path(PROCESSED_PATH, "technical_outliers.rds"))
gene_panel <- readRDS(file.path(PROCESSED_PATH, "gene_panel.rds"))

#-------------------------------------------------------------------------------
# Building the input dataframe, using only the DT-treated mouse, as the acetone
# treatment did not produce enough mutations for the analysis
# dnds_df <- mutations_target %>%
#   filter(
#     !sample_name %in% outliers$sample_name,
#     treatment == "DT") %>%
#   dplyr::select(sample_name, CHROM, POS, REF, ALT, tissue, treatment, time, category)

# dnds_df <- dnds_df %>%
#     filter(!sample_name %in% outliers$sample_name) %>%
#     mutate(tissue = case_when(
#         grepl("ds", sample_name) ~ "Skin_ds",
#         TRUE ~ tissue))

dnds_all <- mutations_unique_dp %>%
  filter(
    !sample_name %in% outliers$sample_name,
    gt_AF >= 0.01,
    treatment == "DT",
    tissue == "Skin") %>%
  dplyr::select(sample_name, CHROM, POS, REF, ALT, tissue, treatment, time, category)

#-------------------------------------------------------------------------------
# Generating a dndscv object for the skin of the DT-treated mouse
dnds_results <- dndscv(
  dnds_all,
  refdb = REF_PATH,
  gene_list = unique(gene_panel$gene),
  max_muts_per_gene_per_sample = Inf,
  cv = NULL
)

#-------------------------------------------------------------------------------
# Generating a dndscv object for time comparison in hair follicle samples
# dnds_time_results <- dnds_all %>%
#   filter(tissue == "Hair follicle") %>%
#   group_by(time) %>%
#   group_split(.keep = TRUE) %>%
#   set_names(dnds_all %>% distinct(time) %>% pull(time)) %>%
#   map(~ dndscv(
#     .x,
#     refdb = REF_PATH,
#     gene_list = unique(gene_panel$gene),
#     max_muts_per_gene_per_sample = Inf),
#     cv = NULL)

#-------------------------------------------------------------------------------
# Generating a dndscv object for category comparison in the skin
# dnds_category_df <- dnds_all %>%
#   filter(category != "Acetone") %>%
#   group_by(category)

# category_names <- dnds_category_df %>%
#   group_keys() %>%
#   pull(category)

# dnds_category_results <- dnds_category_df %>%
#   group_split(.keep = TRUE) %>%
#   set_names(category_names) %>%
#   map(~ dndscv(
#     .x,
#     refdb = REF_PATH,
#     gene_list = unique(gene_panel$gene),
#     max_muts_per_gene_per_sample = Inf,
#     cv = NULL
#   ))

#-------------------------------------------------------------------------------
# Saving dndscv objects
saveRDS(dnds_results,
  file.path(PROCESSED_PATH, "/dnds_results.rds"))

# saveRDS(dnds_time_results,
#   file.path(PROCESSED_PATH, "/dnds_time_results.rds"))

# saveRDS(dnds_category_results,
#   file.path(PROCESSED_PATH, "/dnds_category_results.rds"))


#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "10_dnds_analysis_session.txt"))

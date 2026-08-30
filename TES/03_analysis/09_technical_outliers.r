################################################################################
# @Project - TES Hair Follicle
# @Date - 07/01/2026
# @Author - Yoav Avi-Guy
# @Description - This file identifies hypermutator samples
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "TES/01_source/00_config.r"))
source(file.path(SOURCE_PATH, "01_plotting.r"))
source(file.path(HF_SCC_ROOT, "R", "shared", "sample_manifest.R"))
source(file.path(SOURCE_PATH, "02_sample_manifest.r"))
source(file.path(HF_SCC_ROOT, "R", "shared", "technical_outliers.R"))

#-------------------------------------------------------------------------------
# Loading work files
mutations_unique_dp <- readRDS(
  file.path(PROCESSED_PATH, "mutations_unique_dp.rds"))
sample_manifest <- tes_sample_manifest(file.path(INPUT_PATH, "samples.txt"))
outlier_sample_ids <- select_manifest_stage(
  sample_manifest, "technical_outliers")$sample_name
assert_sample_ids(outlier_sample_ids, unique(mutations_unique_dp$sample_name),
  label = "TES technical-outlier universe")

#-------------------------------------------------------------------------------
# Constants
hypermutator_cols <- c(
  "FALSE" = "black",
  "TRUE" = "red")

#-------------------------------------------------------------------------------
# Generating mutations per sample table
technical_outlier_results <- identify_technical_outliers(
  mutations_unique_dp %>%
    count(sample_name, tissue, treatment, callable_mbp, name = "mutations"))
hypermutator_thresholds <- technical_outlier_results$hypermutator_thresholds
hypermutator_samples <- technical_outlier_results$hypermutator_samples
technical_outliers <- technical_outlier_results$technical_outliers

#-------------------------------------------------------------------------------
# Visualisation
ggplot(
  hypermutator_samples,
  aes(x = tissue, y = mutations_per_mbp, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = col_palette$tissue) +
  geom_jitter(
    aes(color = outlier),
    width = 0.2,
    size = 1.5,
    alpha = 0.8) +
  scale_color_manual(values = hypermutator_cols) +
  labs(
    x = "",
    y = "# of mutations per Mb") +
  my_theme +
  theme(legend.title = element_blank())

ggplot(
  hypermutator_samples,
  aes(x = tissue, y = log_mutations_per_mbp, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = col_palette$tissue) +
  geom_jitter(
    aes(color = outlier),
    width = 0.2,
    size = 1.5,
    alpha = 0.8) +
  scale_color_manual(values = hypermutator_cols) +
  labs(
    x = "",
    y = "log((# of mutations + 0.05) per Mb)") +
  my_theme +
  theme(legend.title = element_blank())

#-------------------------------------------------------------------------------
# Saving hypermutator samples
saveRDS(technical_outliers,
file.path(PROCESSED_PATH, "/technical_outliers.rds"))


#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "09_technical_outliers_session.txt"))

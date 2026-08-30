################################################################################
# @Project - WES Hair Follicle
# @Date - 02/03/2026
# @Author - Yoav Avi-Guy
# @Description - This file identifies publication-authoritative technical outliers
################################################################################
#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "WES/01_source/00_config.r"))
source(file.path(SOURCE_PATH, "01_plotting.r"))
source(file.path(HF_SCC_ROOT, "R", "shared", "sample_manifest.R"))
source(file.path(SOURCE_PATH, "03_sample_manifest.r"))
source(file.path(HF_SCC_ROOT, "R", "shared", "technical_outliers.R"))

#-------------------------------------------------------------------------------
# Loading work files
mutations_unique_dp <- readRDS(
  file.path(PROCESSED_PATH, "mutations_unique_dp.rds"))
sample_metadata <- readRDS(file.path(METADATA_PATH, "sample_metadata.rds"))
callable_summary <- read.delim(file.path(PROCESSED_PATH, "callable_summary.tsv"))
sample_manifest <- wes_sample_manifest(
  file.path(INPUT_PATH, "samples.txt"),
  file.path(SOURCE_PATH, "wes_matched_normal_controls.csv"))

#-------------------------------------------------------------------------------
# Constants
hypermutator_cols <- c(
  "FALSE" = "black",
  "TRUE" = "red")

#-------------------------------------------------------------------------------
# Generating mutations per sample table
eligible_sample_ids <- select_manifest_stage(
  sample_manifest, "technical_outliers")$sample_name
assert_sample_ids(eligible_sample_ids, sample_metadata$sample_name,
  label = "WES technical-outlier metadata universe")
metadata_rows <- match(eligible_sample_ids, sample_metadata$sample_name)
eligible_samples <- sample_metadata[metadata_rows, ] %>%
  select(sample_name, tissue, treatment) %>%
  left_join(
    callable_summary %>% transmute(
      sample_name,
      callable_mbp = overall_callable_10x / 1e6),
    by = "sample_name") %>%
  left_join(
    mutations_unique_dp %>% count(sample_name, name = "mutations"),
    by = "sample_name")
technical_outlier_results <- identify_technical_outliers(
  eligible_samples)
hypermutator_thresholds <- technical_outlier_results$hypermutator_thresholds
hypermutator_samples <- technical_outlier_results$hypermutator_samples
technical_outliers <- technical_outlier_results$technical_outliers
valid_hypermutator_samples <- hypermutator_samples %>%
  filter(!is.na(callable_mbp), callable_mbp > 0)

#-------------------------------------------------------------------------------
# Visualisation
invisible(ggplot(
  valid_hypermutator_samples,
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
  theme(legend.title = element_blank()))

invisible(ggplot(
  valid_hypermutator_samples,
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
  theme(legend.title = element_blank()))

#-------------------------------------------------------------------------------
# Saving hypermutator samples
saveRDS(technical_outliers,
file.path(PROCESSED_PATH, "/technical_outliers.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "08_technical_outliers_session.txt"))

################################################################################
# @Project - WES Hair Follicle
# @Date - 04/03/2026
# @Author - Yoav Avi-Guy
# @Description - This file generates dataframes for tracing clones between the 
# hair follicle and skin samples
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "WES/01_source/00_config.r"))
source(file.path(SOURCE_PATH, "02_grid_adjacency.r"))

#-------------------------------------------------------------------------------
# Loading work files
mutations_unique_dp <- readRDS(file.path(PROCESSED_PATH, "mutations_unique_dp.rds"))
outliers <- readRDS(file.path(PROCESSED_PATH, "technical_outliers.rds"))
sample_metadata <- readRDS(file.path(METADATA_PATH, "sample_metadata.rds"))

#-------------------------------------------------------------------------------
# Generating a unique mutation ID dataframe
mutation_id <- mutations_unique_dp %>%
  filter(
    last_cycle == TRUE,
    !sample_name %in% outliers$sample_name,
    treatment != "Acetone",
    gt_AF >= 0.01) %>%
  transmute(
    sample_name,
    mut_id = paste(mouse, CHROM, POS, REF, ALT, sep = ":")) %>%
  distinct()

hf_meta <- sample_metadata %>% 
  filter(
    last_cycle == TRUE,
    tissue == "Hair follicle",
    treatment != "Acetone",
    !sample_name %in% outliers$sample_name)

sk_meta <- sample_metadata %>% 
  filter(
    last_cycle == TRUE,
    tissue == "Skin",
    treatment != "Acetone",
    !sample_name %in% outliers$sample_name)

hf_muts <- mutation_id %>% semi_join(hf_meta, by = "sample_name")
sk_muts <- mutation_id %>% semi_join(sk_meta, by = "sample_name")

#-------------------------------------------------------------------------------
# Comparing the number of mutations shared between hair follicle and skin 
# samples and the distance between them on the grid
shared_hf_skin <- inner_join(
  hf_muts %>% dplyr::rename(hf_sample_name = sample_name),
  sk_muts %>% dplyr::rename(sk_sample_name = sample_name),
  by = "mut_id",
  relationship = "many-to-many") %>%
  distinct(hf_sample_name, sk_sample_name, mut_id) %>%
  count(hf_sample_name, sk_sample_name, name = "n_shared")

shared_muts_hf_skin <- crossing(
  hf_meta %>%
    dplyr::rename(hf_sample_name = sample_name) %>%
    rename_with(~ paste0("hf_", .x), -hf_sample_name),
  sk_meta %>%
    dplyr::rename(sk_sample_name = sample_name) %>%
    rename_with(~ paste0("sk_", .x), -sk_sample_name)) %>%
  filter(hf_mouse == sk_mouse) %>%
  left_join(shared_hf_skin, by = c("hf_sample_name", "sk_sample_name")) %>%
  mutate(
    n_shared   = coalesce(n_shared, 0L),
    distance   = calculate_distance(hf_grid_x, hf_grid_y, sk_grid_x, sk_grid_y),
    same_mouse = hf_mouse == sk_mouse,
    adjacency  = adjacency_class(distance, same_mouse))

#-------------------------------------------------------------------------------
# Saving dndscv objects
saveRDS(shared_muts_hf_skin,
  file.path(PROCESSED_PATH, "/shared_muts_hf_skin.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "11_grid_tracking_session.txt"))

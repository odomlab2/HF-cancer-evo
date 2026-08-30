################################################################################
# @Project - TES Hair Follicle
# @Date - 11/02/2026
# @Author - Yoav Avi-Guy
# @Description - This file generates dataframes for tracing clones between the 
# hair follicle and skin samples or between the weeks
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "TES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# Loading work files
mutations_target <- readRDS(file.path(PROCESSED_PATH, "mutations_target.rds"))
mutations_unique_dp <- readRDS(file.path(PROCESSED_PATH, "mutations_unique_dp.rds"))
outliers <- readRDS(file.path(PROCESSED_PATH, "technical_outliers.rds"))
sample_metadata <- readRDS(file.path(METADATA_PATH, "sample_metadata.rds"))

#-------------------------------------------------------------------------------
# Helper functions
calculate_distance <- function(x1, y1, x2, y2) pmax(abs(x1 - x2), abs(y1 - y2))

adjacency_class <- function(d) if_else(d <= 1, "Adjacent", "Distal")

## Separating downsampled samples
# mutations_target <- mutations_target %>%
#     filter(!sample_name %in% outliers$sample_name) %>%
#     mutate(tissue = case_when(
#         grepl("ds", sample_name) ~ "Skin_ds",
#         TRUE ~ tissue))

# mutations_unique_dp <- mutations_unique_dp %>%
#     filter(!sample_name %in% outliers$sample_name) %>%
#     mutate(tissue = case_when(
#         grepl("ds", sample_name) ~ "Skin_ds",
#         TRUE ~ tissue))

# sample_metadata <- sample_metadata %>%
#     filter(!sample_name %in% outliers$sample_name) %>%
#     mutate(tissue = case_when(
#         grepl("ds", sample_name) ~ "Skin_ds",
#         TRUE ~ tissue))

#-------------------------------------------------------------------------------
# Generating a unique mutation ID dataframe
mutation_id <- mutations_unique_dp %>%
  filter(
    !sample_name %in% outliers$sample_name,
    treatment != "Acetone",
    gt_AF >= 0.01) %>%
  transmute(
    sample_name,
    mut_id = paste(mouse, CHROM, POS, REF, ALT, sep = ":")) %>%
  distinct()

# Unique mutation ID dataframe for on-target mutations
# mutation_id_ot <- mutations_target %>%
#   filter(
#     !sample_name %in% outliers$sample_name,
#     gt_AF >= 0.05) %>%
#   transmute(
#     sample_name,
#     mut_id = paste(mouse, CHROM, POS, mutation, sep = ":")) %>%
#   distinct()

hf_meta <- sample_metadata %>% 
  filter(
    tissue == "Hair follicle",
    treatment != "Acetone",
    !sample_name %in% outliers$sample_name)

sk_meta <- sample_metadata %>% 
  filter(
    tissue == "Skin",
    treatment != "Acetone",
    !sample_name %in% outliers$sample_name)

# sk_ds_meta <- sample_metadata %>% 
#   filter(
#     tissue == "Skin_ds",
#     !sample_name %in% outliers$sample_name)

hf_muts <- mutation_id %>% semi_join(hf_meta, by = "sample_name")
sk_muts <- mutation_id %>% semi_join(sk_meta, by = "sample_name")
# sk_ds_muts <- mutation_id %>% semi_join(sk_ds_meta, by = "sample_name")

# hf_muts_ot <- mutation_id_ot %>% semi_join(hf_meta, by = "sample_name")
# sk_muts_ot <- mutation_id_ot %>% semi_join(sk_meta, by = "sample_name")
# sk_ds_muts_ot <- mutation_id_ot %>% semi_join(sk_ds_meta, by = "sample_name")

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
  left_join(shared_hf_skin, by = c("hf_sample_name", "sk_sample_name")) %>%
  mutate(
    n_shared   = coalesce(n_shared, 0L),
    distance   = calculate_distance(hf_grid_x, hf_grid_y, sk_grid_x, sk_grid_y),
    adjacency  = adjacency_class(distance))

# Repeating for on-target mutations
# shared_hf_skin_ot <- inner_join(
#   hf_muts_ot %>% rename(hf_sample_name = sample_name),
#   sk_muts_ot %>% rename(sk_sample_name = sample_name),
#   by = "mut_id",
#   relationship = "many-to-many") %>%
#   distinct(hf_sample_name, sk_sample_name, mut_id) %>%
#   count(hf_sample_name, sk_sample_name, name = "n_shared")

# shared_muts_hf_skin_ot <- crossing(
#   hf_meta %>%
#     rename(hf_sample_name = sample_name) %>%
#     rename_with(~ paste0("hf_", .x), -hf_sample_name),
#   sk_meta %>%
#     rename(sk_sample_name = sample_name) %>%
#     rename_with(~ paste0("sk_", .x), -sk_sample_name)) %>%
#   left_join(shared_hf_skin_ot, by = c("hf_sample_name", "sk_sample_name")) %>%
#   mutate(
#     n_shared   = coalesce(n_shared, 0L),
#     distance   = calculate_distance(hf_grid_x, hf_grid_y, sk_grid_x, sk_grid_y),
#     adjacency  = adjacency_class(distance))

#-------------------------------------------------------------------------------
# Comparing the number of mutations shared between hair follicle samples over 
# different weeks, and the distance between them on the grid
shared_hf_hf <- inner_join(
  hf_muts %>% dplyr::rename(hf_ref = sample_name),
  hf_muts %>% dplyr::rename(hf_cmp = sample_name),
  by = "mut_id",
  relationship = "many-to-many") %>%
  filter(hf_ref != hf_cmp) %>%
  distinct(hf_ref, hf_cmp, mut_id) %>%
  count(hf_ref, hf_cmp, name = "n_shared")

shared_muts_hf_time <- crossing(
  hf_meta %>% 
    dplyr::rename(hf_ref = sample_name) %>% 
    rename_with(~ paste0("ref_", .x), -hf_ref),
  hf_meta %>% 
    dplyr::rename(hf_cmp = sample_name) %>% 
    rename_with(~ paste0("cmp_", .x), -hf_cmp)) %>%
  filter(hf_ref != hf_cmp, ref_time != cmp_time) %>%
  left_join(shared_hf_hf, by = c("hf_ref", "hf_cmp")) %>%
  mutate(
    n_shared   = coalesce(n_shared, 0L),
    distance   = calculate_distance(
      ref_grid_x, ref_grid_y, cmp_grid_x, cmp_grid_y),
    adjacency  = adjacency_class(distance))

# Repeating for mutations on-target
# shared_hf_hf_ot <- inner_join(
#   hf_muts_ot %>% rename(hf_ref = sample_name),
#   hf_muts_ot %>% rename(hf_cmp = sample_name),
#   by = "mut_id",
#   relationship = "many-to-many") %>%
#   filter(hf_ref != hf_cmp) %>%
#   distinct(hf_ref, hf_cmp, mut_id) %>%
#   count(hf_ref, hf_cmp, name = "n_shared")

# shared_muts_hf_time_ot <- crossing(
#   hf_meta %>% 
#     rename(hf_ref = sample_name) %>% 
#     rename_with(~ paste0("ref_", .x), -hf_ref),
#   hf_meta %>% 
#     rename(hf_cmp = sample_name) %>% 
#     rename_with(~ paste0("cmp_", .x), -hf_cmp)) %>%
#   filter(hf_ref != hf_cmp, ref_time != cmp_time) %>%
#   left_join(shared_hf_hf_ot, by = c("hf_ref", "hf_cmp")) %>%
#   mutate(
#     n_shared   = coalesce(n_shared, 0L),
#     distance   = calculate_distance(
#       ref_grid_x, ref_grid_y, cmp_grid_x, cmp_grid_y),
#     adjacency  = adjacency_class(distance))

#-------------------------------------------------------------------------------
# Saving dndscv objects
saveRDS(shared_muts_hf_skin,
  file.path(PROCESSED_PATH, "/shared_muts_hf_skin.rds"))

# saveRDS(shared_muts_hf_skin_ot,
#   file.path(PROCESSED_PATH, "/shared_muts_hf_skin_ot.rds"))

saveRDS(shared_muts_hf_time,
  file.path(PROCESSED_PATH, "/shared_muts_hf_time.rds"))

# saveRDS(shared_muts_hf_time_ot,
#   file.path(PROCESSED_PATH, "/shared_muts_hf_time_ot.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "11_grid_tracking_session.txt"))

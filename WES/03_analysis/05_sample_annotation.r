################################################################################
# @Project - WES Hair Follicle
# @Date - 02/03/2025
# @Author - Yoav Avi-Guy
# @Description - This file is responsible for the sample metadata preprocessing 
# for later analysis stratification
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "WES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# Libraries
library(stringr)

#-------------------------------------------------------------------------------
# Loading work files
mutations_df <- readRDS(file.path(
  PROCESSED_PATH, 
  "mutations_df.rds"))
category_annotation <- read.csv(file.path(
  METADATA_PATH, 
  "WES_annotation.csv"))
callable_summary <- read.table(file.path(
  PROCESSED_PATH, 
  "callable_summary.tsv"), 
  header = TRUE)

#-------------------------------------------------------------------------------
# Constants
treated_mice <- c(3, 4, 5)
acetone_mice <- c(9)

#-------------------------------------------------------------------------------
# Adding visual category annotation to the samples
mutations_annotated <- mutations_df %>%
  left_join(
    category_annotation %>% select(sample_name, grid, category),
    by = c("sample_name", "grid")
  ) %>%
  mutate(
    category = case_when(
      mouse %in% acetone_mice ~ "Acetone",
      mouse %in% treated_mice ~ category
    )
  )

#-------------------------------------------------------------------------------
# Adding callable bases per sample
mutations_annotated <- mutations_annotated %>%
  left_join(
    callable_summary %>%
    transmute(
      sample_name,
      callable_mbp = overall_callable_10x / 1e6),
      by = "sample_name")

#-------------------------------------------------------------------------------
# Adding manual annotations for later analysis
mutations_annotated <- mutations_annotated %>%
  mutate(
    treatment = case_when(
      mouse %in% acetone_mice ~ "Acetone",
      mouse %in% treated_mice ~ "DT"
    ),
    grid_y = match(substr(grid, 1, 1), letters[1:7]),
    grid_x = as.integer(substr(grid, 2, 2)),
    mutation = paste0(REF, ">", ALT),
    mutation = recode(mutation,
      "A>C" = "T>G",
      "A>G" = "T>C",
      "A>T" = "T>A",
      "G>A" = "C>T",
      "G>C" = "C>G",
      "G>T" = "C>A"))

#-------------------------------------------------------------------------------
# Adding a counter of the sequencing time
mutations_annotated <- mutations_annotated %>%
  mutate(
    seq_cycle = case_when(
      str_detect(sample_name, "-t2$") ~ "t2",
      str_detect(sample_name, "-t3$") ~ "t3",
      TRUE ~ "t1"),
    sample_base = str_remove(sample_name, "-t[23]$")) %>%
  group_by(sample_base) %>%
  mutate(last_cycle = seq_cycle == max(seq_cycle, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-sample_base)

#-------------------------------------------------------------------------------
# Reordering columns for easier downstream analysis
mutations_annotated_ordered <- mutations_annotated %>%
  relocate(all_of(c(
    "sample_name", "tissue", "grid", "treatment", "category", "CHROM",
    "POS", "REF", "ALT", "mutation"
  )))

#-------------------------------------------------------------------------------
# Generating a sample metadata dataframe
sample_metadata <- mutations_annotated_ordered %>%
  select(
    sample_name, 
    tissue, 
    grid,
    treatment, 
    category, 
    mouse,
    callable_mbp, 
    grid_x,
    grid_y,
    seq_cycle,
    last_cycle) %>%
  distinct(sample_name, .keep_all = TRUE)

#-------------------------------------------------------------------------------
# saving annotated dataframe and sample metadata
saveRDS(mutations_annotated_ordered,
file.path(PROCESSED_PATH, "/mutations_annotated.rds"))

saveRDS(sample_metadata,
file.path(METADATA_PATH, "/sample_metadata.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "05_sample_annotation_session.txt"))

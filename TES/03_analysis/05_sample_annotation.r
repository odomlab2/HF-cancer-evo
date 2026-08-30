################################################################################
# @Project - TES Hair Follicle
# @Date - 29/12/2025
# @Author - Yoav Avi-Guy
# @Description - This file is responsible for the sample metadata preprocessing 
# for later analysis stratification
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "TES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# Loading work files
mutations_df <- readRDS(file.path(
  PROCESSED_PATH, 
  "mutations_df.rds"))
category_annotation <- read.csv(file.path(
  METADATA_PATH, 
  "Mouse21_full_annotation.csv"))
callable_summary <- read.table(file.path(
  PROCESSED_PATH, 
  "callable_summary.tsv"), 
  header = TRUE)
seq_qc_technical <- readRDS(file.path(
  PROCESSED_PATH,
  "seq_qc_technical.rds"))

#-------------------------------------------------------------------------------
# Constants
treated_mice <- c(21)
acetone_mice <- c(19)

#-------------------------------------------------------------------------------
# Adding visual category annotation to the samples
mutations_annotated <- mutations_df %>%
  mutate(
    time = case_when(
      tissue == "Skin" ~ "Skin 19",
      TRUE ~ time)) %>%
  left_join(
    category_annotation %>% dplyr::select(time, grid, category),
    by = c("time", "grid")
  ) %>%
  mutate(
    category = case_when(
      mouse %in% acetone_mice ~ "Acetone",
      mouse %in% treated_mice ~ category
    )
  )

#-------------------------------------------------------------------------------
# Adding coverage and deduplicated read count per sample
mutations_annotated <- mutations_annotated %>%
  left_join(
    seq_qc_technical %>%
    transmute(
      sample_name,
      coverage,
      reads = non_duplicate_mapped_reads),
      by = "sample_name")

#-------------------------------------------------------------------------------
# Adding callable bases per sample
mutations_annotated <- mutations_annotated %>%
  left_join(
    callable_summary %>%
    transmute(
      sample_name,
      ot_callable_mbp = callable_bases_10x / 1e6,
      ot_callable_pct = pct_callable_10x,
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
# Reordering columns for easier downstream analysis
mutations_annotated_ordered <- mutations_annotated %>%
  relocate(all_of(c(
    "sample_name", "tissue", "grid", "time", "treatment", "category", "CHROM",
    "POS", "REF", "ALT", "mutation"
  )))

#-------------------------------------------------------------------------------
# Manual annotation of samples with names that do not fit the strcture
mutations_annotated_ordered <- mutations_annotated_ordered %>%
  mutate(
    sample_base = sub("(_ds\\d+)$", "", sample_name)
  ) %>%
  mutate(
    tissue = case_when(
      sample_base %in% c("skin21-pp-bc-w19", "skin21-pp-de-w19") ~ "Skin",
      TRUE ~ tissue
    ),
    mouse = as.integer(case_when(
      sample_base %in% c("skin21-pp-bc-w19", "skin21-pp-de-w19") ~ 21,
      TRUE ~ mouse
    )),
    grid = case_when(
      sample_base == "skin21-pp-bc-w19" ~ "b6c7",
      sample_base == "skin21-pp-de-w19" ~ "de",
      TRUE ~ grid
    ),
    treatment = case_when(
      sample_base %in% c("skin21-pp-bc-w19", "skin21-pp-de-w19") ~ "DT",
      TRUE ~ treatment
    ),
    category = case_when(
      sample_base %in% c("skin21-pp-bc-w19", "skin21-pp-de-w19") ~ "SCC",
      TRUE ~ category
    ),
    time = case_when(
      sample_base %in% c("skin21-pp-bc-w19", "skin21-pp-de-w19") ~ "Skin 19",
      TRUE ~ time
    )
  ) %>%
  dplyr::select(-sample_base)

#-------------------------------------------------------------------------------
# Generating a sample metadata dataframe
sample_metadata <- mutations_annotated_ordered %>%
  filter(!str_detect(sample_name, "ds")) %>%
  dplyr::select(
    sample_name, 
    tissue, 
    grid, 
    time, 
    treatment, 
    category, 
    mouse,
    ot_callable_mbp, 
    callable_mbp, 
    coverage, 
    reads,
    grid_x,
    grid_y) %>%
  distinct(sample_name, .keep_all = TRUE)

#-------------------------------------------------------------------------------
# Annotating technical QC after metadata exists keeps the dependency acyclic
overlap <- intersect(
  names(seq_qc_technical),
  setdiff(names(sample_metadata), "sample_name"))
seq_qc_summary <- seq_qc_technical %>%
  dplyr::select(-all_of(overlap)) %>%
  left_join(sample_metadata, by = "sample_name")

#-------------------------------------------------------------------------------
# Saving annotated mutations, sample metadata and downstream annotated QC
saveRDS(mutations_annotated_ordered,
file.path(PROCESSED_PATH, "/mutations_annotated.rds"))

saveRDS(sample_metadata,
  file.path(METADATA_PATH, "/sample_metadata.rds"))

saveRDS(seq_qc_summary,
  file.path(PROCESSED_PATH, "/seq_qc_summary.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "05_sample_annotation_session.txt"))

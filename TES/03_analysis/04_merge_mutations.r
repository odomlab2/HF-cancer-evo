################################################################################
# @Project - TES Hair Follicle
# @Date - 19/01/2026
# @Author - Yoav Avi-Guy
# @Description - This file merges the Ensembl VEP output files (annotated 
# mutations) to a single dataframe and adds some basic annotations
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "TES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# PATHS
VCF_PATH <- file.path(PROCESSED_PATH, "/mutation_annotation")

#-------------------------------------------------------------------------------
# Libraries
library(vcfR)
library(stringr)

#-------------------------------------------------------------------------------
# Helpers
## Listing vcf files to merge
vcf_files <- list.files(
  VCF_PATH, 
  pattern = "_annot\\.vcf\\.gz$", 
  full.names = TRUE)

## Parsing INFO column to separate fields
parse_info <- function(info) {
  p <- strsplit(info, ";", fixed = TRUE)
  keys <- unique(unlist(lapply(p, function(x) sub("=.*", "", x))))
  m <- matrix(
    NA_character_, length(p), length(keys), dimnames = list(NULL, keys))
  for (i in seq_along(p)) {
    k <- sub("=.*", "", p[[i]])
    v <- ifelse(grepl("=", p[[i]]), sub("^[^=]+=", "", p[[i]]), "TRUE")
    m[i, match(k, keys)] <- v
  }
  as_tibble(as.data.frame(m, check.names = FALSE, stringsAsFactors = FALSE))
}

## A function to read each VCF file
read_one_vcf <- function(f) {
  v <- read.vcfR(f, verbose = FALSE)
  fix <- as.data.frame(v@fix, stringsAsFactors = FALSE)
  gt  <- v@gt
  fmt_fields <- strsplit(gt[, "FORMAT"][1], ":", fixed = TRUE)[[1]]
  samp_col <- setdiff(colnames(gt), "FORMAT")[1]
  gt_df <- tibble(value = gt[, samp_col]) %>%
    separate(
      value, into = fmt_fields, sep = ":", fill = "right", remove = TRUE) %>%
    rename_with(~ paste0("gt_", .x))

  info_df <- parse_info(fix$INFO)
  csq_fields <- strsplit(
    sub(
      '.*Format: ',
      '', 
      sub('">.*', 
      '', 
      grep("^##INFO=<ID=CSQ", v@meta, value = TRUE)[1])),
      "\\|")[[1]]

  tibble(
    sample_name = sub("_annot\\.vcf\\.gz$", "", basename(f)),
    CHROM = fix$CHROM,
    POS   = as.integer(fix$POS),
    ID    = fix$ID,
    REF   = fix$REF,
    ALT   = fix$ALT,
    QUAL  = suppressWarnings(as.numeric(fix$QUAL)),
    FILTER= fix$FILTER
  ) %>%
    dplyr::bind_cols(info_df, gt_df) %>%
    dplyr::mutate(CSQ = strsplit(CSQ, ",", fixed = TRUE)) %>%
    tidyr::unnest(CSQ) %>%
    tidyr::separate(
      CSQ, 
      into = csq_fields,
      sep = "\\|", 
      fill = "right", 
      extra = "drop")
}

#-------------------------------------------------------------------------------
# Merging the VCF files to a single dataframe
merged_df <- bind_rows(lapply(vcf_files, read_one_vcf))

#-------------------------------------------------------------------------------
# Adding annotation fields for downstream analysis
annotated_df <- merged_df %>%
  tidyr::extract(
    sample_name,
    into   = c("tissue","mouse","grid","weeknum"),
    regex  = "^([A-Za-z]+)(\\d{1,2})-([A-Za-z]\\d{1,2})-w(\\d{1,2})(?:_ds\\d+)?$",
    remove = FALSE
  ) %>%
  mutate(
    tissue = case_when(
      tissue == "skin" ~ "Skin",
      tissue == "hairfollicle" ~ "Hair follicle",
      TRUE ~ tissue
    ),
    mouse = as.integer(mouse),
    time  = case_when(
      tissue == "Hair follicle" ~ paste("Week", as.integer(weeknum)),
      tissue == "Skin" ~ "Skin 19"
  )) %>%
  select(-weeknum)

#-------------------------------------------------------------------------------
# saving annotated dataframe
saveRDS(annotated_df,
file.path(PROCESSED_PATH, "/mutations_df.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "04_merge_mutations_session.txt"))

# 12/05/2026
# Author: Yoav Avi-Guy
# Purpose: Figure 4 panels - candidate genes for tumour initiation and SCC
# predisposition, their skin VAF, and human cSCC recurrence.
#-------------------------------------------------------------------------------
# Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)

# Colours and theme
repo_root <- Sys.getenv("HF_SCC_ROOT", unset = "")
if (!nzchar(repo_root)) stop("Set HF_SCC_ROOT to the repository root.")
repo_root <- normalizePath(repo_root, mustWork = TRUE)
source(file.path(repo_root, "TES/01_source/01_plotting.r"))
source(file.path(repo_root, "scripts/external_data_mode.R"))
EXTERNAL_DATA_MODE <- publication_external_data_mode()
figure4_default_plot_path <- file.path(repo_root, "Figures/Figure 4/Rplots.pdf")
options(device = function(...) grDevices::pdf(
  file = publication_output_path(figure4_default_plot_path), ...))

# Loading files
tes_mutations_df <- readRDS(file.path(
  repo_root, "TES/02_data/02_processed/mutations_unique_dp.rds"))
tes_outliers <- readRDS(file.path(
  repo_root, "TES/02_data/02_processed/technical_outliers.rds"))
wes_mutations_df <- readRDS(file.path(
  repo_root, "WES/02_data/02_processed/mutations_unique_dp.rds"))
wes_outliers <- readRDS(file.path(
  repo_root, "WES/02_data/02_processed/technical_outliers.rds"))

#-------------------------------------------------------------------------------
### Parameters ###
#-------------------------------------------------------------------------------
FIGURE4_SOURCE_DIR <- file.path(repo_root, "Figures/Figure 4")
OUTPUT_DIR_LOGICAL <- Sys.getenv(
  "FIGURE4_OUTPUT_DIR",
  unset = FIGURE4_SOURCE_DIR)
OUTPUT_DIR <- dirname(publication_output_path(
  file.path(OUTPUT_DIR_LOGICAL, ".output-root")))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_GT_AF <- 0.01
LOLLIPOP_VAF_MODE <- "mean" # "mean" or "summed"
VAF_NULL_ITERATIONS <- 10000L
VAF_NULL_SEED <- 20260514L
VAF_NULL_REQUIRE_MIN_GT_AF <- FALSE
VAF_NULL_EXCLUDE_CANDIDATES_FROM_ITERATIONS <- FALSE
VAF_SUPPORT_MATCH_RECURRENCE <- FALSE
VAF_SUPPORT_SEED <- 20260515L
VAF_PLOT_LIMITS <- c(0, 0.1)
VAF_PLOT_BREAKS <- c(0, 0.05, 0.1)
MANUAL_CANDIDATE_GENE_EXCLUSIONS <- c(
  "1110059E24Rik",
  "Gm1989",
  "Vmn2r128")

CBIO_STUDY_IDS <- c(
  "cscc_dfarber_2015",
  "cscc_hgsc_bcm_2014",
  "cscc_ucsf_2021",
  "cscc_ranson_2022")
CBIO_API_BASE_URL <- "https://www.cbioportal.org/api"
ENSEMBL_REST_BASE_URL <- "https://rest.ensembl.org"
CBIO_RESTRICT_TO_NONSYNONYMOUS <- TRUE
CBIO_NONSYNONYMOUS_CLASSES <- c("Missense", "Inframe", "Truncating/splice")
CBIO_ASSUME_NO_PANEL_MEANS_WHOLE_PROFILE <- FALSE
CBIO_MUTATED_GENES_FILE <- file.path(
  FIGURE4_SOURCE_DIR,
  "cBioPortal_Mutated_Genes.txt")
CBIO_INCIDENCE_PERMUTATIONS <- 10000L
CBIO_INCIDENCE_PERMUTATION_SEED <- 20260516L
CBIO_EXTRA_ONCOPLOT_GROUP <- "cSCC drivers"
CBIO_EXTRA_ONCOPLOT_GENES <- c("Hras", "Notch1", "Tp53")
DRIVER_COMUTATION_GENES <- tibble(
  driver_label = c("Hras", "P53", "Notch1"),
  driver_gene_symbol = c("Hras", "Trp53", "Notch1"))

if (identical(EXTERNAL_DATA_MODE, "refresh")) {
  publication_require_live_refresh()
  if (!requireNamespace("cbioportalR", quietly = TRUE)) {
    stop("The explicit refresh path requires cbioportalR.")
  }
  if ("set_cbioportal_db" %in% getNamespaceExports("cbioportalR")) {
    cbioportalR::set_cbioportal_db("public")
  }
}

TUMOUR_CATEGORIES <- c("Papilloma", "SCC")
TIME_LEVELS <- names(col_palette$time)
CANDIDATE_GROUP_LEVELS <- c("Initiation (I)", "I+P", "Predisposition")
CBIO_CANDIDATE_GROUP_LEVELS <- c(
  "Initiation (I)",
  "I+P",
  "Predisposition (P)")
CBIO_ONCOPLOT_GROUP_LEVELS <- c(
  CBIO_CANDIDATE_GROUP_LEVELS,
  CBIO_EXTRA_ONCOPLOT_GROUP)
HF_UNKNOWN_PAPILLOMA_CONTEXT <- "HF Papilloma of unknown fate"
SKIN_UNKNOWN_PAPILLOMA_CONTEXT <- "Skin Papilloma of unknown fate"
HF_REGRESSING_PAPILLOMA_CONTEXT <- "HF Regressing papilloma"
SKIN_REGRESSING_PAPILLOMA_CONTEXT <- "Skin Regressing papilloma"
HF_SCC_CONTEXT <- "TES HF SCC"
SKIN_SCC_CONTEXT <- "TES skin SCC"

CONTEXT_LEVELS <- c(
  "Visually normal",
  "Acetone",
  HF_UNKNOWN_PAPILLOMA_CONTEXT,
  SKIN_UNKNOWN_PAPILLOMA_CONTEXT,
  HF_REGRESSING_PAPILLOMA_CONTEXT,
  SKIN_REGRESSING_PAPILLOMA_CONTEXT,
  HF_SCC_CONTEXT,
  SKIN_SCC_CONTEXT)

SKIN_VAF_CONTEXT_LEVELS <- c(
  "Visually normal",
  "Acetone",
  SKIN_UNKNOWN_PAPILLOMA_CONTEXT,
  SKIN_REGRESSING_PAPILLOMA_CONTEXT,
  SKIN_SCC_CONTEXT)

context_axis_labels <- function(x) {
  recode(as.character(x), "Visually normal" = "Morphologically normal")
}

candidate_matrix_context_levels <- c("MN", "HF Pap", "HF SCC", "SK Pap", "SK SCC")

candidate_matrix_context_label <- function(context) {
  recode(
    as.character(context),
    "Visually normal" = "MN",
    "HF Papilloma of unknown fate" = "HF Pap",
    "HF Regressing papilloma" = "HF Pap",
    "TES HF SCC" = "HF SCC",
    "Skin Papilloma of unknown fate" = "SK Pap",
    "Skin Regressing papilloma" = "SK Pap",
    "TES skin SCC" = "SK SCC",
    .default = as.character(context))
}

MANUAL_HUMAN_ORTHOLOGUES <- tibble(
  gene_symbol = c("C130026I21Rik", "Hras", "Notch1", "Tp53"),
  human_symbols = c("SP140", "HRAS", "NOTCH1", "TP53"))

if (!LOLLIPOP_VAF_MODE %in% c("mean", "summed")) {
  stop('LOLLIPOP_VAF_MODE must be either "mean" or "summed".')
}

#-------------------------------------------------------------------------------
### Small helpers ###
#-------------------------------------------------------------------------------
export_plot_default_panel_width_mm <- as.numeric(formals(export_plot_data)$panel_width)
export_plot_default_panel_height_mm <- as.numeric(formals(export_plot_data)$panel_height)
export_plot_default_tick_size_mm <- as.numeric(formals(export_plot_data)$tick_size)
export_plot_default_extra_size_mm <- 30

export_plot_default_dimension_mm <- function(n_panels, panel_size, tick_size) {
  n_panels * (panel_size + tick_size) + export_plot_default_extra_size_mm
}

set_plot_panel_size <- function(plot, panel_width = NULL, panel_height = NULL) {
  if (is.null(panel_width) && is.null(panel_height)) {
    return(plot)
  }

  panel_size_args <- list()

  if (!is.null(panel_width)) {
    panel_size_args$cols <- grid::unit(panel_width, "mm")
  }

  if (!is.null(panel_height)) {
    panel_size_args$rows <- grid::unit(panel_height, "mm")
  }

  plot + do.call(ggh4x::force_panelsizes, panel_size_args)
}

export_ggplot_data <- function(
    plot,
    data,
    file_name,
    width = NULL,
    height = NULL,
    panel_width = NULL,
    panel_height = NULL,
    tick_size = export_plot_default_tick_size_mm) {

  plot_layout <- ggplot_build(plot)$layout$layout
  n_panel_cols <- max(plot_layout$COL)
  n_panel_rows <- max(plot_layout$ROW)

  plot <- set_plot_panel_size(
    plot = plot,
    panel_width = panel_width,
    panel_height = panel_height)

  if (is.null(width) && !is.null(panel_width)) {
    width <- export_plot_default_dimension_mm(n_panel_cols, panel_width, tick_size)
  }

  if (is.null(height) && !is.null(panel_height)) {
    height <- export_plot_default_dimension_mm(n_panel_rows, panel_height, tick_size)
  }

  data_path <- publication_output_path(paste0(file_name, "_data.csv"))
  pdf_path <- publication_output_path(paste0(file_name, ".pdf"))
  readr::write_csv(data, data_path)
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "mm")

  invisible(data)
}

first_non_empty <- function(x) {
  x <- unique(na.omit(trimws(as.character(x))))
  x <- x[x != ""]

  if (length(x) == 0) NA_character_ else x[[1]]
}

collapse_unique <- function(x, n = 5, sep = " | ") {
  x <- unique(na.omit(trimws(as.character(x))))
  x <- x[x != ""]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(head(x, n), collapse = sep)
}

max_or_na <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) NA_real_ else max(x)
}

effect_label <- function(hgvsp, hgvsc, consequence) {
  hgvsp <- coalesce(as.character(hgvsp), "")
  hgvsc <- coalesce(as.character(hgvsc), "")
  consequence <- coalesce(as.character(consequence), "")
  short_consequence <- str_replace(consequence, "&.*", "")

  if_else(
    hgvsp != "",
    hgvsp,
    if_else(hgvsc != "", hgvsc, short_consequence))
}

earliest_time <- function(time_values) {
  time_values <- as.character(time_values)
  time_values <- time_values[time_values %in% TIME_LEVELS]

  if (length(time_values) == 0) {
    return(NA_character_)
  }

  TIME_LEVELS[min(match(time_values, TIME_LEVELS), na.rm = TRUE)]
}

is_priority_nonsynonymous <- function(consequence, impact) {
  priority_pattern <- paste(
    c(
      "missense_variant",
      "frameshift_variant",
      "protein_altering_variant",
      "inframe_insertion",
      "inframe_deletion",
      "splice_acceptor_variant",
      "splice_donor_variant",
      "splice_region_variant",
      "start_lost",
      "stop_gained",
      "stop_lost",
      "coding_sequence_variant"),
    collapse = "|")

  impact %in% c("HIGH", "MODERATE") |
    str_detect(coalesce(consequence, ""), priority_pattern)
}

active_vaf_col <- function() {
  if (LOLLIPOP_VAF_MODE == "summed") {
    return("mean_summed_vaf")
  }

  "mean_vaf"
}

active_vaf_label <- function() {
  if (LOLLIPOP_VAF_MODE == "summed") {
    return("Mean summed skin VAF")
  }

  "Median skin VAF"
}

vaf_label <- function(x) {
  formatC(x, format = "f", digits = 2)
}

candidate_group_levels <- function(extra_levels = character()) {
  unique(c(CANDIDATE_GROUP_LEVELS, as.character(extra_levels)))
}

candidate_group_display <- function(candidate_group) {
  recode(
    as.character(candidate_group),
    "Predisposition" = "Predisposition (P)",
    .default = as.character(candidate_group))
}

candidate_fill_values <- c(
  col_palette$candidate)

cbio_candidate_fill_values <- c(
  "Initiation (I)" = unname(col_palette$candidate["Initiation (I)"]),
  "I+P" = unname(col_palette$candidate["I+P"]),
  "Predisposition (P)" = unname(col_palette$candidate["Predisposition"]))

long_gene_panel_row_height_mm <- 3.6
long_gene_panel_extra_height_mm <- 44
long_gene_panel_min_height_mm <- 135

long_gene_panel_height_mm <- function(gene_symbols) {
  max(
    long_gene_panel_min_height_mm,
    long_gene_panel_row_height_mm * n_distinct(gene_symbols) +
      long_gene_panel_extra_height_mm)
}

sample_count_size_range <- c(1.9, 4.4)
sample_count_size_limits <- c(1, 4)
sample_count_size_breaks <- 1:4

vaf_null_group_levels <- c(
  "All candidates",
  "Initiation-only",
  "Predisposition-only")

drop_empty_contexts <- function(plot_df, count_col) {
  kept_contexts <- plot_df %>%
    group_by(context) %>%
    summarise(
      has_signal = any(.data[[count_col]] > 0, na.rm = TRUE),
      .groups = "drop") %>%
    filter(has_signal) %>%
    pull(context) %>%
    as.character()

  plot_df %>%
    filter(as.character(context) %in% kept_contexts) %>%
    mutate(context = factor(as.character(context), levels = kept_contexts))
}

validate_outlier_exclusion <- function(tbl, outlier_tbl, dataset_label, tbl_name) {
  leaked_samples <- tbl %>%
    filter(
      dataset == dataset_label,
      sample_name %in% outlier_tbl$sample_name) %>%
    distinct(sample_name) %>%
    pull(sample_name)

  if (length(leaked_samples) > 0) {
    stop(
      tbl_name,
      " contains technical outlier samples for ",
      dataset_label,
      ": ",
      paste(leaked_samples, collapse = ", "))
  }

  invisible(TRUE)
}

#-------------------------------------------------------------------------------
### Mutation preprocessing ###
#-------------------------------------------------------------------------------
add_identity_columns <- function(df, dataset_label) {
  df %>%
    mutate(
      dataset = dataset_label,
      sample_uid = paste(dataset_label, sample_name, sep = "::"),
      mouse_uid = paste(dataset_label, mouse, sep = "::"))
}

add_context_columns <- function(df) {
  df %>%
    mutate(
      context = case_when(
        treatment == "DT" & category == "Visually normal" ~
          "Visually normal",
        treatment == "Acetone" ~
          "Acetone",
        dataset == "WES" &
          treatment == "DT" &
          tissue == "Hair follicle" &
          category == "Papilloma" ~
          HF_UNKNOWN_PAPILLOMA_CONTEXT,
        dataset == "WES" &
          treatment == "DT" &
          tissue == "Skin" &
          category == "Papilloma" ~
          SKIN_UNKNOWN_PAPILLOMA_CONTEXT,
        dataset == "TES" &
          treatment == "DT" &
          tissue == "Hair follicle" &
          category == "Papilloma" ~
          HF_REGRESSING_PAPILLOMA_CONTEXT,
        dataset == "TES" &
          treatment == "DT" &
          tissue == "Skin" &
          category == "Papilloma" ~
          SKIN_REGRESSING_PAPILLOMA_CONTEXT,
        dataset == "TES" &
          treatment == "DT" &
          tissue == "Hair follicle" &
          category == "SCC" ~
          HF_SCC_CONTEXT,
        dataset == "TES" &
          treatment == "DT" &
          tissue == "Skin" &
          category == "SCC" ~
          SKIN_SCC_CONTEXT,
        TRUE ~ NA_character_))
}

prepare_mutations <- function(
    raw_tbl,
    outliers,
    dataset_label,
    use_last_cycle = FALSE,
    min_gt_af = MIN_GT_AF) {
  tbl <- raw_tbl %>%
    mutate(gt_AF = suppressWarnings(as.numeric(gt_AF)))

  if (!"time" %in% names(tbl)) {
    tbl$time <- NA_character_
  }

  if (use_last_cycle && "last_cycle" %in% names(tbl)) {
    tbl <- tbl %>% filter(last_cycle == TRUE)
  }

  tbl <- tbl %>%
    filter(
      !sample_name %in% outliers$sample_name,
      !is.na(sample_name),
      !is.na(gt_AF))

  if (!is.null(min_gt_af)) {
    tbl <- tbl %>% filter(gt_AF >= min_gt_af)
  }

  tbl %>%
    mutate(
      time = case_when(
        dataset_label == "WES" ~ "Week 12",
        as.character(time) == "Skin 19" ~ "Week 19",
        TRUE ~ as.character(time))) %>%
    mutate(
      mutation_id = paste(CHROM, POS, REF, ALT, sep = ":"),
      effect_short = effect_label(HGVSp, HGVSc, Consequence),
      mutation_brief = paste0(
        coalesce(SYMBOL, "NA"),
        " ",
        CHROM,
        ":",
        POS,
        " ",
        REF,
        ">",
        ALT,
        " ",
        effect_short)) %>%
    add_identity_columns(dataset_label) %>%
    add_context_columns() %>%
    filter(
      !(
        context %in% c("Visually normal", "Acetone") &
          IMPACT %in% c("LOW", "MODIFIER")
      ))
}

build_nonsyn_base <- function(presence_tbl) {
  presence_tbl %>%
    filter(
      !is.na(SYMBOL),
      SYMBOL != "",
      Consequence != "intergenic_variant",
      Feature_type != "RegulatoryFeature",
      !str_detect(coalesce(Consequence, ""), "synonymous_variant"),
      is_priority_nonsynonymous(Consequence, IMPACT)) %>%
    distinct(sample_uid, mutation_id, SYMBOL, .keep_all = TRUE)
}

tes_presence <- prepare_mutations(
  raw_tbl = tes_mutations_df,
  outliers = tes_outliers,
  dataset_label = "TES")

wes_presence <- prepare_mutations(
  raw_tbl = wes_mutations_df,
  outliers = wes_outliers,
  dataset_label = "WES",
  use_last_cycle = TRUE)

presence_base <- bind_rows(tes_presence, wes_presence)
nonsyn_base <- build_nonsyn_base(presence_base)

validate_outlier_exclusion(presence_base, tes_outliers, "TES", "presence_base")
validate_outlier_exclusion(presence_base, wes_outliers, "WES", "presence_base")
validate_outlier_exclusion(nonsyn_base, tes_outliers, "TES", "nonsyn_base")
validate_outlier_exclusion(nonsyn_base, wes_outliers, "WES", "nonsyn_base")

#-------------------------------------------------------------------------------
### Candidate genes ###
#-------------------------------------------------------------------------------
assign_candidate_group <- function(initiation, progression) {
  case_when(
    initiation & progression ~ "I+P",
    initiation ~ "Initiation (I)",
    progression ~ "Predisposition",
    TRUE ~ NA_character_)
}

split_human_symbols <- function(orthologue_tbl) {
  orthologue_tbl %>%
    select(gene_symbol, human_symbols) %>%
    mutate(human_symbols = coalesce(human_symbols, "")) %>%
    separate_rows(human_symbols, sep = ";") %>%
    transmute(
      gene_symbol,
      human_symbol = trimws(human_symbols)) %>%
    filter(human_symbol != "")
}

query_babelgene <- function(gene_symbols, min_support) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    return(tibble(gene_symbol = character(), human_symbol = character()))
  }

  gene_symbols <- sort(unique(na.omit(gene_symbols)))
  if (length(gene_symbols) == 0) {
    return(tibble(gene_symbol = character(), human_symbol = character()))
  }

  ortholog_args <- list(genes = gene_symbols, species = "mouse", human = FALSE)
  ortholog_formals <- names(formals(babelgene::orthologs))

  if ("min_support" %in% ortholog_formals) {
    ortholog_args$min_support <- min_support
  }

  if ("top" %in% ortholog_formals) {
    ortholog_args$top <- TRUE
  }

  tryCatch(
    do.call(babelgene::orthologs, ortholog_args) %>%
      as_tibble() %>%
      filter(!is.na(symbol), !is.na(human_symbol), human_symbol != "") %>%
      transmute(gene_symbol = symbol, human_symbol) %>%
      distinct(),
    error = function(e) {
      warning("babelgene orthologue lookup failed: ", conditionMessage(e))
      tibble(gene_symbol = character(), human_symbol = character())
    })
}

build_human_orthologues <- function(gene_symbols) {
  if (identical(EXTERNAL_DATA_MODE, "frozen")) {
    frozen <- readr::read_csv(
      publication_repo_input(
        "Figures/Figure 4/Figure4_cbioportal_candidate_oncoplot_data.csv",
        "eb63bd4ee5a87c8b0ca06e782a8e2482634f72671764851644e20a65066c6ea4"),
      show_col_types = FALSE) %>%
      select(gene_symbol, human_symbols) %>%
      distinct()
    return(frozen %>%
      right_join(
        tibble(gene_symbol = sort(unique(na.omit(gene_symbols)))),
        by = "gene_symbol") %>%
      mutate(
        human_symbols = coalesce(human_symbols, ""),
        has_human_orthologue = human_symbols != ""))
  }

  publication_require_live_refresh()
  manual_hits <- MANUAL_HUMAN_ORTHOLOGUES %>%
    filter(gene_symbol %in% gene_symbols) %>%
    split_human_symbols()

  high_support <- query_babelgene(
    setdiff(gene_symbols, manual_hits$gene_symbol),
    min_support = 3)

  low_support <- query_babelgene(
    setdiff(gene_symbols, c(manual_hits$gene_symbol, high_support$gene_symbol)),
    min_support = 1)

  bind_rows(manual_hits, high_support, low_support) %>%
    distinct() %>%
    right_join(
      tibble(gene_symbol = sort(unique(na.omit(gene_symbols)))),
      by = "gene_symbol") %>%
    group_by(gene_symbol) %>%
    summarise(
      human_symbols = paste(sort(unique(na.omit(human_symbol))), collapse = ";"),
      .groups = "drop") %>%
    mutate(
      human_symbols = coalesce(human_symbols, ""),
      has_human_orthologue = human_symbols != "")
}

format_gene_display_label <- function(gene_symbol, human_symbols) {
  human_display <- str_replace_all(coalesce(human_symbols, ""), ";", "/")

  if_else(
    human_display != "",
    paste(gene_symbol, human_display, sep = " / "),
    gene_symbol)
}

add_gene_display_labels <- function(tbl) {
  tbl %>%
    mutate(
      row_label = format_gene_display_label(gene_symbol, human_symbols),
      gene_label = row_label)
}

escape_plotmath_string <- function(x) {
  str_replace_all(as.character(x), c("\\\\" = "\\\\\\\\", "\"" = "\\\\\""))
}

italic_gene_axis_labels <- function(labels) {
  label_text <- vapply(as.character(labels), function(label) {
    parts <- str_split(label, " / ", n = 2, simplify = TRUE)

    if (ncol(parts) >= 2 && parts[[1, 2]] != "") {
      return(sprintf(
        'italic("%s")~"/"~italic("%s")',
        escape_plotmath_string(parts[[1, 1]]),
        escape_plotmath_string(parts[[1, 2]])))
    }

    sprintf('italic("%s")', escape_plotmath_string(label))
  }, character(1))

  parse(text = label_text)
}

annotate_human_orthologues <- function(tbl, orthologue_tbl) {
  tbl %>%
    left_join(orthologue_tbl, by = "gene_symbol") %>%
    mutate(
      human_symbols = coalesce(human_symbols, ""),
      has_human_orthologue = replace_na(has_human_orthologue, FALSE)) %>%
    add_gene_display_labels()
}

summarise_candidate_genes <- function(evidence_tbl, presence_tbl) {
  gene_absence <- presence_tbl %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    group_by(gene_symbol = SYMBOL) %>%
    summarise(
      visually_normal_modhigh_samples = n_distinct(
        sample_uid[
          context == "Visually normal" &
            IMPACT %in% c("HIGH", "MODERATE")]),
      tes_papilloma_modhigh_samples = n_distinct(
        sample_uid[
          context %in% c(
            HF_REGRESSING_PAPILLOMA_CONTEXT,
            SKIN_REGRESSING_PAPILLOMA_CONTEXT) &
            IMPACT %in% c("HIGH", "MODERATE")]),
      .groups = "drop")

  evidence_tbl %>%
    group_by(gene_symbol = SYMBOL) %>%
    summarise(
      consequence = collapse_unique(str_replace(Consequence, "&.*", "")),
      protein_change = collapse_unique(effect_short),
      representative_mutations = collapse_unique(
        mutation_brief[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      hf_tumour_samples = n_distinct(
        sample_uid[
          treatment == "DT" &
            tissue == "Hair follicle" &
            category %in% TUMOUR_CATEGORIES]),
      skin_tumour_samples = n_distinct(
        sample_uid[
          treatment == "DT" &
            tissue == "Skin" &
            category %in% TUMOUR_CATEGORIES]),
      tumour_samples = n_distinct(
        sample_uid[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      n_mutations = n_distinct(
        mutation_id[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      max_gt_AF = max_or_na(
        gt_AF[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      first_time = earliest_time(
        time[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      tes_hf_scc_samples = n_distinct(sample_uid[context == HF_SCC_CONTEXT]),
      tes_skin_scc_samples = n_distinct(sample_uid[context == SKIN_SCC_CONTEXT]),
      wes_hf_papilloma_samples = n_distinct(
        sample_uid[context == HF_UNKNOWN_PAPILLOMA_CONTEXT]),
      wes_skin_papilloma_samples = n_distinct(
        sample_uid[context == SKIN_UNKNOWN_PAPILLOMA_CONTEXT]),
      .groups = "drop") %>%
    left_join(gene_absence, by = "gene_symbol") %>%
    mutate(
      across(
        c(visually_normal_modhigh_samples, tes_papilloma_modhigh_samples),
        ~ replace_na(.x, 0L)),
      initiation = hf_tumour_samples >= 1 &
        skin_tumour_samples >= 1 &
        visually_normal_modhigh_samples == 0,
      progression = hf_tumour_samples >= 1 &
        skin_tumour_samples >= 1 &
        tes_papilloma_modhigh_samples == 0 &
        (
          (tes_hf_scc_samples >= 1 & tes_skin_scc_samples >= 1) |
            (tes_hf_scc_samples >= 1 & wes_skin_papilloma_samples >= 1) |
            (wes_hf_papilloma_samples >= 1 & tes_skin_scc_samples >= 1)
        ),
      candidate_group = assign_candidate_group(initiation, progression)) %>%
    filter(!is.na(candidate_group)) %>%
    mutate(
      row_label = gene_symbol,
      gene_label = row_label,
      candidate_group = factor(candidate_group, levels = CANDIDATE_GROUP_LEVELS))
}

build_gene_table_for_symbols <- function(
    gene_symbols,
    group_label,
    evidence_tbl,
    presence_tbl,
    orthologue_tbl) {

  gene_symbols <- unique(na.omit(gene_symbols))

  evidence_summary <- evidence_tbl %>%
    filter(SYMBOL %in% gene_symbols) %>%
    group_by(gene_symbol = SYMBOL) %>%
    summarise(
      consequence = collapse_unique(str_replace(Consequence, "&.*", "")),
      protein_change = collapse_unique(effect_short),
      representative_mutations = collapse_unique(
        mutation_brief[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      hf_tumour_samples = n_distinct(
        sample_uid[
          treatment == "DT" &
            tissue == "Hair follicle" &
            category %in% TUMOUR_CATEGORIES]),
      skin_tumour_samples = n_distinct(
        sample_uid[
          treatment == "DT" &
            tissue == "Skin" &
            category %in% TUMOUR_CATEGORIES]),
      tumour_samples = n_distinct(
        sample_uid[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      n_mutations = n_distinct(
        mutation_id[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      max_gt_AF = max_or_na(
        gt_AF[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      first_time = earliest_time(
        time[treatment == "DT" & category %in% TUMOUR_CATEGORIES]),
      .groups = "drop")

  absence_summary <- presence_tbl %>%
    filter(SYMBOL %in% gene_symbols) %>%
    group_by(gene_symbol = SYMBOL) %>%
    summarise(
      visually_normal_modhigh_samples = n_distinct(
        sample_uid[
          context == "Visually normal" &
            IMPACT %in% c("HIGH", "MODERATE")]),
      tes_papilloma_modhigh_samples = n_distinct(
        sample_uid[
          context %in% c(
            HF_REGRESSING_PAPILLOMA_CONTEXT,
            SKIN_REGRESSING_PAPILLOMA_CONTEXT) &
            IMPACT %in% c("HIGH", "MODERATE")]),
      .groups = "drop")

  tibble(gene_symbol = gene_symbols) %>%
    left_join(evidence_summary, by = "gene_symbol") %>%
    left_join(absence_summary, by = "gene_symbol") %>%
    mutate(
      candidate_group = group_label,
      row_label = gene_symbol,
      gene_label = row_label,
      initiation = FALSE,
      progression = FALSE,
      across(
        c(
          hf_tumour_samples,
          skin_tumour_samples,
          tumour_samples,
          n_mutations,
          visually_normal_modhigh_samples,
          tes_papilloma_modhigh_samples),
        ~ replace_na(.x, 0L)),
      consequence = coalesce(consequence, ""),
      protein_change = coalesce(protein_change, ""),
      representative_mutations = coalesce(representative_mutations, "")) %>%
    annotate_human_orthologues(orthologue_tbl)
}

make_gene_matrix <- function(gene_tbl, context_tbl) {
  observed <- context_tbl %>%
    filter(
      SYMBOL %in% gene_tbl$gene_symbol,
      context %in% CONTEXT_LEVELS) %>%
    group_by(gene_symbol = SYMBOL, context) %>%
    summarise(
      first_time = earliest_time(time),
      datasets = collapse_unique(dataset, n = Inf, sep = " + "),
      n_samples = n_distinct(sample_uid),
      n_mutations = n_distinct(mutation_id),
      max_gt_AF = max_or_na(gt_AF),
      .groups = "drop")

  expand_grid(
    gene_symbol = gene_tbl$gene_symbol,
    context = CONTEXT_LEVELS) %>%
    left_join(
      gene_tbl %>%
        select(
          gene_symbol,
          row_label,
          gene_label,
          candidate_group,
          human_symbols,
          has_human_orthologue,
          consequence,
          protein_change,
          representative_mutations,
          initiation,
          progression,
          hf_tumour_samples,
          skin_tumour_samples,
          visually_normal_modhigh_samples,
          tes_papilloma_modhigh_samples),
      by = "gene_symbol") %>%
    left_join(observed, by = c("gene_symbol", "context")) %>%
    mutate(
      n_samples = replace_na(n_samples, 0L),
      n_mutations = replace_na(n_mutations, 0L),
      detected = n_samples > 0,
      context = factor(context, levels = CONTEXT_LEVELS),
      first_time = factor(first_time, levels = TIME_LEVELS))
}

format_candidate_matrix_contexts <- function(matrix_tbl) {
  matrix_tbl %>%
    filter(as.character(context) != "Acetone") %>%
    mutate(context = candidate_matrix_context_label(context)) %>%
    group_by(
      gene_symbol,
      context,
      row_label,
      gene_label,
      candidate_group,
      human_symbols,
      has_human_orthologue,
      consequence,
      protein_change,
      representative_mutations,
      initiation,
      progression,
      hf_tumour_samples,
      skin_tumour_samples,
      visually_normal_modhigh_samples,
      tes_papilloma_modhigh_samples) %>%
    summarise(
      first_time = earliest_time(first_time),
      datasets = collapse_unique(datasets, n = Inf, sep = " + "),
      n_samples = sum(n_samples, na.rm = TRUE),
      n_mutations = sum(n_mutations, na.rm = TRUE),
      max_gt_AF = max_or_na(max_gt_AF),
      detected = any(detected, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(
      context = factor(context, levels = candidate_matrix_context_levels),
      first_time = factor(first_time, levels = TIME_LEVELS)) %>%
    arrange(gene_symbol, context)
}

mouse_human_orthologues <- build_human_orthologues(
  c(nonsyn_base$SYMBOL, CBIO_EXTRA_ONCOPLOT_GENES))

candidate_genes <- summarise_candidate_genes(
  evidence_tbl = nonsyn_base,
  presence_tbl = presence_base) %>%
  annotate_human_orthologues(mouse_human_orthologues)

manually_excluded_candidate_genes <- candidate_genes %>%
  filter(gene_symbol %in% MANUAL_CANDIDATE_GENE_EXCLUSIONS) %>%
  arrange(gene_symbol)

candidate_genes <- candidate_genes %>%
  filter(!gene_symbol %in% MANUAL_CANDIDATE_GENE_EXCLUSIONS)

if (any(
  candidate_genes$initiation &
    candidate_genes$visually_normal_modhigh_samples > 0)) {
  stop("Initiation candidates include visually-normal MODERATE/HIGH hits.")
}

if (any(
  candidate_genes$progression &
    candidate_genes$tes_papilloma_modhigh_samples > 0)) {
  stop("Predisposition candidates include TES papilloma MODERATE/HIGH hits.")
}

#-------------------------------------------------------------------------------
### Driver co-mutated genes in skin tumour regions ###
#-------------------------------------------------------------------------------
build_driver_comutation_table <- function(evidence_tbl, driver_tbl) {
  skin_tumour_hits <- evidence_tbl %>%
    filter(
      tissue == "Skin",
      treatment == "DT",
      category %in% TUMOUR_CATEGORIES,
      !is.na(SYMBOL),
      SYMBOL != "") %>%
    distinct(
      sample_uid,
      sample_name,
      mouse_uid,
      mouse,
      dataset,
      context,
      category,
      time,
      SYMBOL,
      mutation_id,
      gt_AF,
      mutation_brief,
      Consequence,
      IMPACT,
      effect_short,
      .keep_all = TRUE)

  visually_normal_skin_genes <- evidence_tbl %>%
    filter(
      tissue == "Skin",
      context == "Visually normal",
      !is.na(SYMBOL),
      SYMBOL != "") %>%
    distinct(SYMBOL) %>%
    pull(SYMBOL)

  driver_hits <- skin_tumour_hits %>%
    inner_join(driver_tbl, by = c("SYMBOL" = "driver_gene_symbol")) %>%
    transmute(
      sample_uid,
      driver_label,
      driver_gene_symbol = SYMBOL,
      driver_mutation_id = mutation_id,
      driver_gt_AF = gt_AF,
      driver_mutation = mutation_brief)

  partner_hits <- skin_tumour_hits %>%
    filter(!SYMBOL %in% visually_normal_skin_genes) %>%
    transmute(
      sample_uid,
      sample_name,
      mouse_uid,
      mouse,
      dataset,
      context,
      category,
      time,
      gene_symbol = SYMBOL,
      mutation_id,
      gt_AF,
      consequence = str_replace(Consequence, "&.*", ""),
      protein_change = effect_short,
      mutation = mutation_brief)

  driver_hits %>%
    inner_join(partner_hits, by = "sample_uid") %>%
    filter(gene_symbol != driver_gene_symbol) %>%
    group_by(driver_label, gene_symbol) %>%
    summarise(
      n_comutated_samples = n_distinct(sample_uid),
      samples = collapse_unique(sample_name, n = Inf, sep = ";"),
      mice = collapse_unique(mouse, n = Inf, sep = ";"),
      datasets = collapse_unique(dataset, n = Inf, sep = ";"),
      contexts = collapse_unique(context, n = Inf, sep = ";"),
      categories = collapse_unique(category, n = Inf, sep = ";"),
      first_time = earliest_time(time),
      partner_consequence = collapse_unique(consequence, n = Inf, sep = ";"),
      partner_protein_change = collapse_unique(protein_change, n = Inf, sep = ";"),
      n_partner_mutations = n_distinct(mutation_id),
      max_partner_gt_AF = max_or_na(gt_AF),
      partner_mutations = collapse_unique(mutation, n = Inf),
      driver_gene_symbols = collapse_unique(driver_gene_symbol, n = Inf, sep = ";"),
      n_driver_mutations = n_distinct(driver_mutation_id),
      max_driver_gt_AF = max_or_na(driver_gt_AF),
      driver_mutations = collapse_unique(driver_mutation, n = Inf),
      .groups = "drop") %>%
    filter(n_comutated_samples >= 2) %>%
    arrange(driver_label, desc(n_comutated_samples), gene_symbol)
}

driver_comutated_skin_tumour_genes <- build_driver_comutation_table(
  evidence_tbl = nonsyn_base,
  driver_tbl = DRIVER_COMUTATION_GENES) %>%
  left_join(mouse_human_orthologues, by = "gene_symbol") %>%
  mutate(
    human_symbols = coalesce(human_symbols, ""),
    has_human_orthologue = replace_na(has_human_orthologue, FALSE)) %>%
  relocate(human_symbols, has_human_orthologue, .after = gene_symbol)

write_csv(
  driver_comutated_skin_tumour_genes,
  file.path(OUTPUT_DIR, "Figure4_driver_comutated_skin_tumour_genes_data.csv"))

format_p_value <- function(x) {
  if_else(x < 0.001, "P<0.001", paste0("P=", formatC(x, format = "f", digits = 3)))
}

#-------------------------------------------------------------------------------
### Skin VAF summaries ###
#-------------------------------------------------------------------------------
skin_sample_meta <- presence_base %>%
  filter(
    tissue == "Skin",
    context %in% SKIN_VAF_CONTEXT_LEVELS) %>%
  distinct(
    sample_uid,
    sample_name,
    mouse_uid,
    mouse,
    dataset,
    context,
    category) %>%
  mutate(
    context = factor(context, levels = SKIN_VAF_CONTEXT_LEVELS),
    category = factor(
      category,
      levels = c("Visually normal", "Acetone", "Papilloma", "SCC"))) %>%
  arrange(context, sample_uid)

build_skin_gene_sample_vaf <- function(gene_tbl, evidence_tbl, sample_meta) {
  gene_tbl <- gene_tbl %>%
    distinct(
      gene_symbol,
      gene_label,
      candidate_group,
      human_symbols,
      has_human_orthologue)

  evidence_tbl %>%
    filter(
      tissue == "Skin",
      context %in% SKIN_VAF_CONTEXT_LEVELS,
      SYMBOL %in% gene_tbl$gene_symbol) %>%
    distinct(sample_uid, mutation_id, SYMBOL, .keep_all = TRUE) %>%
    group_by(
      sample_uid,
      sample_name,
      mouse_uid,
      mouse,
      dataset,
      context,
      category,
      gene_symbol = SYMBOL) %>%
    summarise(
      summed_vaf = sum(gt_AF, na.rm = TRUE),
      mean_vaf = mean(gt_AF, na.rm = TRUE),
      n_mutations = n_distinct(mutation_id),
      max_gt_AF = max_or_na(gt_AF),
      .groups = "drop") %>%
    right_join(
      expand_grid(
        sample_uid = sample_meta$sample_uid,
        gene_symbol = gene_tbl$gene_symbol),
      by = c("sample_uid", "gene_symbol")) %>%
    left_join(
      sample_meta,
      by = "sample_uid",
      suffix = c("", ".sample")) %>%
    left_join(gene_tbl, by = "gene_symbol") %>%
    mutate(
      summed_vaf = replace_na(summed_vaf, 0),
      mean_vaf = replace_na(mean_vaf, 0),
      n_mutations = replace_na(n_mutations, 0L),
      sample_name = coalesce(sample_name, sample_name.sample),
      mouse_uid = coalesce(mouse_uid, mouse_uid.sample),
      mouse = coalesce(mouse, mouse.sample),
      dataset = coalesce(dataset, dataset.sample),
      context = coalesce(context, context.sample),
      category = coalesce(category, category.sample),
      context = factor(context, levels = SKIN_VAF_CONTEXT_LEVELS),
      category = factor(
        category,
        levels = c("Visually normal", "Acetone", "Papilloma", "SCC")))
}

summarise_skin_gene_vaf <- function(sample_vaf_tbl) {
  sample_vaf_tbl %>%
    group_by(
      gene_symbol,
      gene_label,
      candidate_group,
      human_symbols,
      has_human_orthologue) %>%
    summarise(
      n_samples_mutated = sum(summed_vaf > 0, na.rm = TRUE),
      mean_summed_vaf = if_else(
        n_samples_mutated > 0,
        mean(summed_vaf[summed_vaf > 0], na.rm = TRUE),
        0),
      mean_vaf = if_else(
        n_samples_mutated > 0,
        median(mean_vaf[summed_vaf > 0], na.rm = TRUE),
        0),
      variance_summed_vaf = if_else(
        n_samples_mutated > 1,
        var(summed_vaf[summed_vaf > 0], na.rm = TRUE),
        NA_real_),
      max_summed_vaf = max_or_na(summed_vaf),
      n_mutations = sum(n_mutations, na.rm = TRUE),
      .groups = "drop")
}

set_lollipop_order <- function(vaf_tbl, extra_group_levels = character()) {
  vaf_col <- active_vaf_col()

  vaf_tbl %>%
    mutate(
      candidate_group = factor(
        as.character(candidate_group),
        levels = candidate_group_levels(extra_group_levels))) %>%
    arrange(candidate_group, .data[[vaf_col]], gene_symbol) %>%
    mutate(gene_label = factor(gene_label, levels = rev(unique(gene_label))))
}

skin_candidate_gene_vaf <- build_skin_gene_sample_vaf(
  gene_tbl = candidate_genes,
  evidence_tbl = nonsyn_base,
  sample_meta = skin_sample_meta) %>%
  summarise_skin_gene_vaf()

#-------------------------------------------------------------------------------
### Skin VAF recurrence-matched null ###
#-------------------------------------------------------------------------------
build_vaf_null_nonsyn_base <- function() {
  if (VAF_NULL_REQUIRE_MIN_GT_AF) {
    return(nonsyn_base)
  }

  bind_rows(
    prepare_mutations(
      raw_tbl = tes_mutations_df,
      outliers = tes_outliers,
      dataset_label = "TES",
      min_gt_af = NULL),
    prepare_mutations(
      raw_tbl = wes_mutations_df,
      outliers = wes_outliers,
      dataset_label = "WES",
      use_last_cycle = TRUE,
      min_gt_af = NULL)) %>%
    build_nonsyn_base()
}

summarise_skin_gene_mean_vaf <- function(evidence_tbl, sample_meta) {
  evidence_tbl %>%
    filter(
      tissue == "Skin",
      sample_uid %in% sample_meta$sample_uid,
      context %in% SKIN_VAF_CONTEXT_LEVELS,
      !is.na(SYMBOL),
      SYMBOL != "") %>%
    distinct(sample_uid, mutation_id, SYMBOL, .keep_all = TRUE) %>%
    group_by(gene_symbol = SYMBOL, sample_uid) %>%
    summarise(
      sample_gene_mean_vaf = mean(gt_AF, na.rm = TRUE),
      .groups = "drop") %>%
    group_by(gene_symbol) %>%
    summarise(
      recurrence_n_skin_samples = n_distinct(sample_uid),
      gene_level_mean_vaf = mean(sample_gene_mean_vaf, na.rm = TRUE),
      .groups = "drop")
}

build_vaf_null_candidate_sets <- function(candidate_tbl) {
  bind_rows(
    candidate_tbl %>%
      transmute(
        gene_symbol,
        candidate_group = "All candidates"),
    candidate_tbl %>%
      filter(candidate_group == "Initiation (I)") %>%
      transmute(
        gene_symbol,
        candidate_group = "Initiation-only"),
    candidate_tbl %>%
      filter(candidate_group == "Predisposition") %>%
      transmute(
        gene_symbol,
        candidate_group = "Predisposition-only")) %>%
    mutate(candidate_group = factor(candidate_group, levels = vaf_null_group_levels))
}

collapse_recurrence_profile <- function(recurrence_values) {
  recurrence_values <- recurrence_values[!is.na(recurrence_values)]

  if (length(recurrence_values) == 0) {
    return("")
  }

  profile_tbl <- tibble(recurrence_n_skin_samples = recurrence_values) %>%
    count(recurrence_n_skin_samples, name = "n_genes") %>%
    arrange(recurrence_n_skin_samples)

  paste0(
    profile_tbl$n_genes,
    " genes with ",
    profile_tbl$recurrence_n_skin_samples,
    " skin samples",
    collapse = "; ")
}

summarise_vaf_null_candidates <- function(candidate_sets, gene_vaf_tbl) {
  candidate_sets %>%
    left_join(gene_vaf_tbl, by = "gene_symbol") %>%
    group_by(candidate_group) %>%
    summarise(
      n_genes = n(),
      observed_median_gene_mean_VAF = median(gene_level_mean_vaf, na.rm = TRUE),
      observed_mean_gene_mean_VAF = mean(gene_level_mean_vaf, na.rm = TRUE),
      observed_IQR_gene_mean_VAF = IQR(gene_level_mean_vaf, na.rm = TRUE),
      recurrence_profile = collapse_recurrence_profile(recurrence_n_skin_samples),
      .groups = "drop")
}

summarise_vaf_recurrence_pools <- function(
    gene_vaf_tbl,
    candidate_tbl,
    candidate_recurrence_tbl = gene_vaf_tbl) {
  all_candidates <- unique(candidate_tbl$gene_symbol)

  pool_summary <- gene_vaf_tbl %>%
    group_by(recurrence_n_skin_samples) %>%
    summarise(
      n_eligible_genes_total = n_distinct(gene_symbol),
      n_eligible_candidate_genes = n_distinct(
        gene_symbol[gene_symbol %in% all_candidates]),
      n_noncandidate_genes = n_distinct(
        gene_symbol[!gene_symbol %in% all_candidates]),
      genes_in_pool = paste(sort(unique(gene_symbol)), collapse = ";"),
      candidate_genes_in_null_bin = paste(
        sort(unique(gene_symbol[gene_symbol %in% all_candidates])),
        collapse = ";"),
      .groups = "drop") %>%
    arrange(recurrence_n_skin_samples)

  candidate_bin_summary <- candidate_tbl %>%
    distinct(gene_symbol) %>%
    left_join(
      candidate_recurrence_tbl %>%
        select(gene_symbol, recurrence_n_skin_samples),
      by = "gene_symbol") %>%
    group_by(recurrence_n_skin_samples) %>%
    summarise(
      n_candidate_genes = n_distinct(gene_symbol),
      candidate_genes_in_bin = paste(sort(unique(gene_symbol)), collapse = ";"),
      .groups = "drop")

  pool_summary %>%
    full_join(candidate_bin_summary, by = "recurrence_n_skin_samples") %>%
    mutate(
      n_eligible_genes_total = replace_na(n_eligible_genes_total, 0L),
      n_eligible_candidate_genes = replace_na(n_eligible_candidate_genes, 0L),
      n_noncandidate_genes = replace_na(n_noncandidate_genes, 0L),
      n_candidate_genes = replace_na(n_candidate_genes, 0L),
      genes_in_pool = replace_na(genes_in_pool, ""),
      candidate_genes_in_null_bin = replace_na(candidate_genes_in_null_bin, ""),
      candidate_genes_in_bin = replace_na(candidate_genes_in_bin, "")) %>%
    arrange(recurrence_n_skin_samples)
}

check_vaf_null_pool_sizes <- function(
    candidate_sets,
    gene_vaf_tbl,
    exclude_candidates_from_null,
    candidate_tbl,
    candidate_recurrence_tbl = gene_vaf_tbl) {

  all_candidates <- unique(candidate_tbl$gene_symbol)

  pool_tbl <- gene_vaf_tbl
  if (exclude_candidates_from_null) {
    pool_tbl <- pool_tbl %>% filter(!gene_symbol %in% all_candidates)
  }

  required_tbl <- candidate_sets %>%
    left_join(
      candidate_recurrence_tbl %>%
        select(gene_symbol, recurrence_n_skin_samples),
      by = "gene_symbol") %>%
    count(candidate_group, recurrence_n_skin_samples, name = "n_required")

  available_tbl <- pool_tbl %>%
    count(recurrence_n_skin_samples, name = "n_available")

  size_check <- required_tbl %>%
    left_join(available_tbl, by = "recurrence_n_skin_samples") %>%
    mutate(n_available = replace_na(n_available, 0L))

  undersized <- size_check %>%
    filter(n_available < n_required)

  if (nrow(undersized) > 0) {
    warning(
      "Some recurrence bins have fewer null-pool genes than requested; ",
      "sampling will still use replacement. Bins: ",
      paste(
        paste0(
          undersized$candidate_group,
          " recurrence ",
          undersized$recurrence_n_skin_samples,
          " needs ",
          undersized$n_required,
          ", has ",
          undersized$n_available),
        collapse = "; "),
      call. = FALSE)
  }

  empty_bins <- size_check %>% filter(n_available == 0)

  if (nrow(empty_bins) > 0) {
    stop(
      "Cannot run VAF null: at least one recurrence bin has zero eligible genes.")
  }

  invisible(size_check)
}

run_vaf_recurrence_matched_null <- function(
    candidate_sets,
    gene_vaf_tbl,
    n_iterations,
    seed,
    exclude_candidates_from_null,
    candidate_tbl,
    candidate_recurrence_tbl = gene_vaf_tbl) {

  all_candidates <- unique(candidate_tbl$gene_symbol)
  pool_tbl <- gene_vaf_tbl

  if (exclude_candidates_from_null) {
    pool_tbl <- pool_tbl %>% filter(!gene_symbol %in% all_candidates)
  }

  check_vaf_null_pool_sizes(
    candidate_sets = candidate_sets,
    gene_vaf_tbl = gene_vaf_tbl,
    exclude_candidates_from_null = exclude_candidates_from_null,
    candidate_tbl = candidate_tbl,
    candidate_recurrence_tbl = candidate_recurrence_tbl)

  required_tbl <- candidate_sets %>%
    left_join(
      candidate_recurrence_tbl %>%
        select(gene_symbol, recurrence_n_skin_samples),
      by = "gene_symbol") %>%
    count(candidate_group, recurrence_n_skin_samples, name = "n_to_draw")

  active_group_levels <- candidate_sets %>%
    count(candidate_group) %>%
    filter(n > 0) %>%
    pull(candidate_group) %>%
    as.character()

  draw_plan <- lapply(active_group_levels, function(group_label) {
    group_required <- required_tbl %>%
      filter(candidate_group == group_label)

    split(group_required, seq_len(nrow(group_required))) %>%
      lapply(function(bin_tbl) {
        recurrence_bin <- bin_tbl$recurrence_n_skin_samples[[1]]
        list(
          n_to_draw = bin_tbl$n_to_draw[[1]],
          vaf_values = pool_tbl %>%
            filter(recurrence_n_skin_samples == recurrence_bin) %>%
            pull(gene_level_mean_vaf))
      })
  })
  names(draw_plan) <- active_group_levels

  set.seed(seed + as.integer(exclude_candidates_from_null))

  bind_rows(lapply(active_group_levels, function(group_label) {
    tibble(
      iteration = seq_len(n_iterations),
      candidate_group = group_label,
      exclude_candidates_from_null = exclude_candidates_from_null,
      null_statistic = replicate(n_iterations, {
        sampled_vaf <- unlist(lapply(draw_plan[[group_label]], function(bin_plan) {
          sample(
            bin_plan$vaf_values,
            bin_plan$n_to_draw,
            replace = TRUE)
        }), use.names = FALSE)

        median(sampled_vaf, na.rm = TRUE)
      }))
  })) %>%
    mutate(candidate_group = factor(candidate_group, levels = vaf_null_group_levels))
}

add_vaf_null_observed_difference <- function(null_tbl, observed_tbl) {
  null_tbl %>%
    select(-any_of(c("observed_statistic", "observed_minus_null_statistic"))) %>%
    left_join(
      observed_tbl %>%
        select(
          candidate_group,
          observed_statistic = observed_median_gene_mean_VAF),
      by = "candidate_group") %>%
    mutate(
      observed_minus_null_statistic = observed_statistic - null_statistic)
}

summarise_vaf_null_test <- function(null_tbl, observed_tbl) {
  add_vaf_null_observed_difference(
    null_tbl = null_tbl,
    observed_tbl = observed_tbl) %>%
    group_by(candidate_group, exclude_candidates_from_null) %>%
    summarise(
      n_iterations = n_distinct(iteration),
      observed_statistic = dplyr::first(observed_statistic),
      null_mean = mean(null_statistic, na.rm = TRUE),
      null_median = median(null_statistic, na.rm = TRUE),
      `null_q2.5` = unname(quantile(null_statistic, 0.025, na.rm = TRUE)),
      `null_q97.5` = unname(quantile(null_statistic, 0.975, na.rm = TRUE)),
      difference_median = median(observed_minus_null_statistic, na.rm = TRUE),
      `difference_q2.5` = unname(quantile(
        observed_minus_null_statistic,
        0.025,
        na.rm = TRUE)),
      `difference_q97.5` = unname(quantile(
        observed_minus_null_statistic,
        0.975,
        na.rm = TRUE)),
      empirical_p_upper = (sum(null_statistic >= observed_statistic, na.rm = TRUE) + 1) /
        (n_iterations + 1),
      observed_percentile = 100 * mean(null_statistic <= observed_statistic, na.rm = TRUE),
      .groups = "drop")
}

summarise_candidate_gene_vaf_null_tests <- function(
    candidate_gene_vaf_tbl,
    null_tbl,
    candidate_tbl,
    exclude_candidates_from_null) {
  observed_tbl <- candidate_gene_vaf_tbl %>%
    distinct(gene_symbol, recurrence_n_skin_samples, gene_level_mean_vaf) %>%
    arrange(gene_symbol)
  all_genes_null <- null_tbl %>%
    filter(
      as.character(candidate_group) == "All candidates",
      exclude_candidates_from_null ==
        .env$exclude_candidates_from_null) %>%
    pull(null_statistic)

  if (length(all_genes_null) == 0 || any(!is.finite(all_genes_null))) {
    stop("The All genes VAF violin null distribution is missing or non-finite.")
  }

  observed_tbl %>%
    mutate(
      `Empiric P-value` = vapply(gene_level_mean_vaf, function(observed_vaf) {
        (sum(all_genes_null >= observed_vaf) + 1) /
          (length(all_genes_null) + 1)
      }, numeric(1))) %>%
    left_join(
      candidate_tbl %>% select(gene_symbol, candidate_group),
      by = "gene_symbol") %>%
    arrange(candidate_group, gene_level_mean_vaf, gene_symbol) %>%
    transmute(
      Gene = gene_symbol,
      `Candidate group` = recode(
        as.character(candidate_group),
        `Initiation (I)` = "Initiation"),
      `Median VAF` = gene_level_mean_vaf,
      `Empiric P-value`)
}

plot_vaf_null_violin <- function(
    null_tbl,
    summary_tbl,
    exclude_candidates_from_null = VAF_NULL_EXCLUDE_CANDIDATES_FROM_ITERATIONS) {
  plot_group_levels <- c(
    "All candidates",
    "Predisposition-only",
    "Initiation-only")
  plot_group_labels <- c(
    "All candidates" = "All genes",
    "Predisposition-only" = "Predisposition",
    "Initiation-only" = "Initiation")
  main_null <- null_tbl %>%
    filter(exclude_candidates_from_null == .env$exclude_candidates_from_null) %>%
    mutate(candidate_group = factor(candidate_group, levels = plot_group_levels))
  observed_tbl <- summary_tbl %>%
    filter(exclude_candidates_from_null == .env$exclude_candidates_from_null) %>%
    mutate(candidate_group = factor(candidate_group, levels = plot_group_levels)) %>%
    distinct(candidate_group, observed_statistic)
  if (nrow(observed_tbl) != length(plot_group_levels) ||
      any(!is.finite(observed_tbl$observed_statistic))) {
    stop("Expected one finite observed VAF statistic per candidate group.")
  }
  candidate_higher_pct <- 100 * mean(
    main_null$observed_minus_null_statistic > 0,
    na.rm = TRUE)
  candidate_higher_label <- paste0(
    "Candidate genes > random selection in ",
    formatC(candidate_higher_pct, format = "f", digits = 1),
    " % iterations")
  candidate_higher_label <- str_wrap(candidate_higher_label, width = 30)
  violin_fill_values <- c(
    "All candidates" = "grey72",
    "Initiation-only" = "grey58",
    "Predisposition-only" = "grey84")

  ggplot(main_null, aes(x = candidate_group, y = null_statistic)) +
    geom_violin(
      aes(fill = candidate_group),
      colour = "grey25",
      linewidth = 0.25,
      trim = TRUE,
      scale = "width") +
    geom_boxplot(
      fill = "white",
      outlier.shape = NA) +
    geom_crossbar(
      data = observed_tbl,
      aes(
        x = candidate_group,
        y = observed_statistic,
        ymin = observed_statistic,
        ymax = observed_statistic),
      inherit.aes = FALSE,
      width = 0.6,
      colour = "red",
      linewidth = 0.5) +
    scale_fill_manual(values = violin_fill_values, guide = "none") +
    scale_x_discrete(drop = FALSE, labels = plot_group_labels) +
    scale_y_continuous(
      labels = vaf_label,
      limits = VAF_PLOT_LIMITS,
      breaks = VAF_PLOT_BREAKS,
      expand = expansion(mult = c(0.02, 0.16))) +
    labs(
      title = candidate_higher_label,
      x = NULL,
      y = "Median VAF") +
    my_theme +
    theme(
      text = element_text(size = 8),
      plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      strip.text = element_text(size = 8))
}

vaf_null_nonsyn_base <- build_vaf_null_nonsyn_base()
vaf_null_gene_vaf <- summarise_skin_gene_mean_vaf(
  evidence_tbl = vaf_null_nonsyn_base,
  sample_meta = skin_sample_meta)

#-------------------------------------------------------------------------------
### Skin VAF support histograms ###
#-------------------------------------------------------------------------------
summarise_gene_tissue_recurrence <- function(evidence_tbl) {
  evidence_tbl %>%
    filter(
      tissue %in% c("Skin", "Hair follicle"),
      context %in% CONTEXT_LEVELS,
      !is.na(SYMBOL),
      SYMBOL != "") %>%
    distinct(gene_symbol = SYMBOL, tissue, sample_uid) %>%
    group_by(gene_symbol) %>%
    summarise(
      skin_n_samples = n_distinct(sample_uid[tissue == "Skin"]),
      hf_n_samples = n_distinct(sample_uid[tissue == "Hair follicle"]),
      .groups = "drop")
}

make_vaf_support_pool <- function(gene_vaf_tbl, tissue_recurrence_tbl, pool_name) {
  pool_tbl <- gene_vaf_tbl %>%
    left_join(tissue_recurrence_tbl, by = "gene_symbol") %>%
    mutate(
      skin_n_samples = replace_na(skin_n_samples, 0L),
      hf_n_samples = replace_na(hf_n_samples, 0L))

  if (pool_name == "skin") {
    return(pool_tbl %>% filter(skin_n_samples >= 1))
  }

  if (pool_name == "skin_hf") {
    return(pool_tbl %>% filter(skin_n_samples >= 1, hf_n_samples >= 1))
  }

  stop("Unknown VAF support pool: ", pool_name)
}

draw_vaf_support_stat <- function(pool_tbl, recurrence_profile, set_size) {
  if (VAF_SUPPORT_MATCH_RECURRENCE) {
    sampled_vaf <- unlist(lapply(seq_len(nrow(recurrence_profile)), function(i) {
      recurrence_bin <- recurrence_profile$recurrence_n_skin_samples[[i]]
      n_to_draw <- recurrence_profile$n_to_draw[[i]]
      values <- pool_tbl %>%
        filter(recurrence_n_skin_samples == recurrence_bin) %>%
        pull(gene_level_mean_vaf)

      if (length(values) == 0) {
        stop("No genes available in recurrence bin ", recurrence_bin)
      }

      sample(values, n_to_draw, replace = length(values) < n_to_draw)
    }), use.names = FALSE)

    return(median(sampled_vaf, na.rm = TRUE))
  }

  values <- pool_tbl$gene_level_mean_vaf
  if (length(values) < set_size) {
    warning(
      "VAF support pool has fewer genes than requested; sampling with replacement.",
      call. = FALSE)
  }

  median(sample(values, set_size, replace = length(values) < set_size), na.rm = TRUE)
}

run_vaf_support_null <- function(
    pool_tbl,
    pool_label,
    candidate_tbl,
    gene_vaf_tbl,
    n_iterations,
    seed) {

  set_size <- n_distinct(candidate_tbl$gene_symbol)
  recurrence_profile <- candidate_tbl %>%
    distinct(gene_symbol) %>%
    left_join(
      gene_vaf_tbl %>%
        select(gene_symbol, recurrence_n_skin_samples),
      by = "gene_symbol") %>%
    count(recurrence_n_skin_samples, name = "n_to_draw")

  if (VAF_SUPPORT_MATCH_RECURRENCE) {
    missing_bins <- recurrence_profile %>%
      anti_join(
        pool_tbl %>% distinct(recurrence_n_skin_samples),
        by = "recurrence_n_skin_samples")

    if (nrow(missing_bins) > 0) {
      stop(
        pool_label,
        " is missing recurrence bins required by candidate genes: ",
        paste(missing_bins$recurrence_n_skin_samples, collapse = ", "))
    }
  }

  set.seed(seed)

  tibble(
    iteration = seq_len(n_iterations),
    pool = pool_label,
    match_recurrence = VAF_SUPPORT_MATCH_RECURRENCE,
    random_set_size = set_size,
    pool_size = n_distinct(pool_tbl$gene_symbol),
    null_statistic = replicate(
      n_iterations,
      draw_vaf_support_stat(pool_tbl, recurrence_profile, set_size)))
}

summarise_vaf_support_null <- function(null_tbl) {
  null_tbl %>%
    group_by(pool, match_recurrence) %>%
    summarise(
      n_iterations = n_distinct(iteration),
      random_set_size = dplyr::first(random_set_size),
      pool_size = dplyr::first(pool_size),
      null_mean = mean(null_statistic, na.rm = TRUE),
      null_median = median(null_statistic, na.rm = TRUE),
      null_q2.5 = unname(quantile(null_statistic, 0.025, na.rm = TRUE)),
      null_q97.5 = unname(quantile(null_statistic, 0.975, na.rm = TRUE)),
      .groups = "drop")
}

plot_skin_mutated_gene_vaf_histogram <- function(gene_vaf_tbl, support_summary) {
  line_tbl <- support_summary %>%
    transmute(
      pool,
      null_median)

  ggplot(gene_vaf_tbl, aes(x = gene_level_mean_vaf)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 40,
      fill = "grey70",
      colour = "white",
      linewidth = 0.2) +
    geom_vline(
      data = line_tbl,
      aes(xintercept = null_median, colour = pool),
      linewidth = 0.55) +
    scale_colour_manual(
      values = c(
        "Skin mutated genes" = "#D62728",
        "Skin and HF mutated genes" = "#7F3C8D"),
      name = "Null median") +
    scale_x_continuous(
      labels = vaf_label,
      breaks = VAF_PLOT_BREAKS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(xlim = VAF_PLOT_LIMITS) +
    labs(
      x = "Mean VAF per skin-mutated gene",
      y = "Density") +
    my_theme +
    theme(legend.position = "right")
}

vaf_support_tissue_recurrence <- summarise_gene_tissue_recurrence(vaf_null_nonsyn_base)
vaf_support_skin_pool <- make_vaf_support_pool(
  gene_vaf_tbl = vaf_null_gene_vaf,
  tissue_recurrence_tbl = vaf_support_tissue_recurrence,
  pool_name = "skin")
vaf_support_skin_hf_pool <- make_vaf_support_pool(
  gene_vaf_tbl = vaf_null_gene_vaf,
  tissue_recurrence_tbl = vaf_support_tissue_recurrence,
  pool_name = "skin_hf")

vaf_support_skin_null <- run_vaf_support_null(
  pool_tbl = vaf_support_skin_pool,
  pool_label = "Skin mutated genes",
  candidate_tbl = candidate_genes,
  gene_vaf_tbl = vaf_null_gene_vaf,
  n_iterations = VAF_NULL_ITERATIONS,
  seed = VAF_SUPPORT_SEED)
vaf_support_skin_hf_null <- run_vaf_support_null(
  pool_tbl = vaf_support_skin_hf_pool,
  pool_label = "Skin and HF mutated genes",
  candidate_tbl = candidate_genes,
  gene_vaf_tbl = vaf_null_gene_vaf,
  n_iterations = VAF_NULL_ITERATIONS,
  seed = VAF_SUPPORT_SEED + 1L)

vaf_support_null <- bind_rows(
  vaf_support_skin_null,
  vaf_support_skin_hf_null)
vaf_support_summary <- summarise_vaf_support_null(
  vaf_support_null)
skin_mutated_gene_vaf_histogram_data <- summarise_skin_gene_mean_vaf(
  evidence_tbl = nonsyn_base,
  sample_meta = skin_sample_meta)

skin_mutated_gene_vaf_histogram <- plot_skin_mutated_gene_vaf_histogram(
  gene_vaf_tbl = skin_mutated_gene_vaf_histogram_data,
  support_summary = vaf_support_summary)

export_ggplot_data(
  plot = skin_mutated_gene_vaf_histogram,
  data = skin_mutated_gene_vaf_histogram_data,
  file_name = file.path(OUTPUT_DIR, "Figure4_skin_mutated_gene_mean_VAF_histogram"),
  panel_width = export_plot_default_panel_width_mm,
  panel_height = export_plot_default_panel_height_mm)

candidate_vaf_threshold <- max(
  vaf_support_summary$null_median,
  na.rm = TRUE)

candidate_vaf_rule_tbl <- vaf_null_gene_vaf %>%
  select(gene_symbol) %>%
  left_join(
    skin_candidate_gene_vaf %>%
      transmute(
        gene_symbol,
        skin_mean_vaf_for_candidate_rule = mean_vaf),
    by = "gene_symbol") %>%
  filter(!is.na(skin_mean_vaf_for_candidate_rule)) %>%
  distinct(
    gene_symbol,
    skin_mean_vaf_for_candidate_rule)

candidate_genes_prefilter <- candidate_genes
candidate_genes <- candidate_genes %>%
  left_join(candidate_vaf_rule_tbl, by = "gene_symbol") %>%
  filter(skin_mean_vaf_for_candidate_rule > candidate_vaf_threshold) %>%
  mutate(
    candidate_vaf_threshold = candidate_vaf_threshold,
    candidate_group = factor(candidate_group, levels = CANDIDATE_GROUP_LEVELS))

if (nrow(candidate_genes) == 0) {
  stop("The support-histogram VAF threshold removed all candidate genes.")
}

candidate_vaf_rule_summary <- tibble(
  candidate_vaf_threshold = candidate_vaf_threshold,
  n_candidates_before_manual_exclusion = n_distinct(candidate_genes_prefilter$gene_symbol) +
    n_distinct(manually_excluded_candidate_genes$gene_symbol),
  n_manually_excluded_candidate_genes = n_distinct(
    manually_excluded_candidate_genes$gene_symbol),
  manually_excluded_candidate_genes = paste(
    manually_excluded_candidate_genes$gene_symbol,
    collapse = ";"),
  n_candidates_before_vaf_rule = n_distinct(candidate_genes_prefilter$gene_symbol),
  n_candidates_after_vaf_rule = n_distinct(candidate_genes$gene_symbol),
  vaf_support_match_recurrence = VAF_SUPPORT_MATCH_RECURRENCE,
  vaf_null_require_min_gt_af = VAF_NULL_REQUIRE_MIN_GT_AF)

write_csv(
  candidate_vaf_rule_summary,
  file.path(OUTPUT_DIR, "Figure4_candidate_gene_VAF_rule_summary_data.csv"))

skin_candidate_gene_vaf <- build_skin_gene_sample_vaf(
  gene_tbl = candidate_genes,
  evidence_tbl = nonsyn_base,
  sample_meta = skin_sample_meta) %>%
  summarise_skin_gene_vaf()

skin_gene_vaf <- skin_candidate_gene_vaf %>%
  set_lollipop_order()

vaf_null_candidate_observed_vaf <- skin_candidate_gene_vaf %>%
  transmute(
    gene_symbol,
    recurrence_n_skin_samples = n_samples_mutated,
    gene_level_mean_vaf = mean_vaf)

vaf_null_candidate_sets <- build_vaf_null_candidate_sets(candidate_genes)
missing_vaf_null_candidates <- setdiff(
  unique(vaf_null_candidate_sets$gene_symbol),
  vaf_null_gene_vaf$gene_symbol)
missing_vaf_null_observed_candidates <- setdiff(
  unique(vaf_null_candidate_sets$gene_symbol),
  vaf_null_candidate_observed_vaf$gene_symbol)

if (length(missing_vaf_null_candidates) > 0) {
  stop(
    "Some candidate genes are absent from the skin VAF null table: ",
    paste(missing_vaf_null_candidates, collapse = ", "))
}

if (length(missing_vaf_null_observed_candidates) > 0) {
  stop(
    "Some candidate genes are absent from the filtered candidate skin VAF table: ",
    paste(missing_vaf_null_observed_candidates, collapse = ", "))
}

vaf_null_candidate_summary <- summarise_vaf_null_candidates(
  candidate_sets = vaf_null_candidate_sets,
  gene_vaf_tbl = vaf_null_candidate_observed_vaf)
vaf_null_recurrence_pool_summary <- summarise_vaf_recurrence_pools(
  gene_vaf_tbl = vaf_null_gene_vaf,
  candidate_tbl = candidate_genes,
  candidate_recurrence_tbl = vaf_null_candidate_observed_vaf)
vaf_null_long <- bind_rows(
  run_vaf_recurrence_matched_null(
    candidate_sets = vaf_null_candidate_sets,
    gene_vaf_tbl = vaf_null_gene_vaf,
    n_iterations = VAF_NULL_ITERATIONS,
    seed = VAF_NULL_SEED,
    exclude_candidates_from_null = TRUE,
    candidate_tbl = candidate_genes,
    candidate_recurrence_tbl = vaf_null_candidate_observed_vaf),
  run_vaf_recurrence_matched_null(
    candidate_sets = vaf_null_candidate_sets,
    gene_vaf_tbl = vaf_null_gene_vaf,
    n_iterations = VAF_NULL_ITERATIONS,
    seed = VAF_NULL_SEED,
    exclude_candidates_from_null = FALSE,
    candidate_tbl = candidate_genes,
    candidate_recurrence_tbl = vaf_null_candidate_observed_vaf)) %>%
  add_vaf_null_observed_difference(
    observed_tbl = vaf_null_candidate_summary)
vaf_null_summary <- summarise_vaf_null_test(
  null_tbl = vaf_null_long,
  observed_tbl = vaf_null_candidate_summary)
vaf_null_candidate_gene_table <- summarise_candidate_gene_vaf_null_tests(
  candidate_gene_vaf_tbl = vaf_null_candidate_observed_vaf,
  null_tbl = vaf_null_long,
  candidate_tbl = candidate_genes,
  exclude_candidates_from_null =
    VAF_NULL_EXCLUDE_CANDIDATES_FROM_ITERATIONS)

if (nrow(vaf_null_candidate_gene_table) != n_distinct(candidate_genes$gene_symbol) ||
    anyDuplicated(vaf_null_candidate_gene_table$Gene) > 0) {
  stop("Candidate-gene VAF null table does not contain one row per candidate gene.")
}

vaf_null_violin_plot <- plot_vaf_null_violin(
  null_tbl = vaf_null_long,
  summary_tbl = vaf_null_summary,
  exclude_candidates_from_null = VAF_NULL_EXCLUDE_CANDIDATES_FROM_ITERATIONS)

write_csv(
  vaf_null_candidate_summary,
  file.path(OUTPUT_DIR, "Figure4_skin_gene_VAF_null_candidate_summary_data.csv"))
write_csv(
  vaf_null_recurrence_pool_summary,
  file.path(OUTPUT_DIR, "Figure4_skin_gene_VAF_null_recurrence_pool_summary_data.csv"))
write_csv(
  vaf_null_summary,
  file.path(OUTPUT_DIR, "Figure4_skin_gene_VAF_null_randomization_summary_data.csv"))
write_csv(
  vaf_null_candidate_gene_table,
  file.path(
    OUTPUT_DIR,
    "Figure4_skin_gene_VAF_null_violin_candidate_gene_table.csv"))
write_csv(
  vaf_null_long,
  file.path(OUTPUT_DIR, "Figure4_skin_gene_VAF_null_randomization_long_data.csv"))

vaf_null_violin_panel_width_mm <- 30
vaf_null_violin_plot
ggsave(
  filename = file.path(OUTPUT_DIR, "Figure4_skin_gene_VAF_null_violin.pdf"),
  plot = set_plot_panel_size(
    vaf_null_violin_plot,
    panel_width = vaf_null_violin_panel_width_mm,
    panel_height = export_plot_default_panel_height_mm),
  width = export_plot_default_dimension_mm(
    1,
    vaf_null_violin_panel_width_mm,
    export_plot_default_tick_size_mm) + 35,
  height = export_plot_default_dimension_mm(
    1,
    export_plot_default_panel_height_mm,
    export_plot_default_tick_size_mm) + 30,
  units = "mm")

#-------------------------------------------------------------------------------
### Candidate gene matrices ###
#-------------------------------------------------------------------------------
candidate_gene_matrix <- make_gene_matrix(
  gene_tbl = candidate_genes,
  context_tbl = nonsyn_base) %>%
  format_candidate_matrix_contexts()

order_matrix_like_lollipop <- function(matrix_tbl, lollipop_tbl) {
  row_order_levels <- unique(as.character(lollipop_tbl$gene_label))

  matrix_tbl %>%
    mutate(
      lollipop_order = match(as.character(row_label), row_order_levels),
      row_label = factor(
        as.character(row_label),
        levels = rev(row_order_levels)),
      candidate_group = factor(
        as.character(candidate_group),
        levels = levels(lollipop_tbl$candidate_group)),
      mean_n_samples = as.numeric(n_samples)) %>%
    filter(!is.na(row_label)) %>%
    arrange(lollipop_order, context)
}

candidate_gene_matrix_ordered <- candidate_gene_matrix %>%
  order_matrix_like_lollipop(skin_gene_vaf)

plot_candidate_matrix <- function(plot_df) {
  plot_df <- drop_empty_contexts(plot_df, "n_samples")

  ggplot(plot_df, aes(x = context, y = row_label)) +
    geom_tile(fill = "white", colour = "grey86", linewidth = 0.25) +
    geom_point(
      data = plot_df %>% filter(detected),
      aes(size = n_samples, fill = first_time),
      shape = 21,
      stroke = 0.55) +
    facet_grid(
      candidate_group ~ .,
      scales = "free_y",
      space = "free_y") +
    scale_fill_manual(
      values = col_palette$time,
      drop = FALSE,
      name = "Earliest evidence") +
    scale_size_continuous(
      range = sample_count_size_range,
      limits = sample_count_size_limits,
      breaks = sample_count_size_breaks,
      name = "# samples") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(labels = italic_gene_axis_labels) +
    labs(x = "", y = "") +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      strip.text.y = element_text(angle = 0),
      legend.position = "right")
}

candidate_genes_plot <- plot_candidate_matrix(candidate_gene_matrix_ordered)
skin_gene_vaf_export_height_mm <- max(
  long_gene_panel_height_mm(skin_gene_vaf$gene_symbol),
  long_gene_panel_height_mm(candidate_gene_matrix_ordered$gene_symbol))

export_ggplot_data(
  plot = candidate_genes_plot,
  data = candidate_gene_matrix_ordered,
  file_name = file.path(OUTPUT_DIR, "Figure4_candidate_genes"),
  height = skin_gene_vaf_export_height_mm,
  panel_width = export_plot_default_panel_width_mm)

#-------------------------------------------------------------------------------
### Skin gene VAF lollipops ###
#-------------------------------------------------------------------------------
plot_skin_gene_lollipop <- function(plot_df) {
  vaf_col <- active_vaf_col()
  null_line_tbl <- vaf_support_summary %>%
    transmute(
      pool,
      null_median)

  ggplot(plot_df, aes(x = .data[[vaf_col]], y = gene_label)) +
    geom_vline(
      data = null_line_tbl,
      aes(xintercept = null_median, colour = pool),
      linetype = "longdash",
      linewidth = 0.4) +
    geom_segment(
      aes(x = 0, xend = .data[[vaf_col]], yend = gene_label),
      colour = "grey58",
      linewidth = 0.45) +
    geom_point(
      aes(size = n_samples_mutated, fill = candidate_group),
      shape = 21,
      colour = "black",
      stroke = 0.45) +
    facet_grid(
      candidate_group ~ .,
      scales = "free_y",
      space = "free_y") +
    scale_fill_manual(
      values = candidate_fill_values,
      drop = FALSE,
      name = "Candidate group") +
    scale_colour_manual(
      values = c(
        "Skin mutated genes" = "#D62728",
        "Skin and HF mutated genes" = "#7F3C8D"),
      name = "Null median") +
    scale_size_continuous(
      range = sample_count_size_range,
      limits = sample_count_size_limits,
      breaks = sample_count_size_breaks,
      name = "# mutated skin samples") +
    scale_x_continuous(
      labels = vaf_label,
      breaks = VAF_PLOT_BREAKS,
      expand = expansion(mult = c(0, 0.04))) +
    scale_y_discrete(labels = italic_gene_axis_labels) +
    coord_cartesian(xlim = VAF_PLOT_LIMITS) +
    labs(x = active_vaf_label(), y = "") +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      strip.text.y = element_text(angle = 0),
      legend.position = "right")
}

skin_gene_vaf_plot <- plot_skin_gene_lollipop(skin_gene_vaf)

export_ggplot_data(
  plot = skin_gene_vaf_plot,
  data = skin_gene_vaf,
  file_name = file.path(OUTPUT_DIR, "Figure4_skin_gene_VAF"),
  height = skin_gene_vaf_export_height_mm,
  panel_width = export_plot_default_panel_width_mm)

#-------------------------------------------------------------------------------
### cBioPortal cSCC oncoplots ###
#-------------------------------------------------------------------------------
first_existing_name <- function(tbl, candidates) {
  hits <- intersect(candidates, names(tbl))

  if (length(hits) == 0) NA_character_ else hits[[1]]
}

pull_existing_column <- function(tbl, candidates, default = NA_character_) {
  hit <- first_existing_name(tbl, candidates)

  if (is.na(hit)) {
    return(rep(default, nrow(tbl)))
  }

  tbl[[hit]]
}

standardise_cbio_mutations <- function(tbl, study_id) {
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(tibble(
      study_id = character(),
      sample_id = character(),
      patient_id = character(),
      hugo_symbol = character(),
      mutation_type = character(),
      protein_change = character(),
      molecular_profile_id = character(),
      cbio_sample_uid = character()))
  }

  tibble(
    study_id = coalesce(
      as.character(pull_existing_column(tbl, c("studyId", "study_id"))),
      study_id),
    sample_id = as.character(pull_existing_column(
      tbl,
      c("sampleId", "sample_id", "Tumor_Sample_Barcode"))),
    patient_id = as.character(pull_existing_column(
      tbl,
      c("patientId", "patient_id"))),
    hugo_symbol = as.character(pull_existing_column(
      tbl,
      c(
        "hugoGeneSymbol",
        "HugoGeneSymbol",
        "Hugo_Symbol",
        "gene.hugoGeneSymbol"))),
    mutation_type = as.character(pull_existing_column(
      tbl,
      c("mutationType", "Variant_Classification", "variantClassification"))),
    protein_change = as.character(pull_existing_column(
      tbl,
      c("proteinChange", "Protein_Change", "HGVSp_Short"))),
    molecular_profile_id = as.character(pull_existing_column(
      tbl,
      c("molecularProfileId", "molecular_profile_id")))) %>%
    mutate(cbio_sample_uid = paste(study_id, sample_id, sep = "::")) %>%
    filter(
      !is.na(sample_id),
      sample_id != "",
      !is.na(hugo_symbol),
      hugo_symbol != "")
}

cbio_mutation_class <- function(mutation_type) {
  mutation_type <- str_to_lower(coalesce(as.character(mutation_type), ""))

  case_when(
    str_detect(mutation_type, "missense") ~ "Missense",
    str_detect(mutation_type, "in[_ ]?frame|inframe") ~ "Inframe",
    str_detect(
      mutation_type,
      "frame|nonsense|nonstop|nonstart|splice|start|stop") ~
      "Truncating/splice",
    TRUE ~ "Other")
}

load_cbio_available_profiles <- function(study_id) {
  if ("available_profiles" %in% getNamespaceExports("cbioportalR")) {
    return(cbioportalR::available_profiles(study_id = study_id))
  }

  cbioportalR::available_molecular_profiles(study_id = study_id)
}

pick_cbio_mutation_profile_id <- function(study_id) {
  profiles <- load_cbio_available_profiles(study_id)
  profile_id_col <- first_existing_name(
    profiles,
    c("molecularProfileId", "molecular_profile_id", "profileId", "profile_id"))

  mutation_profiles <- profiles %>%
    mutate(
      profile_text = str_to_lower(paste(
        pull_existing_column(
          .,
          c(
            "molecularAlterationType",
            "molecular_alteration_type",
            "datatype",
            "profileType",
            "profile_type"),
          ""),
        pull_existing_column(., c("name"), ""),
        pull_existing_column(., c("description"), ""),
        sep = " "))) %>%
    filter(
      str_detect(profile_text, "mutation"),
      !str_detect(profile_text, "copy|cna|fusion|structural|rna|mrna"))

  if (nrow(mutation_profiles) == 0) {
    mutation_profiles <- profiles %>%
      filter(str_detect(as.character(.data[[profile_id_col]]), "mutations?$"))
  }

  if (nrow(mutation_profiles) == 0) {
    stop("Could not resolve the mutation molecular profile for ", study_id, ".")
  }

  mutation_profiles %>%
    arrange(as.character(.data[[profile_id_col]])) %>%
    dplyr::slice(1) %>%
    pull(all_of(profile_id_col)) %>%
    as.character()
}

get_cbio_mutation_profile_info <- function(study_id) {
  profiles <- load_cbio_available_profiles(study_id)
  profile_id_col <- first_existing_name(
    profiles,
    c("molecularProfileId", "molecular_profile_id", "profileId", "profile_id"))
  mutation_profile_id <- pick_cbio_mutation_profile_id(study_id)
  panel_col <- first_existing_name(
    profiles,
    c("genePanelId", "gene_panel_id", "genePanel", "gene_panel"))

  profiles %>%
    filter(as.character(.data[[profile_id_col]]) == mutation_profile_id) %>%
    dplyr::slice(1) %>%
    transmute(
      study_id = study_id,
      molecular_profile_id = mutation_profile_id,
      profile_gene_panel = if (is.na(panel_col)) {
        NA_character_
      } else {
        na_if(as.character(.data[[panel_col]]), "NA")
      }) %>%
    mutate(
      profile_gene_panel = if_else(
        profile_gene_panel == "",
        NA_character_,
        profile_gene_panel))
}

call_cbio_mutations_by_sample <- function(sample_study_pairs, gene_symbols) {
  mutation_fun <- cbioportalR::get_mutations_by_sample
  mutation_formals <- names(formals(mutation_fun))
  args <- list(sample_study_pairs = sample_study_pairs)

  if ("genes" %in% mutation_formals) {
    args$genes <- gene_symbols
  } else if ("hugo_gene_symbols" %in% mutation_formals) {
    args$hugo_gene_symbols <- gene_symbols
  } else if ("gene_list" %in% mutation_formals) {
    args$gene_list <- gene_symbols
  } else {
    stop(
      "cbioportalR::get_mutations_by_sample() does not expose a supported ",
      "gene-symbol argument.")
  }

  if ("add_hugo" %in% mutation_formals) {
    args$add_hugo <- TRUE
  }

  do.call(mutation_fun, args)
}

safe_cbio_mutations_by_sample <- function(sample_study_pairs, gene_symbols, study_id) {
  tryCatch(
    call_cbio_mutations_by_sample(sample_study_pairs, gene_symbols),
    error = function(err) {
      if (length(gene_symbols) > 1) {
        return(bind_rows(lapply(gene_symbols, function(gene_symbol) {
          safe_cbio_mutations_by_sample(
            sample_study_pairs,
            gene_symbol,
            study_id)
        })))
      }

      warning(
        "cBioPortal mutation query failed for ",
        study_id,
        " gene ",
        gene_symbols[[1]],
        ": ",
        conditionMessage(err),
        call. = FALSE)
      tibble()
    })
}

load_cbio_profiled_samples <- function(study_ids) {
  bind_rows(lapply(study_ids, function(study_id) {
    profile_info <- get_cbio_mutation_profile_info(study_id)
    sample_lists <- cbioportalR::available_sample_lists(study_id = study_id)
    sample_list_id_col <- first_existing_name(
      sample_lists,
      c("sampleListId", "sample_list_id"))

    mutation_lists <- sample_lists %>%
      mutate(
        list_text = str_to_lower(paste(
          pull_existing_column(., c("category"), ""),
          pull_existing_column(., c("name"), ""),
          pull_existing_column(., c("description"), ""),
          sep = " "))) %>%
      filter(str_detect(list_text, "mutation"))

    if (nrow(mutation_lists) == 0) {
      stop("Could not resolve mutation-profiled samples for ", study_id, ".")
    }

    sample_list_id <- mutation_lists %>%
      mutate(
        preferred_rank = if_else(
          str_to_lower(pull_existing_column(., c("category"), "")) ==
            "all_cases_with_mutation_data",
          0L,
          1L)) %>%
      arrange(preferred_rank) %>%
      dplyr::slice(1) %>%
      pull(all_of(sample_list_id_col))

    samples <- cbioportalR::available_samples(sample_list_id = sample_list_id)
    sample_id_col <- first_existing_name(samples, c("sampleId", "sample_id"))

    tibble(
      study_id = study_id,
      sample_id = as.character(samples[[sample_id_col]]),
      sample_list_id = sample_list_id,
      molecular_profile_id = profile_info$molecular_profile_id[[1]],
      profile_gene_panel = profile_info$profile_gene_panel[[1]],
      cbio_sample_uid = paste(study_id, sample_id, sep = "::"))
  })) %>%
    distinct()
}

cbio_api_post_json <- function(path, body) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The httr and jsonlite packages are required for this cBioPortal query.")
  }

  response <- httr::POST(
    url = paste0(CBIO_API_BASE_URL, path),
    httr::add_headers(
      Accept = "application/json",
      `Content-Type` = "application/json"),
    body = jsonlite::toJSON(body, auto_unbox = TRUE, null = "null"),
    encode = "raw")

  if (httr::http_error(response)) {
    stop(
      "cBioPortal API request failed: ",
      httr::status_code(response),
      " ",
      httr::content(response, as = "text", encoding = "UTF-8"))
  }

  content_text <- httr::content(response, as = "text", encoding = "UTF-8")
  if (identical(content_text, "") || is.na(content_text)) {
    return(tibble())
  }

  jsonlite::fromJSON(content_text, flatten = TRUE) %>%
    as_tibble()
}

query_ensembl_human_gene_lengths <- function(human_symbols) {
  human_symbols <- sort(unique(na.omit(human_symbols)))
  human_symbols <- human_symbols[human_symbols != ""]

  empty_tbl <- tibble(human_symbol = human_symbols) %>%
    mutate(
      human_gene_length_bp = NA_integer_,
      human_gene_assembly = NA_character_,
      human_gene_seq_region = NA_character_,
      human_gene_start = NA_integer_,
      human_gene_end = NA_integer_,
      human_gene_length_source = NA_character_)

  if (length(human_symbols) == 0) {
    return(empty_tbl)
  }

  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    warning(
      "httr and jsonlite are required for Ensembl human gene length lookup; ",
      "human_gene_length_bp will be NA.",
      call. = FALSE)
    return(empty_tbl)
  }

  bind_rows(lapply(human_symbols, function(human_symbol) {
    tryCatch({
      response <- httr::GET(
        url = paste0(
          ENSEMBL_REST_BASE_URL,
          "/lookup/symbol/homo_sapiens/",
          utils::URLencode(human_symbol, reserved = TRUE)),
        httr::add_headers(
          Accept = "application/json",
          `Content-Type` = "application/json"))

      if (httr::status_code(response) == 404) {
        warning(
          "No Ensembl human gene length found for ",
          human_symbol,
          ".",
          call. = FALSE)
        return(empty_tbl %>% filter(human_symbol == .env$human_symbol))
      }

      if (httr::http_error(response)) {
        stop(
          "Ensembl REST request failed: ",
          httr::status_code(response),
          " ",
          httr::content(response, as = "text", encoding = "UTF-8"))
      }

      gene_info <- jsonlite::fromJSON(
        httr::content(response, as = "text", encoding = "UTF-8"),
        flatten = TRUE)

      gene_start <- as.integer(gene_info$start)
      gene_end <- as.integer(gene_info$end)

      tibble(
        human_symbol = human_symbol,
        human_gene_length_bp = abs(gene_end - gene_start) + 1L,
        human_gene_assembly = as.character(gene_info$assembly_name),
        human_gene_seq_region = as.character(gene_info$seq_region_name),
        human_gene_start = gene_start,
        human_gene_end = gene_end,
        human_gene_length_source = "Ensembl REST lookup/symbol")
    }, error = function(err) {
      warning(
        "Ensembl human gene length lookup failed for ",
        human_symbol,
        ": ",
        conditionMessage(err),
        call. = FALSE)
      empty_tbl %>% filter(human_symbol == .env$human_symbol)
    })
  })) %>%
    distinct(human_symbol, .keep_all = TRUE)
}

summarise_human_gene_lengths_for_mouse_genes <- function(gene_tbl, human_gene_lengths) {
  gene_tbl %>%
    select(gene_symbol, human_symbols) %>%
    mutate(human_symbols = coalesce(human_symbols, "")) %>%
    separate_rows(human_symbols, sep = ";") %>%
    transmute(
      gene_symbol,
      human_symbol = trimws(human_symbols)) %>%
    filter(human_symbol != "") %>%
    left_join(human_gene_lengths, by = "human_symbol") %>%
    group_by(gene_symbol) %>%
    summarise(
      human_gene_length_bp = if_else(
        n_distinct(human_symbol) == 1L,
        dplyr::first(human_gene_length_bp),
        NA_integer_),
      human_gene_lengths_bp = paste(
        paste0(
          human_symbol,
          "=",
          if_else(
            is.na(human_gene_length_bp),
            "NA",
            as.character(human_gene_length_bp))),
        collapse = ";"),
      human_gene_length_source = collapse_unique(
        human_gene_length_source,
        n = Inf,
        sep = ";"),
      human_gene_assembly = collapse_unique(
        human_gene_assembly,
        n = Inf,
        sep = ";"),
      .groups = "drop")
}

normalise_cbio_sample_panels <- function(panel_tbl, profiled_samples) {
  if (is.null(panel_tbl) || nrow(panel_tbl) == 0) {
    return(profiled_samples %>%
      mutate(
        gene_panel = profile_gene_panel,
        gene_profiled = TRUE,
        gene_panel_status = if_else(
          is.na(gene_panel) | gene_panel == "",
          "whole_profile_no_panel_matrix",
          "profile_gene_panel_fallback")))
  }

  sample_id_col <- first_existing_name(panel_tbl, c("sampleId", "sample_id"))
  study_id_col <- first_existing_name(panel_tbl, c("studyId", "study_id"))
  panel_col <- first_existing_name(
    panel_tbl,
    c("genePanel", "gene_panel", "genePanelId", "gene_panel_id"))
  profiled_col <- first_existing_name(panel_tbl, c("profiled", "isProfiled"))

  if (is.na(sample_id_col) || is.na(study_id_col) || is.na(panel_col)) {
    return(profiled_samples %>%
      mutate(
        gene_panel = profile_gene_panel,
        gene_profiled = TRUE,
        gene_panel_status = if_else(
          is.na(gene_panel) | gene_panel == "",
          "whole_profile_panel_lookup_unparsed",
          "profile_gene_panel_fallback")))
  }

  panel_tbl %>%
    transmute(
      study_id = as.character(.data[[study_id_col]]),
      sample_id = as.character(.data[[sample_id_col]]),
      sample_gene_panel = na_if(as.character(.data[[panel_col]]), "NA"),
      gene_profiled = if (is.na(profiled_col)) {
        TRUE
      } else {
        as.logical(.data[[profiled_col]])
      }) %>%
    right_join(profiled_samples, by = c("study_id", "sample_id")) %>%
    mutate(
      gene_panel = coalesce(sample_gene_panel, profile_gene_panel),
      gene_profiled = replace_na(gene_profiled, TRUE),
      gene_panel_status = if_else(
        is.na(gene_panel) | gene_panel == "",
        "whole_profile_no_gene_panel_reported",
        if_else(
          is.na(sample_gene_panel) | sample_gene_panel == "",
          "profile_gene_panel_fallback",
          "gene_panel_reported")),
      cbio_sample_uid = paste(study_id, sample_id, sep = "::")) %>%
    distinct()
}

load_cbio_sample_gene_panels <- function(profiled_samples) {
  direct_panel_data <- tryCatch({
    sample_pairs <- profiled_samples %>%
      transmute(
        molecularProfileId = molecular_profile_id,
        sampleId = sample_id) %>%
      distinct()

    bind_rows(lapply(
      split(sample_pairs, ceiling(seq_len(nrow(sample_pairs)) / 250)),
      function(sample_chunk) {
        cbio_api_post_json(
          path = "/gene-panel-data/fetch?projection=DETAILED",
          body = list(sampleMolecularIdentifiers = sample_chunk))
      }))
  }, error = function(err) tibble())

  if (nrow(direct_panel_data) > 0) {
    return(normalise_cbio_sample_panels(
      direct_panel_data,
      profiled_samples) %>%
        mutate(
          gene_panel_status = case_when(
            is.na(gene_panel) | gene_panel == "" ~
              "whole_profile_gene_panel_data_direct",
            TRUE ~ "gene_panel_data_direct")))
  }

  if (!"get_panel_by_sample" %in% getNamespaceExports("cbioportalR")) {
    return(profiled_samples %>%
      mutate(
        gene_panel = profile_gene_panel,
        gene_profiled = TRUE,
        gene_panel_status = if_else(
          is.na(gene_panel) | gene_panel == "",
          "whole_profile_get_panel_by_sample_unavailable",
          "profile_gene_panel_fallback")))
  }

  sample_study_pairs <- profiled_samples %>%
    transmute(study_id, sample_id, molecular_profile_id) %>%
    distinct()

  tryCatch(
    normalise_cbio_sample_panels(
      cbioportalR::get_panel_by_sample(sample_study_pairs = sample_study_pairs),
      profiled_samples),
    error = function(err) {
      profiled_samples %>%
        mutate(
          gene_panel = profile_gene_panel,
          gene_profiled = TRUE,
          gene_panel_status = if_else(
            is.na(gene_panel) | gene_panel == "",
            "whole_profile_panel_lookup_failed",
            "profile_gene_panel_fallback"))
    })
}

extract_cbio_panel_genes <- function(panel_tbl, panel_id) {
  symbol_col <- first_existing_name(
    panel_tbl,
    c("hugoGeneSymbol", "hugo_gene_symbol", "geneSymbol", "gene_symbol", "symbol"))

  if (!is.na(symbol_col)) {
    return(panel_tbl %>%
      transmute(
        panel_id = panel_id,
        human_symbol = as.character(.data[[symbol_col]])) %>%
      filter(!is.na(human_symbol), human_symbol != "") %>%
      distinct())
  }

  gene_list_col <- first_existing_name(panel_tbl, c("geneList", "gene_list", "genes"))
  if (is.na(gene_list_col)) {
    return(tibble(panel_id = character(), human_symbol = character()))
  }

  panel_tbl %>%
    transmute(
      panel_id = panel_id,
      human_symbol = str_split(as.character(.data[[gene_list_col]]), "\\s+|;|,")) %>%
    unnest(human_symbol) %>%
    mutate(human_symbol = trimws(human_symbol)) %>%
    filter(!is.na(human_symbol), human_symbol != "") %>%
    distinct()
}

load_cbio_panel_gene_members <- function(panel_ids, query_human_symbols) {
  panel_ids <- sort(unique(na.omit(panel_ids)))

  if (length(panel_ids) == 0 ||
      !"get_gene_panel" %in% getNamespaceExports("cbioportalR")) {
    return(tibble(panel_id = character(), human_symbol = character()))
  }

  bind_rows(lapply(panel_ids, function(panel_id) {
    tryCatch(
      extract_cbio_panel_genes(
        cbioportalR::get_gene_panel(panel_id = panel_id),
        panel_id),
      error = function(err) tibble(panel_id = character(), human_symbol = character()))
  })) %>%
    filter(human_symbol %in% query_human_symbols) %>%
    distinct()
}

build_cbio_mouse_gene_map <- function(gene_symbols, orthologue_tbl) {
  tibble(gene_symbol = unique(gene_symbols)) %>%
    left_join(orthologue_tbl, by = "gene_symbol") %>%
    mutate(human_symbols = coalesce(human_symbols, "")) %>%
    split_human_symbols()
}

build_cbio_profiled_sample_map <- function(
    mouse_gene_map,
    sample_gene_panels,
    panel_gene_members) {

  query_human_symbols <- sort(unique(mouse_gene_map$human_symbol))

  whole_profile_samples <- sample_gene_panels %>%
    filter(
      gene_profiled,
      is.na(gene_panel) | gene_panel == "",
      CBIO_ASSUME_NO_PANEL_MEANS_WHOLE_PROFILE |
        str_detect(gene_panel_status, "^whole_profile")) %>%
    select(study_id, sample_id, cbio_sample_uid) %>%
    distinct()

  whole_profile_map <- whole_profile_samples %>%
    crossing(human_symbol = query_human_symbols)

  panel_profile_map <- sample_gene_panels %>%
    filter(gene_profiled, !is.na(gene_panel), gene_panel != "") %>%
    inner_join(
      panel_gene_members,
      by = c("gene_panel" = "panel_id"),
      relationship = "many-to-many") %>%
    select(study_id, sample_id, cbio_sample_uid, human_symbol) %>%
    distinct()

  bind_rows(whole_profile_map, panel_profile_map) %>%
    inner_join(
      mouse_gene_map,
      by = "human_symbol",
      relationship = "many-to-many") %>%
    select(gene_symbol, human_symbol, study_id, sample_id, cbio_sample_uid) %>%
    distinct()
}

load_cbio_mutations <- function(profiled_samples, query_human_symbols) {
  if ("set_cbioportal_db" %in% getNamespaceExports("cbioportalR")) {
    try(cbioportalR::set_cbioportal_db("public"), silent = TRUE)
  }

  query_human_symbols <- sort(unique(na.omit(query_human_symbols)))

  bind_rows(lapply(split(profiled_samples, profiled_samples$study_id), function(samples) {
    study_id <- unique(samples$study_id)[[1]]
    message("Loading cBioPortal mutations for ", study_id)

    sample_study_pairs <- samples %>%
      transmute(study_id, sample_id, molecular_profile_id) %>%
      distinct()

    gene_chunks <- split(
      query_human_symbols,
      ceiling(seq_along(query_human_symbols) / 50))

    bind_rows(lapply(gene_chunks, function(gene_chunk) {
      safe_cbio_mutations_by_sample(
        sample_study_pairs = sample_study_pairs,
        gene_symbols = gene_chunk,
        study_id = study_id) %>%
        standardise_cbio_mutations(study_id)
    }))
  })) %>%
    mutate(mutation_class = cbio_mutation_class(mutation_type))
}

load_cbio_all_study_mutations <- function(
    profiled_samples,
    query_human_symbols) {

  if ("set_cbioportal_db" %in% getNamespaceExports("cbioportalR")) {
    try(cbioportalR::set_cbioportal_db("public"), silent = TRUE)
  }

  query_human_symbols <- sort(unique(na.omit(query_human_symbols)))
  profiled_sample_uids <- profiled_samples %>%
    select(study_id, sample_id, cbio_sample_uid) %>%
    distinct()

  bind_rows(lapply(
    split(profiled_samples, profiled_samples$study_id),
    function(samples) {
      study_id <- unique(samples$study_id)[[1]]
      molecular_profile_id <- unique(samples$molecular_profile_id)[[1]]
      sample_list_id <- unique(samples$sample_list_id)[[1]]
      message("Loading all cBioPortal mutations for ", study_id)

      study_mutations <- tryCatch(
        cbio_api_post_json(
          path = paste0(
            "/molecular-profiles/",
            molecular_profile_id,
            "/mutations/fetch?projection=DETAILED&pageSize=10000000"),
          body = list(sampleListId = sample_list_id)),
        error = function(err) {
          stop(
            "Could not load all mutations for ",
            study_id,
            ": ",
            conditionMessage(err),
            call. = FALSE)
        })

      standardise_cbio_mutations(study_mutations, study_id)
    })) %>%
    semi_join(
      profiled_sample_uids,
      by = c("study_id", "sample_id", "cbio_sample_uid")) %>%
    filter(hugo_symbol %in% query_human_symbols) %>%
    mutate(mutation_class = cbio_mutation_class(mutation_type))
}

summarise_cbio_occurrence <- function(mouse_gene_map, cbio_mutations, profiled_map) {
  denominators <- profiled_map %>%
    distinct(gene_symbol, cbio_sample_uid) %>%
    count(gene_symbol, name = "n_profiled_samples") %>%
    right_join(
      tibble(gene_symbol = unique(mouse_gene_map$gene_symbol)),
      by = "gene_symbol") %>%
    mutate(n_profiled_samples = replace_na(n_profiled_samples, 0L))

  altered_samples <- cbio_mutations %>%
    inner_join(
      mouse_gene_map,
      by = c("hugo_symbol" = "human_symbol"),
      relationship = "many-to-many") %>%
    semi_join(
      profiled_map,
      by = c("gene_symbol", "hugo_symbol" = "human_symbol", "cbio_sample_uid")) %>%
    distinct(gene_symbol, cbio_sample_uid) %>%
    count(gene_symbol, name = "n_altered_samples")

  denominators %>%
    left_join(altered_samples, by = "gene_symbol") %>%
    mutate(
      n_altered_samples = replace_na(n_altered_samples, 0L),
      occurrence_pct = if_else(
        n_profiled_samples > 0,
        100 * n_altered_samples / n_profiled_samples,
        NA_real_))
}

load_cbio_mutated_gene_incidence <- function(file_name) {
  read_tsv(file_name, show_col_types = FALSE) %>%
    transmute(
      human_symbol = Gene,
      n_mutations = `# Mut`,
      n_altered_samples = `#`,
      n_profiled_samples = `Profiled Samples`,
      incidence_pct = parse_number(Freq),
      is_cancer_gene = `Is Cancer Gene (source: OncoKB)`) %>%
    filter(!is.na(human_symbol), human_symbol != "", !is.na(incidence_pct)) %>%
    distinct(human_symbol, .keep_all = TRUE)
}

summarise_cbio_candidate_incidence <- function(candidate_tbl, incidence_tbl) {
  candidate_human_symbols <- candidate_tbl %>%
    select(gene_symbol, gene_label, candidate_group, human_symbols) %>%
    mutate(human_symbols = coalesce(human_symbols, "")) %>%
    separate_rows(human_symbols, sep = ";") %>%
    transmute(
      gene_symbol,
      gene_label,
      candidate_group,
      human_symbol = trimws(human_symbols)) %>%
    filter(human_symbol != "")

  missing_symbols <- candidate_human_symbols %>%
    anti_join(incidence_tbl, by = "human_symbol")

  if (nrow(missing_symbols) > 0) {
    warning(
      "Some candidate human orthologues were absent from ",
      basename(CBIO_MUTATED_GENES_FILE),
      ": ",
      paste(sort(unique(missing_symbols$human_symbol)), collapse = ", "),
      call. = FALSE)
  }

  candidate_human_symbols %>%
    inner_join(incidence_tbl, by = "human_symbol") %>%
    group_by(gene_symbol, gene_label, candidate_group) %>%
    summarise(
      matched_human_symbols = paste(sort(unique(human_symbol)), collapse = ";"),
      incidence_pct = mean(incidence_pct, na.rm = TRUE),
      max_incidence_pct = max(incidence_pct, na.rm = TRUE),
      n_matched_human_symbols = n_distinct(human_symbol),
      .groups = "drop") %>%
    mutate(
      candidate_group_plot = factor(
        candidate_group_display(candidate_group),
        levels = CBIO_CANDIDATE_GROUP_LEVELS))
}

run_cbio_incidence_permutation <- function(
    incidence_tbl,
    candidate_incidence_tbl,
    n_iterations,
    seed) {

  set.seed(seed)

  incidence_values <- incidence_tbl$incidence_pct
  n_candidates <- n_distinct(candidate_incidence_tbl$gene_symbol)
  observed_median <- median(candidate_incidence_tbl$incidence_pct, na.rm = TRUE)

  null_medians <- vapply(seq_len(n_iterations), function(iteration) {
    median(sample(incidence_values, n_candidates, replace = FALSE), na.rm = TRUE)
  }, numeric(1))

  permutation_long <- tibble(
    iteration = seq_len(n_iterations),
    n_candidate_genes = n_candidates,
    median_incidence_pct = null_medians)

  permutation_summary <- tibble(
    n_iterations = n_iterations,
    n_candidate_genes = n_candidates,
    universe_genes = nrow(incidence_tbl),
    observed_median_incidence_pct = observed_median,
    null_mean = mean(null_medians, na.rm = TRUE),
    null_median = median(null_medians, na.rm = TRUE),
    null_q2.5 = unname(quantile(null_medians, 0.025, na.rm = TRUE)),
    null_q97.5 = unname(quantile(null_medians, 0.975, na.rm = TRUE)),
    empirical_p_upper = (sum(null_medians >= observed_median, na.rm = TRUE) + 1) /
      (n_iterations + 1),
    empirical_p_lower = (sum(null_medians <= observed_median, na.rm = TRUE) + 1) /
      (n_iterations + 1),
    empirical_p_two_sided = pmin(1, 2 * pmin(empirical_p_upper, empirical_p_lower)))

  list(
    summary = permutation_summary,
    long = permutation_long)
}

calculate_cbio_candidate_incidence_empirical_p_values <- function(
    candidate_incidence_tbl,
    incidence_tbl) {

  incidence_values <- incidence_tbl$incidence_pct
  incidence_values <- incidence_values[!is.na(incidence_values)]
  n_universe_genes <- length(incidence_values)

  candidate_incidence_tbl %>%
    mutate(
      incidence_empirical_p_upper = vapply(
        incidence_pct,
        function(observed_incidence) {
          (sum(incidence_values >= observed_incidence, na.rm = TRUE) + 1) /
            (n_universe_genes + 1)
        },
        numeric(1)),
      incidence_empirical_p_lower = vapply(
        incidence_pct,
        function(observed_incidence) {
          (sum(incidence_values <= observed_incidence, na.rm = TRUE) + 1) /
            (n_universe_genes + 1)
        },
        numeric(1)),
      incidence_empirical_p_two_sided = pmin(
        1,
        2 * pmin(incidence_empirical_p_upper, incidence_empirical_p_lower)),
      incidence_neg_log10_empirical_p = -log10(incidence_empirical_p_upper),
      n_universe_genes = n_universe_genes) %>%
    arrange(incidence_empirical_p_upper, desc(incidence_pct), gene_label) %>%
    mutate(incidence_empirical_rank = row_number())
}

plot_cbio_incidence_permutation_histogram <- function(
    permutation_long,
    permutation_summary) {

  histogram_counts <- hist(
    permutation_long$median_incidence_pct,
    breaks = 40,
    plot = FALSE)$counts
  y_tick_max <- max(200, ceiling(max(histogram_counts, na.rm = TRUE) / 1000) * 1000)

  label_df <- permutation_summary %>%
    mutate(
      label = paste0(
        "Median=",
        formatC(observed_median_incidence_pct, format = "f", digits = 1),
        "%\n",
        format_p_value(empirical_p_upper)))

  ggplot(permutation_long, aes(x = median_incidence_pct)) +
    geom_histogram(
      bins = 40,
      fill = "grey70",
      colour = "white",
      linewidth = 0.2) +
    geom_vline(
      data = permutation_summary,
      aes(xintercept = observed_median_incidence_pct),
      colour = "#D62728",
      linewidth = 0.5) +
    geom_text(
      data = label_df,
      aes(x = observed_median_incidence_pct, y = Inf, label = label),
      inherit.aes = FALSE,
      colour = "#D62728",
      hjust = 1.05,
      vjust = 1.2,
      size = 8 / .pt) +
    scale_x_continuous(
      breaks = c(0, 15),
      labels = function(x) formatC(x, format = "f", digits = 0)) +
    scale_y_continuous(
      breaks = c(0, y_tick_max),
      labels = function(x) formatC(x, format = "f", digits = 0),
      limits = c(0, y_tick_max),
      expand = expansion(mult = c(0.04, 0.12))) +
    coord_cartesian(xlim = c(0, 15), clip = "off") +
    labs(
      x = "Incidence (%)",
      y = "Permutation iterations") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.position = "none")
}

cbio_histogram_breaks <- function(x, bins = 40L) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(seq(0, 1, length.out = bins + 1L))
  }

  x_range <- range(x)
  if (diff(x_range) == 0) {
    bin_width <- 1
  } else {
    bin_width <- diff(x_range) / (bins - 1L)
  }

  seq(
    x_range[[1]] - bin_width / 2,
    x_range[[2]] + bin_width / 2,
    length.out = bins + 1L)
}

prepare_cbio_candidate_incidence_density <- function(
    permutation_long,
    candidate_incidence_tbl,
    bins = 40L,
    density_n = 512L) {

  x_axis_max <- ceiling(max(
    permutation_long$median_incidence_pct,
    candidate_incidence_tbl$incidence_pct,
    na.rm = TRUE) / 5) * 5

  candidate_density <- density(
    candidate_incidence_tbl$incidence_pct,
    from = 0,
    to = x_axis_max,
    n = density_n,
    na.rm = TRUE)
  candidate_density_max <- max(candidate_density$y, na.rm = TRUE)

  histogram_breaks <- cbio_histogram_breaks(
    permutation_long$median_incidence_pct,
    bins = bins)
  histogram_counts <- hist(
    permutation_long$median_incidence_pct,
    breaks = histogram_breaks,
    plot = FALSE)$counts
  max_histogram_count <- max(histogram_counts, na.rm = TRUE)
  y_tick_max <- if_else(
    max_histogram_count <= 1800 && candidate_density_max <= 0.04,
    1800,
    max(200, ceiling(max_histogram_count / 1000) * 1000))
  density_axis_max <- if_else(
    max_histogram_count <= 1800 && candidate_density_max <= 0.04,
    0.04,
    ceiling(candidate_density_max * 100) / 100)

  density_scale <- if_else(
    density_axis_max > 0,
    y_tick_max / density_axis_max,
    1)

  tibble(
    incidence_pct = candidate_density$x,
    candidate_density = candidate_density$y,
    density_scale = density_scale,
    density_scaled = candidate_density * density_scale,
    x_axis_max = x_axis_max,
    y_tick_max = y_tick_max,
    density_axis_max = density_axis_max)
}

plot_cbio_incidence_permutation_histogram_with_candidate_density <- function(
    permutation_long,
    permutation_summary,
    candidate_density_tbl) {

  x_axis_max <- dplyr::first(candidate_density_tbl$x_axis_max)
  y_tick_max <- dplyr::first(candidate_density_tbl$y_tick_max)
  density_scale <- dplyr::first(candidate_density_tbl$density_scale)
  density_axis_max <- dplyr::first(candidate_density_tbl$density_axis_max)
  histogram_breaks <- cbio_histogram_breaks(permutation_long$median_incidence_pct)

  label_df <- permutation_summary %>%
    mutate(
      median_y = approx(
        x = candidate_density_tbl$incidence_pct,
        y = candidate_density_tbl$density_scaled,
        xout = observed_median_incidence_pct,
        rule = 2)$y,
      label_x = observed_median_incidence_pct + 8,
      label_y = 0.9 * y_tick_max,
      label = format_p_value(empirical_p_upper))

  ggplot(permutation_long, aes(x = median_incidence_pct)) +
    geom_area(
      data = candidate_density_tbl,
      aes(x = incidence_pct, y = density_scaled),
      inherit.aes = FALSE,
      fill = "#BC667E",
      alpha = 0.35,
      colour = "#8C3F55",
      linewidth = 0.25) +
    geom_histogram(
      breaks = histogram_breaks,
      fill = "grey55",
      colour = "grey55",
      linewidth = 0.1) +
    geom_segment(
      data = label_df,
      aes(
        x = observed_median_incidence_pct,
        xend = observed_median_incidence_pct,
        y = 0,
        yend = median_y),
      inherit.aes = FALSE,
      colour = "#D62728",
      linewidth = 0.5) +
    geom_text(
      data = label_df,
      aes(x = label_x, y = label_y, label = label),
      inherit.aes = FALSE,
      colour = "#D62728",
      hjust = 0.5,
      vjust = 0.5,
      size = 8 / .pt) +
    scale_x_continuous(
      breaks = c(0, 30, 60),
      labels = function(x) formatC(x, format = "f", digits = 0),
      expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      breaks = c(0, y_tick_max),
      labels = function(x) formatC(x, format = "f", digits = 0),
      limits = c(0, y_tick_max),
      expand = expansion(mult = c(0, 0.08)),
      sec.axis = sec_axis(
        ~ . / density_scale,
        name = "Candidate gene density",
        breaks = c(0, density_axis_max),
        labels = function(x) formatC(x, format = "f", digits = 2))) +
    coord_cartesian(xlim = c(0, x_axis_max), clip = "off") +
    labs(
      x = "Incidence (%)",
      y = "Permutation iterations") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text.y.right = element_text(size = 8),
      axis.title.y.right = element_text(size = 8),
      legend.position = "none")
}

plot_cbio_candidate_p_values <- function(
    candidate_p_value_tbl,
    x_col,
    x_label,
    x_breaks,
    x_labels = waiver()) {

  plot_df <- candidate_p_value_tbl %>%
    filter(!is.na(.data[[x_col]]), !is.na(incidence_neg_log10_empirical_p))

  label_df <- plot_df %>%
    filter(incidence_empirical_p_upper <= 0.05) %>%
    arrange(incidence_empirical_p_upper, desc(.data[[x_col]]), gene_label)

  label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene_symbol),
      colour = "black",
      size = 7 / .pt,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.size = 0.2,
      force = 3,
      force_pull = 0.05,
      max.time = 2,
      max.iter = 10000,
      max.overlaps = Inf,
      seed = 1)
  } else {
    geom_text(
      data = label_df,
      aes(label = gene_symbol),
      colour = "black",
      size = 7 / .pt,
      hjust = -0.08,
      vjust = 0.5)
  }

  y_max <- ceiling(max(plot_df$incidence_neg_log10_empirical_p, na.rm = TRUE))

  ggplot(
    plot_df,
    aes(x = .data[[x_col]], y = incidence_neg_log10_empirical_p)) +
    geom_hline(
      yintercept = -log10(0.05),
      colour = "grey55",
      linewidth = 0.25,
      linetype = "dashed") +
    geom_point(
      aes(fill = candidate_group_plot),
      shape = 21,
      colour = "black",
      alpha = 0.85,
      stroke = 0.35,
      size = 2.2) +
    label_layer +
    scale_fill_manual(values = cbio_candidate_fill_values, drop = FALSE) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0.04, 0.2))) +
    scale_y_continuous(
      breaks = 0:y_max,
      labels = function(x) formatC(x, format = "f", digits = 0),
      expand = expansion(mult = c(0.08, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(
      x = x_label,
      y = expression(-log[10]("Empirical P value")),
      fill = "Candidate group") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.position = "right")
}

plot_cbio_candidate_incidence_p_values <- function(candidate_p_value_tbl) {
  plot_cbio_candidate_p_values(
    candidate_p_value_tbl = candidate_p_value_tbl,
    x_col = "incidence_pct",
    x_label = "Incidence (%)",
    x_breaks = c(0, 30, 60),
    x_labels = function(x) formatC(x, format = "f", digits = 0))
}

build_cbio_oncoplot_hits <- function(mouse_gene_map, cbio_mutations) {
  mutation_priority <- c("Truncating/splice", "Inframe", "Missense", "Other")

  cbio_mutations %>%
    inner_join(
      mouse_gene_map,
      by = c("hugo_symbol" = "human_symbol"),
      relationship = "many-to-many") %>%
    mutate(mutation_class = factor(mutation_class, levels = mutation_priority)) %>%
    group_by(gene_symbol, cbio_sample_uid) %>%
    summarise(
      human_symbols_mutated = paste(sort(unique(hugo_symbol)), collapse = ";"),
      mutation_classes = paste(sort(unique(as.character(mutation_class))), collapse = ";"),
      mutation_class = as.character(
        mutation_class[which.min(as.integer(mutation_class))]),
      n_cell_mutations = n(),
      .groups = "drop")
}

build_cbio_extra_gene_table <- function(gene_symbols, orthologue_tbl) {
  tibble(
    gene_symbol = gene_symbols,
    candidate_group = CBIO_EXTRA_ONCOPLOT_GROUP) %>%
    annotate_human_orthologues(orthologue_tbl) %>%
    select(
      gene_symbol,
      gene_label,
      candidate_group,
      human_symbols,
      has_human_orthologue)
}

build_cbio_gene_labels <- function(
    gene_tbl,
    cbio_occurrence,
    human_gene_lengths) {
  length_summary <- summarise_human_gene_lengths_for_mouse_genes(
    gene_tbl,
    human_gene_lengths)

  gene_tbl %>%
    select(
      gene_symbol,
      gene_label,
      candidate_group,
      human_symbols,
      has_human_orthologue) %>%
    filter(has_human_orthologue) %>%
    left_join(cbio_occurrence, by = "gene_symbol") %>%
    left_join(length_summary, by = "gene_symbol") %>%
    mutate(
      candidate_group_plot = factor(
        candidate_group_display(candidate_group),
        levels = CBIO_ONCOPLOT_GROUP_LEVELS),
      cbio_human_label = str_replace_all(human_symbols, ";", "/"),
      cbio_row_label = paste0(
        gene_symbol,
        " / ",
        cbio_human_label,
        " (",
        if_else(
          is.na(occurrence_pct),
          "not profiled",
          paste0(formatC(occurrence_pct, format = "f", digits = 1), "%")),
        ")"))
}

build_cbio_oncoplot_df <- function(gene_labels, cbio_hits, profiled_map, sample_order) {
  expand_grid(
    gene_symbol = gene_labels$gene_symbol,
    cbio_sample_uid = sample_order$cbio_sample_uid) %>%
    left_join(sample_order, by = "cbio_sample_uid") %>%
    left_join(
      profiled_map %>%
        distinct(gene_symbol, cbio_sample_uid) %>%
        mutate(is_gene_profiled = TRUE),
      by = c("gene_symbol", "cbio_sample_uid")) %>%
    left_join(cbio_hits, by = c("gene_symbol", "cbio_sample_uid")) %>%
    left_join(gene_labels, by = "gene_symbol") %>%
    mutate(
      cbio_sample_uid = factor(
        cbio_sample_uid,
        levels = sample_order$cbio_sample_uid),
      is_gene_profiled = replace_na(is_gene_profiled, FALSE),
      mutation_class_plot = case_when(
        !is_gene_profiled ~ "Not profiled",
        !is.na(mutation_class) ~ mutation_class,
        TRUE ~ "Not altered"),
      n_cell_mutations = replace_na(n_cell_mutations, 0L))
}

plot_cbio_oncoplot <- function(plot_df, row_levels) {
  plot_df <- plot_df %>%
    mutate(
      cbio_row_label = factor(cbio_row_label, levels = row_levels),
      mutation_class_plot = factor(
        mutation_class_plot,
        levels = c(
          "Truncating/splice",
          "Inframe",
          "Missense",
          "Other",
          "Not altered",
          "Not profiled"))) %>%
    filter(!is.na(cbio_row_label))

  ggplot(plot_df, aes(x = cbio_sample_uid, y = cbio_row_label)) +
    geom_tile(
      aes(fill = mutation_class_plot),
      colour = "white",
      linewidth = 0.2,
      height = 0.85) +
    geom_point(
      data = plot_df %>% filter(mutation_class_plot == "Not profiled"),
      shape = 95,
      colour = "grey55",
      size = 1.6,
      stroke = 0.25) +
    facet_grid(
      candidate_group_plot ~ study_label,
      scales = "free",
      space = "free") +
    scale_fill_manual(
      values = c(
        "Truncating/splice" = "#000000",
        Inframe = "#8B5A2B",
        Missense = "#3D8B37",
        Other = "#7EA3CF",
        "Not altered" = "#D9D9D9",
        "Not profiled" = "#FFFFFF"),
      breaks = c("Truncating/splice", "Inframe", "Missense", "Other"),
      drop = FALSE) +
    labs(x = "", y = "", fill = "") +
    my_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text.y = element_text(angle = 0),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.35, "cm"),
      panel.spacing.x = grid::unit(0.08, "lines")) +
    guides(fill = guide_legend(nrow = 1))
}

cbio_mutated_gene_incidence <- load_cbio_mutated_gene_incidence(
  CBIO_MUTATED_GENES_FILE)
cbio_candidate_incidence <- summarise_cbio_candidate_incidence(
  candidate_tbl = candidate_genes,
  incidence_tbl = cbio_mutated_gene_incidence)

if (nrow(cbio_candidate_incidence) == 0) {
  stop("No candidate genes could be matched to cBioPortal_Mutated_Genes.txt.")
}

cbio_incidence_permutation <- run_cbio_incidence_permutation(
  incidence_tbl = cbio_mutated_gene_incidence,
  candidate_incidence_tbl = cbio_candidate_incidence,
  n_iterations = CBIO_INCIDENCE_PERMUTATIONS,
  seed = CBIO_INCIDENCE_PERMUTATION_SEED)
cbio_incidence_permutation_summary <- cbio_incidence_permutation$summary
cbio_incidence_permutation_long <- cbio_incidence_permutation$long
cbio_candidate_incidence_p_values <- calculate_cbio_candidate_incidence_empirical_p_values(
  candidate_incidence_tbl = cbio_candidate_incidence,
  incidence_tbl = cbio_mutated_gene_incidence)

cbio_incidence_permutation_histogram <- plot_cbio_incidence_permutation_histogram(
  permutation_long = cbio_incidence_permutation_long,
  permutation_summary = cbio_incidence_permutation_summary)
cbio_incidence_permutation_candidate_density <- prepare_cbio_candidate_incidence_density(
  permutation_long = cbio_incidence_permutation_long,
  candidate_incidence_tbl = cbio_candidate_incidence)
cbio_incidence_permutation_histogram_candidate_density <- plot_cbio_incidence_permutation_histogram_with_candidate_density(
  permutation_long = cbio_incidence_permutation_long,
  permutation_summary = cbio_incidence_permutation_summary,
  candidate_density_tbl = cbio_incidence_permutation_candidate_density)
cbio_candidate_incidence_p_value_plot <- plot_cbio_candidate_incidence_p_values(
  cbio_candidate_incidence_p_values)

write_csv(
  cbio_mutated_gene_incidence,
  file.path(OUTPUT_DIR, "Figure4_cbioportal_mutated_gene_incidence_data.csv"))
write_csv(
  cbio_candidate_incidence,
  file.path(OUTPUT_DIR, "Figure4_cbioportal_candidate_incidence_data.csv"))
write_csv(
  cbio_incidence_permutation_summary,
  file.path(OUTPUT_DIR, "Figure4_cbioportal_incidence_permutation_summary_data.csv"))
write_csv(
  cbio_incidence_permutation_long,
  file.path(OUTPUT_DIR, "Figure4_cbioportal_incidence_permutation_long_data.csv"))
write_csv(
  cbio_candidate_incidence_p_values,
  file.path(OUTPUT_DIR, "Figure4_cbioportal_candidate_incidence_empirical_p_values_data.csv"))

ggsave(
  filename = file.path(OUTPUT_DIR, "Figure4_cbioportal_incidence_permutation_histogram.pdf"),
  plot = set_plot_panel_size(
    cbio_incidence_permutation_histogram,
    panel_width = export_plot_default_panel_width_mm,
    panel_height = export_plot_default_panel_height_mm),
  width = export_plot_default_dimension_mm(
    1,
    export_plot_default_panel_width_mm,
    export_plot_default_tick_size_mm),
  height = export_plot_default_dimension_mm(
    1,
    export_plot_default_panel_height_mm,
    export_plot_default_tick_size_mm),
  units = "mm")
ggplot2::set_last_plot(cbio_incidence_permutation_histogram_candidate_density)
export_plot_data(
  data = cbio_incidence_permutation_candidate_density,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_incidence_permutation_histogram_candidate_density"))
ggplot2::set_last_plot(cbio_candidate_incidence_p_value_plot)
export_plot_data(
  data = cbio_candidate_incidence_p_values,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_incidence_empirical_p_values"))

if (identical(EXTERNAL_DATA_MODE, "frozen")) {
  publication_copy_frozen(
    file.path(FIGURE4_SOURCE_DIR,
      "Figure4_cbioportal_candidate_oncoplot_data.csv"),
    file.path(FIGURE4_SOURCE_DIR,
      "Figure4_cbioportal_candidate_oncoplot_data.csv"),
    "eb63bd4ee5a87c8b0ca06e782a8e2482634f72671764851644e20a65066c6ea4")
  publication_copy_frozen(
    file.path(FIGURE4_SOURCE_DIR,
      "Figure4_cbioportal_candidate_oncoplot.pdf"),
    file.path(FIGURE4_SOURCE_DIR,
      "Figure4_cbioportal_candidate_oncoplot.pdf"),
    "4315fd9a2425861fa4cd8be3fb1e8ffcc83ba767e8694660327f0a1a5caccec1")
} else {
cbio_profiled_samples <- load_cbio_profiled_samples(CBIO_STUDY_IDS)
cbio_sample_order <- cbio_profiled_samples %>%
  mutate(
    study_id = factor(study_id, levels = CBIO_STUDY_IDS),
    study_label = as.character(study_id)) %>%
  arrange(study_id, sample_id) %>%
  mutate(sample_rank = row_number())

cbio_extra_oncoplot_genes <- build_cbio_extra_gene_table(
  gene_symbols = CBIO_EXTRA_ONCOPLOT_GENES,
  orthologue_tbl = mouse_human_orthologues)
cbio_oncoplot_gene_tbl <- bind_rows(
  candidate_genes %>%
    filter(!gene_symbol %in% CBIO_EXTRA_ONCOPLOT_GENES) %>%
    select(
      gene_symbol,
      gene_label,
      candidate_group,
      human_symbols,
      has_human_orthologue),
  cbio_extra_oncoplot_genes)

cbio_all_mouse_gene_map <- build_cbio_mouse_gene_map(
  unique(cbio_oncoplot_gene_tbl$gene_symbol),
  mouse_human_orthologues)

cbio_query_human_symbols <- sort(unique(cbio_all_mouse_gene_map$human_symbol))
cbio_human_gene_lengths <- query_ensembl_human_gene_lengths(
  cbio_query_human_symbols)
cbio_sample_gene_panels <- load_cbio_sample_gene_panels(cbio_profiled_samples)
cbio_panel_gene_members <- load_cbio_panel_gene_members(
  panel_ids = cbio_sample_gene_panels$gene_panel,
  query_human_symbols = cbio_query_human_symbols)
cbio_profiled_sample_map <- build_cbio_profiled_sample_map(
  mouse_gene_map = cbio_all_mouse_gene_map,
  sample_gene_panels = cbio_sample_gene_panels,
  panel_gene_members = cbio_panel_gene_members)

if (nrow(cbio_profiled_sample_map) == 0) {
  stop("No gene-specific cBioPortal denominators were resolved.")
}

cbio_mutations <- load_cbio_mutations(
  profiled_samples = cbio_profiled_samples,
  query_human_symbols = cbio_query_human_symbols)

cbio_mutations_for_occurrence <- cbio_mutations %>%
  filter(
    !CBIO_RESTRICT_TO_NONSYNONYMOUS |
      mutation_class %in% CBIO_NONSYNONYMOUS_CLASSES)

cbio_occurrence <- summarise_cbio_occurrence(
  mouse_gene_map = cbio_all_mouse_gene_map,
  cbio_mutations = cbio_mutations_for_occurrence,
  profiled_map = cbio_profiled_sample_map)

cbio_hits <- build_cbio_oncoplot_hits(
  mouse_gene_map = cbio_all_mouse_gene_map,
  cbio_mutations = cbio_mutations_for_occurrence)

cbio_candidate_labels <- build_cbio_gene_labels(
  gene_tbl = cbio_oncoplot_gene_tbl,
  cbio_occurrence = cbio_occurrence,
  human_gene_lengths = cbio_human_gene_lengths)

cbio_candidate_row_levels <- tibble(gene_label = levels(skin_gene_vaf$gene_label)) %>%
  left_join(
    cbio_candidate_labels %>%
      select(gene_label, cbio_row_label),
    by = "gene_label") %>%
  filter(!is.na(cbio_row_label)) %>%
  pull(cbio_row_label) %>%
  c(
    cbio_candidate_labels %>%
      filter(gene_symbol %in% CBIO_EXTRA_ONCOPLOT_GENES) %>%
      arrange(match(gene_symbol, CBIO_EXTRA_ONCOPLOT_GENES)) %>%
      pull(cbio_row_label)) %>%
  unique()

cbio_candidate_oncoplot_df <- build_cbio_oncoplot_df(
  gene_labels = cbio_candidate_labels,
  cbio_hits = cbio_hits,
  profiled_map = cbio_profiled_sample_map,
  sample_order = cbio_sample_order) %>%
  mutate(
    cbio_sample_uid = factor(
      cbio_sample_uid,
      levels = cbio_sample_order$cbio_sample_uid))

cbio_candidate_oncoplot <- plot_cbio_oncoplot(
  cbio_candidate_oncoplot_df,
  row_levels = cbio_candidate_row_levels)

export_ggplot_data(
  plot = cbio_candidate_oncoplot,
  data = cbio_candidate_oncoplot_df,
  file_name = file.path(OUTPUT_DIR, "Figure4_cbioportal_candidate_oncoplot"),
  width = max(260, min(660, 150 + 3 * nrow(cbio_sample_order))),
  height = max(
    130,
    min(400, 12 * n_distinct(cbio_candidate_oncoplot_df$cbio_row_label) + 50)))

#-------------------------------------------------------------------------------
### cBioPortal clinical associations ###
#-------------------------------------------------------------------------------
CBIO_CLINICAL_ENDPOINT_DEFINITIONS <- tribble(
  ~endpoint_id,
  ~analysis_family,
  ~analysis_scope,
  ~endpoint_label,
  ~study_ids,
  ~clinical_attribute,
  ~reference_group,
  ~case_group,
  ~endpoint_order,
  "sample_type_pooled",
  "Primary / metastatic",
  "Pooled primary analysis",
  "Primary vs metastatic (all cohorts)",
  paste(CBIO_STUDY_IDS, collapse = ";"),
  "SAMPLE_TYPE plus study design",
  "Primary",
  "Metastasis",
  1L,
  "sample_type_explicit",
  "Primary / metastatic",
  "Explicit-field sensitivity",
  "Primary vs metastatic (explicit field)",
  "cscc_hgsc_bcm_2014;cscc_ranson_2022",
  "SAMPLE_TYPE",
  "Primary",
  "Metastasis",
  2L,
  "sample_type_hgsc",
  "Primary / metastatic",
  "Within-study sensitivity",
  "Primary vs metastatic (HGSC)",
  "cscc_hgsc_bcm_2014",
  "SAMPLE_TYPE",
  "Primary",
  "Metastasis",
  3L,
  "hgsc_t_stage",
  "Disease aggressiveness",
  "Within-study",
  "HGSC T stage",
  "cscc_hgsc_bcm_2014",
  "T_STAGE",
  "T1/T2",
  "T3/T4",
  4L,
  "hgsc_n_stage",
  "Disease aggressiveness",
  "Within-study",
  "HGSC N stage",
  "cscc_hgsc_bcm_2014",
  "N_STAGE",
  "N0",
  "N1",
  5L,
  "hgsc_m_stage",
  "Disease aggressiveness",
  "Within-study",
  "HGSC M stage",
  "cscc_hgsc_bcm_2014",
  "M_STAGE",
  "M0",
  "M1",
  6L,
  "hgsc_recurrence_persistence",
  "Disease aggressiveness",
  "Within-study",
  "HGSC recurrence / persistence",
  "cscc_hgsc_bcm_2014",
  "DISEASE_RECURRENCE_OR_PERSISTENCE",
  "No",
  "Yes",
  7L,
  "hgsc_efs_status",
  "Disease aggressiveness",
  "Within-study",
  "HGSC event-free status",
  "cscc_hgsc_bcm_2014",
  "EFS_STATUS",
  "Censored",
  "Event",
  8L,
  "ranson_nodal_stage",
  "Disease aggressiveness",
  "Within-study",
  "Ranson nodal stage",
  "cscc_ranson_2022",
  "NODAL_STAGE",
  "N1/N2",
  "N3",
  9L,
  "ranson_grade",
  "Disease aggressiveness",
  "Within-study",
  "Ranson grade",
  "cscc_ranson_2022",
  "GRADE",
  "Grade 1/2",
  "Grade 3",
  10L,
  "ranson_extracapsular_spread",
  "Disease aggressiveness",
  "Within-study",
  "Ranson extracapsular spread",
  "cscc_ranson_2022",
  "EXTRACAPSULAR_SPREAD",
  "No",
  "Yes",
  11L,
  "dfci_dfs_status",
  "Disease aggressiveness",
  "Within-study",
  "DFCI disease-free status",
  "cscc_dfarber_2015",
  "DFS_STATUS",
  "Disease free",
  "Recurred / progressed",
  12L,
  "recurrence_pooled",
  "Recurrence / progression",
  "Pooled study-adjusted primary analysis",
  "Disease-free vs recurrent / progressed (DFCI + HGSC)",
  "cscc_dfarber_2015;cscc_hgsc_bcm_2014",
  "DFS_STATUS or DISEASE_RECURRENCE_OR_PERSISTENCE",
  "Disease free",
  "Recurred / progressed",
  13L)

CBIO_FOCUSED_PRIMARY_METASTATIC_GENES <- c("Kndc1", "Foxp4")
CBIO_FOCUSED_SENSITIVITY_DEFINITIONS <- tribble(
  ~endpoint_id,
  ~analysis_scope,
  ~endpoint_label,
  ~source_endpoint_id,
  ~excluded_study_id,
  ~included_study_ids,
  ~endpoint_order,
  "focused_all_cohorts",
  "All cohorts",
  "All cohorts",
  "sample_type_pooled",
  NA_character_,
  paste(CBIO_STUDY_IDS, collapse = ";"),
  1L,
  "focused_explicit_field",
  "Explicit SAMPLE_TYPE field",
  "Explicit SAMPLE_TYPE",
  "sample_type_explicit",
  NA_character_,
  "cscc_hgsc_bcm_2014;cscc_ranson_2022",
  2L,
  "focused_hgsc_only",
  "HGSC within-study",
  "HGSC only",
  "sample_type_hgsc",
  NA_character_,
  "cscc_hgsc_bcm_2014",
  3L,
  "focused_leave_out_dfci",
  "Leave-one-cohort-out",
  "All except DFCI",
  "sample_type_pooled",
  "cscc_dfarber_2015",
  paste(setdiff(CBIO_STUDY_IDS, "cscc_dfarber_2015"), collapse = ";"),
  4L,
  "focused_leave_out_hgsc",
  "Leave-one-cohort-out",
  "All except HGSC",
  "sample_type_pooled",
  "cscc_hgsc_bcm_2014",
  paste(setdiff(CBIO_STUDY_IDS, "cscc_hgsc_bcm_2014"), collapse = ";"),
  5L,
  "focused_leave_out_ranson",
  "Leave-one-cohort-out",
  "All except Ranson",
  "sample_type_pooled",
  "cscc_ranson_2022",
  paste(setdiff(CBIO_STUDY_IDS, "cscc_ranson_2022"), collapse = ";"),
  6L,
  "focused_leave_out_ucsf",
  "Leave-one-cohort-out",
  "All except UCSF",
  "sample_type_pooled",
  "cscc_ucsf_2021",
  paste(setdiff(CBIO_STUDY_IDS, "cscc_ucsf_2021"), collapse = ";"),
  7L)
CBIO_PRIMARY_METASTATIC_COHORT_GROUPS <- tribble(
  ~comparison_cohort, ~endpoint_group, ~study_id,
  "DFCI", "Metastasis", "cscc_dfarber_2015",
  "HGSC", "Primary", "cscc_hgsc_bcm_2014",
  "HGSC", "Metastasis", "cscc_hgsc_bcm_2014",
  "Ranson", "Metastasis", "cscc_ranson_2022",
  "UCSF", "Primary", "cscc_ucsf_2021",
  "Pooled", "Primary", paste(CBIO_STUDY_IDS, collapse = ";"),
  "Pooled", "Metastasis", paste(CBIO_STUDY_IDS, collapse = ";"))

CBIO_TREATMENT_RESPONSE_AUDIT <- tibble(
  analysis_family = "Treatment response",
  direct_response_available = FALSE,
  analysis_performed = FALSE,
  available_related_attributes = paste(
    c(
      "DFCI: RADIATION_THERAPY, DFS_STATUS, PFS_MONTHS",
      "HGSC: PRIOR_TREATMENT, RADIATION_OR_CHEMOTHERAPY, EFS_STATUS"),
    collapse = "; "),
  interpretation = paste(
    "No cSCC study supplies a direct response, responder, or RECIST endpoint.",
    "Treatment exposure and recurrence/progression are not relabelled as",
    "treatment response."))

load_cbio_clinical_data <- function(study_ids) {
  bind_rows(lapply(study_ids, function(study_id) {
    message("Loading cBioPortal clinical data for ", study_id)

    clinical_tbl <- cbioportalR::get_clinical_by_study(study_id = study_id)

    clinical_tbl %>%
      transmute(
        study_id = coalesce(
          as.character(pull_existing_column(., c("studyId", "study_id"))),
          study_id),
        patient_id = as.character(pull_existing_column(
          .,
          c("patientId", "patient_id"))),
        sample_id = as.character(pull_existing_column(
          .,
          c("sampleId", "sample_id"))),
        clinical_attribute_id = as.character(pull_existing_column(
          .,
          c("clinicalAttributeId", "clinical_attribute_id"))),
        clinical_value = as.character(pull_existing_column(
          .,
          c("value", "clinical_value"))),
        clinical_level = str_to_upper(as.character(pull_existing_column(
          .,
          c("dataLevel", "data_level")))))
  })) %>%
    mutate(
      sample_id = na_if(sample_id, ""),
      patient_id = na_if(patient_id, ""),
      clinical_value = na_if(trimws(clinical_value), "")) %>%
    filter(
      !is.na(study_id),
      !is.na(patient_id),
      !is.na(clinical_attribute_id),
      clinical_attribute_id != "",
      clinical_level %in% c("SAMPLE", "PATIENT")) %>%
    distinct()
}

build_cbio_sample_clinical_data <- function(clinical_long) {
  sample_map <- clinical_long %>%
    filter(clinical_level == "SAMPLE", !is.na(sample_id)) %>%
    distinct(study_id, sample_id, patient_id)

  ambiguous_sample_map <- sample_map %>%
    count(study_id, sample_id, name = "n_patients") %>%
    filter(n_patients != 1L)

  if (nrow(ambiguous_sample_map) > 0) {
    stop("Some cBioPortal samples do not map uniquely to one patient.")
  }

  sample_values <- clinical_long %>%
    filter(clinical_level == "SAMPLE", !is.na(sample_id)) %>%
    select(
      study_id,
      sample_id,
      patient_id,
      clinical_attribute_id,
      clinical_value,
      clinical_level)

  patient_values <- sample_map %>%
    inner_join(
      clinical_long %>%
        filter(clinical_level == "PATIENT") %>%
        select(
          study_id,
          patient_id,
          clinical_attribute_id,
          clinical_value,
          clinical_level),
      by = c("study_id", "patient_id"),
      relationship = "many-to-many") %>%
    select(
      study_id,
      sample_id,
      patient_id,
      clinical_attribute_id,
      clinical_value,
      clinical_level)

  sample_patient_values <- bind_rows(sample_values, patient_values) %>%
    group_by(
      study_id,
      sample_id,
    patient_id,
      clinical_attribute_id) %>%
    summarise(
      n_distinct_values = n_distinct(na.omit(clinical_value)),
      clinical_value = first_non_empty(clinical_value),
      clinical_level = collapse_unique(clinical_level, n = Inf, sep = ";"),
      .groups = "drop")

  conflicting_values <- sample_patient_values %>%
    filter(n_distinct_values > 1L)

  if (nrow(conflicting_values) > 0) {
    stop(
      "Conflicting sample- and patient-level values were found for cBioPortal ",
      "clinical attributes.")
  }

  sample_patient_values %>%
    select(
      study_id,
      sample_id,
      patient_id,
      clinical_attribute_id,
      clinical_value) %>%
    pivot_wider(
      names_from = clinical_attribute_id,
      values_from = clinical_value)
}

normalise_cbio_t_stage <- function(x) {
  x <- str_to_upper(trimws(as.character(x)))
  x <- str_remove(x, "^[PR]")

  case_when(
    str_detect(x, "^T[12]") ~ "T1/T2",
    str_detect(x, "^T[34]") ~ "T3/T4",
    TRUE ~ NA_character_)
}

normalise_cbio_n_stage <- function(x, low_label = "N0", high_label = "N1") {
  x <- str_to_upper(trimws(as.character(x)))
  x <- str_remove(x, "^R")
  stage_number <- str_extract(x, "[0-4]")

  case_when(
    stage_number == "0" ~ low_label,
    stage_number == "1" ~ high_label,
    TRUE ~ NA_character_)
}

normalise_cbio_m_stage <- function(x) {
  x <- str_to_upper(trimws(as.character(x)))
  x <- str_remove(x, "^R")
  stage_number <- str_extract(x, "[01]")

  case_when(
    stage_number == "0" ~ "M0",
    stage_number == "1" ~ "M1",
    TRUE ~ NA_character_)
}

normalise_cbio_ranson_nodal_stage <- function(x) {
  stage_number <- str_extract(str_to_upper(trimws(as.character(x))), "[1-3]")

  case_when(
    stage_number %in% c("1", "2") ~ "N1/N2",
    stage_number == "3" ~ "N3",
    TRUE ~ NA_character_)
}

normalise_cbio_grade <- function(x) {
  grade_number <- suppressWarnings(parse_number(as.character(x)))

  case_when(
    grade_number %in% c(1, 2) ~ "Grade 1/2",
    grade_number == 3 ~ "Grade 3",
    TRUE ~ NA_character_)
}

make_cbio_endpoint_rows <- function(
    sample_tbl,
    endpoint_id,
    endpoint_raw_value,
    endpoint_group,
    endpoint_value_source = NULL) {

  endpoint_definition <- CBIO_CLINICAL_ENDPOINT_DEFINITIONS %>%
    filter(.data$endpoint_id == .env$endpoint_id)

  if (nrow(endpoint_definition) != 1L) {
    stop("Could not uniquely resolve clinical endpoint ", endpoint_id, ".")
  }

  if (is.null(endpoint_value_source)) {
    endpoint_value_source <- rep(
      paste0(
        "cBioPortal clinical attribute: ",
        endpoint_definition$clinical_attribute[[1]]),
      nrow(sample_tbl))
  }

  sample_tbl %>%
    transmute(
      endpoint_id = endpoint_id,
      analysis_family = endpoint_definition$analysis_family[[1]],
      analysis_scope = endpoint_definition$analysis_scope[[1]],
      endpoint_label = endpoint_definition$endpoint_label[[1]],
      endpoint_order = endpoint_definition$endpoint_order[[1]],
      clinical_attribute = endpoint_definition$clinical_attribute[[1]],
      reference_group = endpoint_definition$reference_group[[1]],
      case_group = endpoint_definition$case_group[[1]],
      study_id,
      sample_id,
      patient_id,
      endpoint_raw_value = as.character(endpoint_raw_value),
      endpoint_value_source = as.character(endpoint_value_source),
      endpoint_group = as.character(endpoint_group),
      endpoint_role = case_when(
        endpoint_group == endpoint_definition$reference_group[[1]] ~ "Reference",
        endpoint_group == endpoint_definition$case_group[[1]] ~ "Higher-risk",
        TRUE ~ NA_character_),
      ajcc_staging_edition = as.character(pull_existing_column(
        sample_tbl,
        c("AJCC_STAGING_EDITION"))),
      tmb_nonsynonymous = suppressWarnings(parse_number(as.character(
        pull_existing_column(sample_tbl, c("TMB_NONSYNONYMOUS"))))))
}

build_cbio_clinical_endpoint_rows <- function(sample_clinical) {
  hgsc <- sample_clinical %>%
    filter(study_id == "cscc_hgsc_bcm_2014")
  ranson <- sample_clinical %>%
    filter(study_id == "cscc_ranson_2022")
  dfci <- sample_clinical %>%
    filter(study_id == "cscc_dfarber_2015")
  ucsf <- sample_clinical %>%
    filter(study_id == "cscc_ucsf_2021")

  sample_type_pooled <- bind_rows(dfci, hgsc, ranson, ucsf) %>%
    mutate(
      analysis_sample_type = case_when(
        study_id == "cscc_dfarber_2015" ~ "Metastasis",
        study_id == "cscc_ucsf_2021" ~ "Primary",
        SAMPLE_TYPE == "Primary" ~ "Primary",
        SAMPLE_TYPE == "Metastasis" ~ "Metastasis",
        TRUE ~ NA_character_),
      analysis_sample_type_raw = case_when(
        study_id == "cscc_dfarber_2015" ~
          "Study description: metastatic cSCC cohort",
        study_id == "cscc_ucsf_2021" ~
          "Study description: primary cSCC tumour cohort",
        TRUE ~ as.character(SAMPLE_TYPE)),
      analysis_sample_type_source = case_when(
        study_id == "cscc_dfarber_2015" ~
          "DFCI study design (metastatic cSCC)",
        study_id == "cscc_ucsf_2021" ~
          "UCSF study design (primary cSCC)",
        TRUE ~ "cBioPortal clinical attribute: SAMPLE_TYPE"))

  sample_type_explicit <- bind_rows(hgsc, ranson)
  recurrence_pooled <- bind_rows(
    dfci %>%
      mutate(
        recurrence_raw_value = as.character(DFS_STATUS),
        recurrence_group = case_when(
          str_detect(DFS_STATUS, "^0:") ~ "Disease free",
          str_detect(DFS_STATUS, "^1:") ~ "Recurred / progressed",
          TRUE ~ NA_character_),
        recurrence_value_source =
          "cBioPortal clinical attribute: DFS_STATUS"),
    hgsc %>%
      mutate(
        recurrence_raw_value =
          as.character(DISEASE_RECURRENCE_OR_PERSISTENCE),
        recurrence_group = case_when(
          DISEASE_RECURRENCE_OR_PERSISTENCE == "No" ~ "Disease free",
          DISEASE_RECURRENCE_OR_PERSISTENCE == "Yes" ~
            "Recurred / progressed",
          TRUE ~ NA_character_),
        recurrence_value_source = paste(
          "cBioPortal clinical attribute:",
          "DISEASE_RECURRENCE_OR_PERSISTENCE")))

  bind_rows(
    make_cbio_endpoint_rows(
      sample_tbl = sample_type_pooled,
      endpoint_id = "sample_type_pooled",
      endpoint_raw_value = sample_type_pooled$analysis_sample_type_raw,
      endpoint_group = sample_type_pooled$analysis_sample_type,
      endpoint_value_source =
        sample_type_pooled$analysis_sample_type_source),
    make_cbio_endpoint_rows(
      sample_tbl = sample_type_explicit,
      endpoint_id = "sample_type_explicit",
      endpoint_raw_value = sample_type_explicit$SAMPLE_TYPE,
      endpoint_group = case_when(
        sample_type_explicit$SAMPLE_TYPE == "Primary" ~ "Primary",
        sample_type_explicit$SAMPLE_TYPE == "Metastasis" ~ "Metastasis",
        TRUE ~ NA_character_)),
    make_cbio_endpoint_rows(
      sample_tbl = hgsc,
      endpoint_id = "sample_type_hgsc",
      endpoint_raw_value = hgsc$SAMPLE_TYPE,
      endpoint_group = case_when(
        hgsc$SAMPLE_TYPE == "Primary" ~ "Primary",
        hgsc$SAMPLE_TYPE == "Metastasis" ~ "Metastasis",
        TRUE ~ NA_character_)),
    make_cbio_endpoint_rows(
      sample_tbl = hgsc,
      endpoint_id = "hgsc_t_stage",
      endpoint_raw_value = hgsc$T_STAGE,
      endpoint_group = normalise_cbio_t_stage(hgsc$T_STAGE)),
    make_cbio_endpoint_rows(
      sample_tbl = hgsc,
      endpoint_id = "hgsc_n_stage",
      endpoint_raw_value = hgsc$N_STAGE,
      endpoint_group = normalise_cbio_n_stage(hgsc$N_STAGE)),
    make_cbio_endpoint_rows(
      sample_tbl = hgsc,
      endpoint_id = "hgsc_m_stage",
      endpoint_raw_value = hgsc$M_STAGE,
      endpoint_group = normalise_cbio_m_stage(hgsc$M_STAGE)),
    make_cbio_endpoint_rows(
      sample_tbl = hgsc,
      endpoint_id = "hgsc_recurrence_persistence",
      endpoint_raw_value = hgsc$DISEASE_RECURRENCE_OR_PERSISTENCE,
      endpoint_group = case_when(
        hgsc$DISEASE_RECURRENCE_OR_PERSISTENCE == "No" ~ "No",
        hgsc$DISEASE_RECURRENCE_OR_PERSISTENCE == "Yes" ~ "Yes",
        TRUE ~ NA_character_)),
    make_cbio_endpoint_rows(
      sample_tbl = hgsc,
      endpoint_id = "hgsc_efs_status",
      endpoint_raw_value = hgsc$EFS_STATUS,
      endpoint_group = case_when(
        str_detect(hgsc$EFS_STATUS, "^0:") ~ "Censored",
        str_detect(hgsc$EFS_STATUS, "^1:") ~ "Event",
        TRUE ~ NA_character_)),
    make_cbio_endpoint_rows(
      sample_tbl = ranson,
      endpoint_id = "ranson_nodal_stage",
      endpoint_raw_value = ranson$NODAL_STAGE,
      endpoint_group = normalise_cbio_ranson_nodal_stage(ranson$NODAL_STAGE)),
    make_cbio_endpoint_rows(
      sample_tbl = ranson,
      endpoint_id = "ranson_grade",
      endpoint_raw_value = ranson$GRADE,
      endpoint_group = normalise_cbio_grade(ranson$GRADE)),
    make_cbio_endpoint_rows(
      sample_tbl = ranson,
      endpoint_id = "ranson_extracapsular_spread",
      endpoint_raw_value = ranson$EXTRACAPSULAR_SPREAD,
      endpoint_group = case_when(
        ranson$EXTRACAPSULAR_SPREAD == "No" ~ "No",
        ranson$EXTRACAPSULAR_SPREAD == "Yes" ~ "Yes",
        TRUE ~ NA_character_)),
    make_cbio_endpoint_rows(
      sample_tbl = dfci,
      endpoint_id = "dfci_dfs_status",
      endpoint_raw_value = dfci$DFS_STATUS,
      endpoint_group = case_when(
        str_detect(dfci$DFS_STATUS, "^0:") ~ "Disease free",
        str_detect(dfci$DFS_STATUS, "^1:") ~ "Recurred / progressed",
        TRUE ~ NA_character_)),
    make_cbio_endpoint_rows(
      sample_tbl = recurrence_pooled,
      endpoint_id = "recurrence_pooled",
      endpoint_raw_value = recurrence_pooled$recurrence_raw_value,
      endpoint_group = recurrence_pooled$recurrence_group,
      endpoint_value_source =
        recurrence_pooled$recurrence_value_source)) %>%
    arrange(endpoint_order, study_id, sample_id)
}

summarise_cbio_endpoint_availability <- function(endpoint_rows) {
  summarise_availability <- function(tbl, availability_level) {
    if (availability_level == "Endpoint") {
      tbl <- tbl %>%
        mutate(study_id = "All included studies")
    }

    tbl %>%
      group_by(
        endpoint_id,
        analysis_family,
        analysis_scope,
        endpoint_label,
        endpoint_order,
        clinical_attribute,
        reference_group,
        case_group,
        study_id) %>%
      summarise(
        n_samples_total = n_distinct(sample_id),
        n_samples_analysed = n_distinct(sample_id[!is.na(endpoint_group)]),
        n_reference = n_distinct(
          sample_id[endpoint_group == first(reference_group)],
          na.rm = TRUE),
        n_case = n_distinct(
          sample_id[endpoint_group == first(case_group)],
          na.rm = TRUE),
        n_excluded = n_samples_total - n_samples_analysed,
        excluded_raw_values = paste(
          sort(unique(na.omit(endpoint_raw_value[is.na(endpoint_group)]))),
          collapse = ";"),
        n_with_ajcc_staging_edition = n_distinct(
          sample_id[!is.na(ajcc_staging_edition)]),
        .groups = "drop") %>%
      mutate(availability_level = availability_level)
  }

  bind_rows(
    summarise_availability(endpoint_rows, "Endpoint"),
    summarise_availability(endpoint_rows, "Study")) %>%
    arrange(endpoint_order, availability_level, study_id)
}

prepare_cbio_candidate_clinical_matrix <- function(
    oncoplot_df,
    candidate_tbl) {

  candidate_metadata <- candidate_tbl %>%
    transmute(
      gene_symbol,
      gene_label = as.character(gene_label),
      candidate_group = as.character(candidate_group),
      human_symbols) %>%
    distinct()

  oncoplot_df %>%
    semi_join(candidate_metadata, by = "gene_symbol") %>%
    select(
      study_id,
      sample_id,
      cbio_sample_uid,
      gene_symbol,
      is_gene_profiled,
      n_cell_mutations) %>%
    left_join(candidate_metadata, by = "gene_symbol") %>%
    mutate(
      is_gene_profiled = replace_na(is_gene_profiled, FALSE),
      is_mutated = is_gene_profiled & replace_na(n_cell_mutations, 0L) > 0L,
      mutation_definition = paste(
        "Any nonsynonymous mutation in the mapped human orthologue;",
        "exact mouse variants are not matched.")) %>%
    distinct(
      study_id,
      sample_id,
      cbio_sample_uid,
      gene_symbol,
      .keep_all = TRUE)
}

p_adjust_with_na <- function(p_values) {
  adjusted <- rep(NA_real_, length(p_values))
  keep <- !is.na(p_values)

  if (any(keep)) {
    adjusted[keep] <- p.adjust(p_values[keep], method = "BH")
  }

  adjusted
}

summarise_cbio_gene_endpoint_test <- function(tbl) {
  reference_group <- first(tbl$reference_group)
  case_group <- first(tbl$case_group)
  test_tbl <- tbl %>%
    filter(
      is_gene_profiled,
      endpoint_group %in% c(reference_group, case_group))

  n_reference_profiled <- sum(test_tbl$endpoint_group == reference_group)
  n_case_profiled <- sum(test_tbl$endpoint_group == case_group)
  n_reference_mutated <- sum(
    test_tbl$endpoint_group == reference_group & test_tbl$is_mutated)
  n_case_mutated <- sum(
    test_tbl$endpoint_group == case_group & test_tbl$is_mutated)
  n_total_profiled <- n_reference_profiled + n_case_profiled
  n_total_mutated <- n_reference_mutated + n_case_mutated

  test_status <- case_when(
    n_reference_profiled == 0L | n_case_profiled == 0L ~
      "not_testable_missing_endpoint_group",
    n_total_mutated == 0L | n_total_mutated == n_total_profiled ~
      "not_testable_no_mutation_variation",
    TRUE ~ "tested")

  contingency_table <- matrix(
    c(
      n_case_mutated,
      n_case_profiled - n_case_mutated,
      n_reference_mutated,
      n_reference_profiled - n_reference_mutated),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      endpoint = c(case_group, reference_group),
      mutation = c("Mutated", "Not mutated")))

  fisher_result <- if (test_status == "tested") {
    fisher.test(contingency_table, alternative = "two.sided")
  } else {
    NULL
  }

  corrected_odds_ratio <- (
    (n_case_mutated + 0.5) *
      (n_reference_profiled - n_reference_mutated + 0.5)) /
    (
      (n_case_profiled - n_case_mutated + 0.5) *
        (n_reference_mutated + 0.5))

  tibble(
    analysis_family = first(tbl$analysis_family),
    analysis_scope = first(tbl$analysis_scope),
    endpoint_label = first(tbl$endpoint_label),
    endpoint_order = first(tbl$endpoint_order),
    clinical_attribute = first(tbl$clinical_attribute),
    reference_group = reference_group,
    case_group = case_group,
    gene_label = first(tbl$gene_label),
    candidate_group = first(tbl$candidate_group),
    human_symbols = first(tbl$human_symbols),
    mutation_definition = first(tbl$mutation_definition),
    test_status = test_status,
    n_reference_profiled = n_reference_profiled,
    n_reference_mutated = n_reference_mutated,
    reference_prevalence_pct = if_else(
      n_reference_profiled > 0,
      100 * n_reference_mutated / n_reference_profiled,
      NA_real_),
    n_case_profiled = n_case_profiled,
    n_case_mutated = n_case_mutated,
    case_prevalence_pct = if_else(
      n_case_profiled > 0,
      100 * n_case_mutated / n_case_profiled,
      NA_real_),
    prevalence_difference_pct = case_prevalence_pct - reference_prevalence_pct,
    odds_ratio = if (is.null(fisher_result)) {
      NA_real_
    } else {
      unname(fisher_result$estimate)
    },
    odds_ratio_conf_low = if (is.null(fisher_result)) {
      NA_real_
    } else {
      unname(fisher_result$conf.int[[1]])
    },
    odds_ratio_conf_high = if (is.null(fisher_result)) {
      NA_real_
    } else {
      unname(fisher_result$conf.int[[2]])
    },
    corrected_odds_ratio = corrected_odds_ratio,
    corrected_log2_odds_ratio = log2(corrected_odds_ratio),
    p_value = if (is.null(fisher_result)) {
      NA_real_
    } else {
      fisher_result$p.value
    })
}

run_cbio_candidate_clinical_associations <- function(
    mutation_matrix,
    endpoint_rows) {

  endpoint_mutation_tbl <- mutation_matrix %>%
    inner_join(
      endpoint_rows %>%
        filter(!is.na(endpoint_group)),
      by = c("study_id", "sample_id"),
      relationship = "many-to-many")

  endpoint_mutation_tbl %>%
    group_by(endpoint_id, gene_symbol) %>%
    group_modify(~ summarise_cbio_gene_endpoint_test(.x)) %>%
    ungroup() %>%
    group_by(endpoint_id) %>%
    mutate(q_value_within_endpoint = p_adjust_with_na(p_value)) %>%
    ungroup() %>%
    group_by(analysis_family) %>%
    mutate(q_value_analysis_family = p_adjust_with_na(p_value)) %>%
    ungroup() %>%
    arrange(endpoint_order, candidate_group, gene_label)
}

build_cbio_candidate_human_metadata <- function(candidate_tbl) {
  candidate_tbl %>%
    transmute(
      candidate_gene_symbol = gene_symbol,
      candidate_gene_label = as.character(gene_label),
      candidate_group = as.character(candidate_group),
      human_symbols = coalesce(human_symbols, "")) %>%
    separate_rows(human_symbols, sep = ";") %>%
    transmute(
      human_symbol = trimws(human_symbols),
      candidate_gene_symbol,
      candidate_gene_label,
      candidate_group) %>%
    filter(human_symbol != "") %>%
    group_by(human_symbol) %>%
    summarise(
      is_candidate = TRUE,
      candidate_gene_symbols = collapse_unique(
        candidate_gene_symbol,
        n = Inf,
        sep = ";"),
      candidate_gene_labels = collapse_unique(
        candidate_gene_label,
        n = Inf,
        sep = ";"),
      candidate_groups = collapse_unique(
        candidate_group,
        n = Inf,
        sep = ";"),
      candidate_group = first(candidate_group),
      .groups = "drop")
}

run_cbio_all_gene_primary_metastatic_associations <- function(
    incidence_tbl,
    profiled_map,
    cbio_mutations,
    endpoint_rows,
    candidate_tbl) {

  sample_type_rows <- endpoint_rows %>%
    filter(
      endpoint_id == "sample_type_pooled",
      !is.na(endpoint_group)) %>%
    transmute(
      study_id,
      sample_id,
      cbio_sample_uid = paste(study_id, sample_id, sep = "::"),
      endpoint_group) %>%
    distinct()

  profiled_counts <- profiled_map %>%
    transmute(
      human_symbol,
      study_id,
      sample_id,
      cbio_sample_uid) %>%
    distinct() %>%
    inner_join(
      sample_type_rows,
      by = c("study_id", "sample_id", "cbio_sample_uid")) %>%
    group_by(human_symbol) %>%
    summarise(
      n_reference_profiled = n_distinct(
        cbio_sample_uid[endpoint_group == "Primary"]),
      n_case_profiled = n_distinct(
        cbio_sample_uid[endpoint_group == "Metastasis"]),
      .groups = "drop")

  mutated_counts <- cbio_mutations %>%
    filter(
      !CBIO_RESTRICT_TO_NONSYNONYMOUS |
        mutation_class %in% CBIO_NONSYNONYMOUS_CLASSES) %>%
    transmute(
      human_symbol = hugo_symbol,
      study_id,
      sample_id,
      cbio_sample_uid) %>%
    distinct() %>%
    semi_join(
      profiled_map %>%
        select(
          human_symbol,
          study_id,
          sample_id,
          cbio_sample_uid) %>%
        distinct(),
      by = c(
        "human_symbol",
        "study_id",
        "sample_id",
        "cbio_sample_uid")) %>%
    inner_join(
      sample_type_rows,
      by = c("study_id", "sample_id", "cbio_sample_uid")) %>%
    group_by(human_symbol) %>%
    summarise(
      n_reference_mutated = n_distinct(
        cbio_sample_uid[endpoint_group == "Primary"]),
      n_case_mutated = n_distinct(
        cbio_sample_uid[endpoint_group == "Metastasis"]),
      .groups = "drop")

  candidate_metadata <- build_cbio_candidate_human_metadata(candidate_tbl)

  association_tbl <- incidence_tbl %>%
    transmute(
      human_symbol,
      incidence_n_mutations = n_mutations,
      incidence_n_altered_samples = n_altered_samples,
      incidence_n_profiled_samples = n_profiled_samples,
      incidence_pct,
      is_cancer_gene) %>%
    left_join(profiled_counts, by = "human_symbol") %>%
    left_join(mutated_counts, by = "human_symbol") %>%
    left_join(candidate_metadata, by = "human_symbol") %>%
    mutate(
      across(
        c(
          n_reference_profiled,
          n_case_profiled,
          n_reference_mutated,
          n_case_mutated),
        ~ replace_na(as.integer(.x), 0L)),
      is_candidate = replace_na(is_candidate, FALSE),
      n_total_profiled = n_reference_profiled + n_case_profiled,
      n_total_mutated = n_reference_mutated + n_case_mutated,
      profile_denominator_matches_incidence =
        n_total_profiled == incidence_n_profiled_samples,
      test_status = case_when(
        n_reference_profiled == 0L | n_case_profiled == 0L ~
          "not_testable_missing_endpoint_group",
        n_total_mutated == 0L | n_total_mutated == n_total_profiled ~
          "not_testable_no_mutation_variation",
        TRUE ~ "tested"),
      reference_prevalence_pct = if_else(
        n_reference_profiled > 0L,
        100 * n_reference_mutated / n_reference_profiled,
        NA_real_),
      case_prevalence_pct = if_else(
        n_case_profiled > 0L,
        100 * n_case_mutated / n_case_profiled,
        NA_real_),
      prevalence_difference_pct =
        case_prevalence_pct - reference_prevalence_pct,
      corrected_odds_ratio = (
        (n_case_mutated + 0.5) *
          (n_reference_profiled - n_reference_mutated + 0.5)) /
        (
          (n_case_profiled - n_case_mutated + 0.5) *
            (n_reference_mutated + 0.5)),
      corrected_log2_odds_ratio = log2(corrected_odds_ratio))

  tested_indices <- which(association_tbl$test_status == "tested")
  fisher_results <- matrix(
    NA_real_,
    nrow = nrow(association_tbl),
    ncol = 2,
    dimnames = list(NULL, c("odds_ratio", "p_value")))

  fisher_results[tested_indices, ] <- t(vapply(
    tested_indices,
    function(row_index) {
      row_data <- association_tbl[row_index, ]
      contingency_table <- matrix(
        c(
          row_data$n_case_mutated,
          row_data$n_case_profiled - row_data$n_case_mutated,
          row_data$n_reference_mutated,
          row_data$n_reference_profiled - row_data$n_reference_mutated),
        nrow = 2,
        byrow = TRUE)
      fisher_result <- fisher.test(
        contingency_table,
        alternative = "two.sided",
        conf.int = FALSE)

      c(
        odds_ratio = unname(fisher_result$estimate),
        p_value = fisher_result$p.value)
    },
    numeric(2)))

  association_tbl %>%
    mutate(
      odds_ratio = fisher_results[, "odds_ratio"],
      p_value = fisher_results[, "p_value"],
      q_value_all_mutated_genes = p_adjust_with_na(p_value),
      n_fdr_tests = sum(!is.na(p_value)),
      association_direction = case_when(
        test_status != "tested" ~ "Not testable",
        q_value_all_mutated_genes <= 0.05 &
          corrected_log2_odds_ratio < 0 ~ "Primary-associated",
        q_value_all_mutated_genes <= 0.05 &
          corrected_log2_odds_ratio > 0 ~ "Metastatic-associated",
        TRUE ~ "Not significant"),
      association_direction = factor(
        association_direction,
        levels = c(
          "Primary-associated",
          "Metastatic-associated",
          "Not significant",
          "Not testable")),
      is_significant_candidate =
        is_candidate &
        association_direction %in% c(
          "Primary-associated",
          "Metastatic-associated"),
      volcano_category = case_when(
        is_significant_candidate &
          association_direction == "Primary-associated" ~
          "Candidate: primary-associated",
        is_significant_candidate &
          association_direction == "Metastatic-associated" ~
          "Candidate: metastatic-associated",
        association_direction %in% c(
          "Primary-associated",
          "Metastatic-associated") ~
          "Other FDR-significant gene",
        test_status == "tested" ~ "Not significant",
        TRUE ~ "Not testable"),
      candidate_group_plot = factor(
        candidate_group_display(candidate_group),
        levels = CBIO_CANDIDATE_GROUP_LEVELS)) %>%
    arrange(p_value, human_symbol)
}

summarise_cbio_primary_metastatic_prevalence <- function(
    mutation_matrix,
    endpoint_rows) {

  sample_type_rows <- endpoint_rows %>%
    filter(
      endpoint_id == "sample_type_pooled",
      !is.na(endpoint_group)) %>%
    select(study_id, sample_id, endpoint_group)

  prevalence_by_study <- mutation_matrix %>%
    inner_join(
      sample_type_rows,
      by = c("study_id", "sample_id")) %>%
    filter(is_gene_profiled) %>%
    group_by(
      study_id,
      endpoint_group,
      gene_symbol,
      gene_label,
      candidate_group,
      human_symbols) %>%
    summarise(
      n_profiled = n(),
      n_mutated = sum(is_mutated),
      prevalence_pct = 100 * n_mutated / n_profiled,
      .groups = "drop") %>%
    mutate(
      comparison_cohort = recode(
        study_id,
        cscc_dfarber_2015 = "DFCI",
        cscc_hgsc_bcm_2014 = "HGSC",
        cscc_ranson_2022 = "Ranson",
        cscc_ucsf_2021 = "UCSF"))

  pooled_prevalence <- mutation_matrix %>%
    inner_join(
      sample_type_rows,
      by = c("study_id", "sample_id")) %>%
    filter(is_gene_profiled) %>%
    group_by(
      endpoint_group,
      gene_symbol,
      gene_label,
      candidate_group,
      human_symbols) %>%
    summarise(
      n_profiled = n(),
      n_mutated = sum(is_mutated),
      prevalence_pct = 100 * n_mutated / n_profiled,
      study_id = paste(
        sort(unique(study_id)),
        collapse = ";"),
      .groups = "drop") %>%
    mutate(comparison_cohort = "Pooled")

  bind_rows(prevalence_by_study, pooled_prevalence) %>%
    mutate(
      endpoint_group = factor(
        endpoint_group,
        levels = c("Primary", "Metastasis")),
      comparison_cohort = factor(
        comparison_cohort,
        levels = c("DFCI", "HGSC", "Ranson", "UCSF", "Pooled"))) %>%
    arrange(comparison_cohort, candidate_group, gene_label, endpoint_group)
}

summarise_cbio_candidate_recurrence_cmh <- function(tbl) {
  test_tbl <- tbl %>%
    filter(is_gene_profiled) %>%
    mutate(
      endpoint_group = factor(
        endpoint_group,
        levels = c("Recurred / progressed", "Disease free")),
      mutation_status = factor(
        if_else(is_mutated, "Mutated", "Not mutated"),
        levels = c("Mutated", "Not mutated")),
      study_id = factor(study_id))

  study_counts <- test_tbl %>%
    count(
      study_id,
      endpoint_group,
      mutation_status,
      .drop = FALSE) %>%
    pivot_wider(
      names_from = c(endpoint_group, mutation_status),
      values_from = n,
      values_fill = 0,
      names_sep = "__")

  case_mutated <- study_counts$`Recurred / progressed__Mutated`
  case_not_mutated <- study_counts$`Recurred / progressed__Not mutated`
  reference_mutated <- study_counts$`Disease free__Mutated`
  reference_not_mutated <- study_counts$`Disease free__Not mutated`
  stratum_totals <- case_mutated + case_not_mutated +
    reference_mutated + reference_not_mutated
  informative_strata <- (
    case_mutated + case_not_mutated > 0L) &
    (reference_mutated + reference_not_mutated > 0L) &
    (case_mutated + reference_mutated > 0L) &
    (case_not_mutated + reference_not_mutated > 0L)

  corrected_common_odds_ratio <- sum(
    (case_mutated + 0.5) *
      (reference_not_mutated + 0.5) /
      (stratum_totals + 2)) /
    sum(
      (case_not_mutated + 0.5) *
        (reference_mutated + 0.5) /
        (stratum_totals + 2))

  contingency_array <- xtabs(
    ~ endpoint_group + mutation_status + study_id,
    data = test_tbl,
    drop.unused.levels = FALSE)
  test_status <- case_when(
    sum(test_tbl$endpoint_group == "Disease free") == 0L |
      sum(test_tbl$endpoint_group == "Recurred / progressed") == 0L ~
        "not_testable_missing_endpoint_group",
    sum(test_tbl$is_mutated) == 0L |
      sum(test_tbl$is_mutated) == nrow(test_tbl) ~
        "not_testable_no_mutation_variation",
    !any(informative_strata) ~ "not_testable_no_informative_study",
    TRUE ~ "tested")
  association_result <- if (
      test_status == "tested" &&
      sum(informative_strata) >= 2L) {
    tryCatch(
      mantelhaen.test(
        contingency_array[, , informative_strata, drop = FALSE],
        correct = FALSE,
        conf.level = 0.95),
      error = function(e) NULL)
  } else if (
      test_status == "tested" &&
      sum(informative_strata) == 1L) {
    informative_index <- which(informative_strata)[[1]]
    single_study_table <- matrix(
      c(
        case_mutated[[informative_index]],
        case_not_mutated[[informative_index]],
        reference_mutated[[informative_index]],
        reference_not_mutated[[informative_index]]),
      nrow = 2,
      byrow = TRUE)
    fisher.test(single_study_table, alternative = "two.sided")
  } else {
    NULL
  }
  if (test_status == "tested" && is.null(association_result)) {
    test_status <- "not_testable_association_failure"
  }
  test_method <- if (test_status != "tested") {
    NA_character_
  } else if (sum(informative_strata) >= 2L) {
    paste(
      "Cochran-Mantel-Haenszel test stratified by study;",
      "continuity correction disabled")
  } else {
    paste0(
      "Two-sided Fisher exact test within ",
      recode(
        as.character(study_counts$study_id[informative_strata][[1]]),
        cscc_dfarber_2015 = "DFCI",
        cscc_hgsc_bcm_2014 = "HGSC"),
      "; only one study profiles this candidate")
  }
  study_adjustment <- if (test_status != "tested") {
    NA_character_
  } else if (sum(informative_strata) >= 2L) {
    "Study-stratified CMH"
  } else {
    "Not applicable: candidate profiled in one informative study"
  }

  tibble(
    analysis_family = "Recurrence / progression",
    analysis_scope = "Study-adjusted primary analysis",
    endpoint_id = "recurrence_pooled",
    endpoint_label =
      "Disease-free vs recurrent / progressed (DFCI + HGSC)",
    reference_group = "Disease free",
    case_group = "Recurred / progressed",
    gene_label = first(tbl$gene_label),
    candidate_group = first(tbl$candidate_group),
    human_symbols = first(tbl$human_symbols),
    mutation_definition = first(tbl$mutation_definition),
    n_reference_profiled =
      sum(test_tbl$endpoint_group == "Disease free"),
    n_reference_mutated = sum(
      test_tbl$endpoint_group == "Disease free" & test_tbl$is_mutated),
    n_case_profiled =
      sum(test_tbl$endpoint_group == "Recurred / progressed"),
    n_case_mutated = sum(
      test_tbl$endpoint_group == "Recurred / progressed" &
        test_tbl$is_mutated),
    n_informative_studies = sum(informative_strata),
    informative_studies = paste(
      as.character(study_counts$study_id[informative_strata]),
      collapse = ";"),
    test_method = test_method,
    study_adjustment = study_adjustment,
    test_status = test_status,
    common_odds_ratio = if (is.null(association_result)) {
      NA_real_
    } else {
      unname(association_result$estimate)
    },
    common_odds_ratio_conf_low = if (is.null(association_result)) {
      NA_real_
    } else {
      unname(association_result$conf.int[[1]])
    },
    common_odds_ratio_conf_high = if (is.null(association_result)) {
      NA_real_
    } else {
      unname(association_result$conf.int[[2]])
    },
    corrected_common_odds_ratio = corrected_common_odds_ratio,
    corrected_log2_common_odds_ratio =
      log2(corrected_common_odds_ratio),
    p_value = if (is.null(association_result)) {
      NA_real_
    } else {
      association_result$p.value
    })
}

run_cbio_candidate_recurrence_adjusted_associations <- function(
    mutation_matrix,
    endpoint_rows) {

  recurrence_rows <- endpoint_rows %>%
    filter(
      endpoint_id == "recurrence_pooled",
      !is.na(endpoint_group)) %>%
    select(study_id, sample_id, endpoint_group) %>%
    distinct()

  mutation_matrix %>%
    inner_join(
      recurrence_rows,
      by = c("study_id", "sample_id")) %>%
    group_by(gene_symbol) %>%
    group_modify(~ summarise_cbio_candidate_recurrence_cmh(.x)) %>%
    ungroup() %>%
    mutate(
      q_value = p_adjust_with_na(p_value),
      association_direction = case_when(
        test_status != "tested" ~ "Not testable",
        q_value <= 0.05 &
          corrected_log2_common_odds_ratio < 0 ~
            "Disease-free-associated",
        q_value <= 0.05 &
          corrected_log2_common_odds_ratio > 0 ~
            "Recurrence-associated",
        TRUE ~ "Not significant"),
      association_direction = factor(
        association_direction,
        levels = c(
          "Disease-free-associated",
          "Recurrence-associated",
          "Not significant",
          "Not testable")),
      candidate_group_plot = factor(
        candidate_group_display(candidate_group),
        levels = CBIO_CANDIDATE_GROUP_LEVELS)) %>%
    arrange(p_value, gene_label)
}

build_cbio_focused_sensitivity_endpoint_rows <- function(endpoint_rows) {
  bind_rows(lapply(
    seq_len(nrow(CBIO_FOCUSED_SENSITIVITY_DEFINITIONS)),
    function(definition_index) {
      endpoint_definition <- CBIO_FOCUSED_SENSITIVITY_DEFINITIONS[
        definition_index,
        ]

      focused_rows <- endpoint_rows %>%
        filter(
          endpoint_id == endpoint_definition$source_endpoint_id[[1]],
          !is.na(endpoint_group))

      excluded_study_id <- endpoint_definition$excluded_study_id[[1]]
      if (!is.na(excluded_study_id)) {
        focused_rows <- focused_rows %>%
          filter(study_id != excluded_study_id)
      }

      focused_rows %>%
        mutate(
          endpoint_id = endpoint_definition$endpoint_id[[1]],
          analysis_family = "Focused primary / metastatic sensitivity",
          analysis_scope = endpoint_definition$analysis_scope[[1]],
          endpoint_label = endpoint_definition$endpoint_label[[1]],
          endpoint_order = endpoint_definition$endpoint_order[[1]])
    })) %>%
    arrange(endpoint_order, study_id, sample_id)
}

run_cbio_focused_primary_metastatic_associations <- function(
    mutation_matrix,
    focused_endpoint_rows) {

  focused_results <- run_cbio_candidate_clinical_associations(
    mutation_matrix = mutation_matrix %>%
      filter(gene_symbol %in% CBIO_FOCUSED_PRIMARY_METASTATIC_GENES),
    endpoint_rows = focused_endpoint_rows)

  focused_results %>%
    left_join(
      CBIO_FOCUSED_SENSITIVITY_DEFINITIONS %>%
        select(
          endpoint_id,
          source_endpoint_id,
          excluded_study_id,
          included_study_ids),
      by = "endpoint_id") %>%
    arrange(endpoint_order, gene_symbol)
}

plot_cbio_focused_primary_metastatic_prevalence <- function(prevalence_tbl) {
  plot_df <- prevalence_tbl %>%
    mutate(
      gene_label = factor(
        gene_label,
        levels = unique(gene_label[match(
          CBIO_FOCUSED_PRIMARY_METASTATIC_GENES,
          gene_symbol)])),
      sample_count_label = if_else(
        n_profiled > 0,
        paste0(n_mutated, "/", n_profiled),
        "not profiled"),
      sample_count_label_y = if_else(
        n_profiled > 0,
        prevalence_pct,
        4))

  dodge <- position_dodge(width = 0.75)

  ggplot(
    plot_df,
    aes(
      x = comparison_cohort,
      y = prevalence_pct,
      fill = endpoint_group,
      group = endpoint_group)) +
    geom_point(
      data = plot_df %>% filter(n_profiled > 0),
      shape = 21,
      colour = "black",
      stroke = 0.3,
      size = 3,
      position = dodge) +
    geom_text(
      aes(y = sample_count_label_y, label = sample_count_label),
      colour = "black",
      size = 7 / .pt,
      vjust = if_else(plot_df$n_profiled > 0, -0.9, 0.5),
      position = dodge) +
    facet_wrap(vars(gene_label), nrow = 1) +
    scale_fill_manual(
      values = c(
        Primary = "#4C78A8",
        Metastasis = "#D62728"),
      drop = FALSE) +
    scale_x_discrete(
      limits = c("DFCI", "HGSC", "Ranson", "UCSF", "Pooled"),
      drop = FALSE) +
    scale_y_continuous(
      breaks = c(0, 50, 100),
      labels = function(x) paste0(x, "%"),
      limits = c(0, 100),
      expand = expansion(mult = c(0.03, 0.14))) +
    labs(
      x = "",
      y = "Samples with a nonsynonymous mutation",
      fill = "Sample type") +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(size = 8),
      legend.position = "bottom")
}

plot_cbio_focused_primary_metastatic_sensitivity <- function(focused_results) {
  plot_df <- focused_results %>%
    mutate(
      endpoint_label = factor(
        endpoint_label,
        levels = rev(CBIO_FOCUSED_SENSITIVITY_DEFINITIONS$endpoint_label)),
      gene_label = factor(
        gene_label,
        levels = unique(gene_label[match(
          CBIO_FOCUSED_PRIMARY_METASTATIC_GENES,
          gene_symbol)])),
      evidence_strength = if_else(
        test_status == "tested",
        pmin(-log10(pmax(p_value, .Machine$double.xmin)), 5),
        NA_real_))

  ggplot(
    plot_df,
    aes(
      x = corrected_log2_odds_ratio,
      y = endpoint_label)) +
    geom_vline(
      xintercept = 0,
      colour = "grey55",
      linetype = "dashed",
      linewidth = 0.3) +
    geom_point(
      data = plot_df %>% filter(test_status == "tested"),
      aes(
        fill = corrected_log2_odds_ratio,
        size = evidence_strength),
      shape = 21,
      colour = "black",
      stroke = 0.3) +
    geom_text(
      data = plot_df %>% filter(test_status != "tested"),
      aes(x = 0, label = "not testable"),
      colour = "grey45",
      size = 7 / .pt) +
    facet_wrap(vars(gene_label), nrow = 1) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      oob = scales::squish,
      name = expression(log[2]("odds ratio"))) +
    scale_size_continuous(
      range = c(1.8, 4.5),
      limits = c(0, 5),
      name = expression(-log[10](italic(P)))) +
    labs(
      x = expression(log[2]("corrected odds ratio")),
      y = "") +
    my_theme +
    theme(
      strip.text = element_text(size = 8),
      legend.position = "right")
}

build_cbio_candidate_burden <- function(
    mutation_matrix,
    sample_clinical) {

  sample_tmb <- sample_clinical %>%
    transmute(
      study_id,
      sample_id,
      patient_id,
      tmb_nonsynonymous = suppressWarnings(parse_number(as.character(
        pull_existing_column(., c("TMB_NONSYNONYMOUS"))))))

  mutation_matrix %>%
    group_by(study_id, sample_id, cbio_sample_uid) %>%
    summarise(
      n_candidate_genes = n_distinct(gene_symbol),
      n_profiled_candidate_genes = sum(is_gene_profiled),
      n_mutated_candidate_genes = sum(is_mutated & is_gene_profiled),
      n_candidate_nonsynonymous_mutations = sum(
        if_else(is_gene_profiled, n_cell_mutations, 0L)),
      candidate_burden_fraction = if_else(
        n_profiled_candidate_genes > 0,
        n_mutated_candidate_genes / n_profiled_candidate_genes,
        NA_real_),
      mutated_candidate_genes = paste(
        sort(unique(gene_symbol[is_mutated & is_gene_profiled])),
        collapse = ";"),
      .groups = "drop") %>%
    left_join(
      sample_tmb,
      by = c("study_id", "sample_id"))
}

summarise_cbio_burden_endpoint_test <- function(tbl) {
  reference_group <- first(tbl$reference_group)
  case_group <- first(tbl$case_group)
  test_tbl <- tbl %>%
    filter(
      !is.na(candidate_burden_fraction),
      endpoint_group %in% c(reference_group, case_group))

  reference_values <- test_tbl$candidate_burden_fraction[
    test_tbl$endpoint_group == reference_group]
  case_values <- test_tbl$candidate_burden_fraction[
    test_tbl$endpoint_group == case_group]

  test_status <- case_when(
    length(reference_values) == 0L | length(case_values) == 0L ~
      "not_testable_missing_endpoint_group",
    n_distinct(c(reference_values, case_values)) < 2L ~
      "not_testable_no_burden_variation",
    TRUE ~ "tested")

  wilcoxon_result <- if (test_status == "tested") {
    suppressWarnings(wilcox.test(
      case_values,
      reference_values,
      alternative = "two.sided",
      exact = FALSE))
  } else {
    NULL
  }

  tibble(
    analysis_family = first(tbl$analysis_family),
    analysis_scope = first(tbl$analysis_scope),
    endpoint_label = first(tbl$endpoint_label),
    endpoint_order = first(tbl$endpoint_order),
    clinical_attribute = first(tbl$clinical_attribute),
    reference_group = reference_group,
    case_group = case_group,
    test_status = test_status,
    n_reference = length(reference_values),
    reference_median_burden = median(reference_values, na.rm = TRUE),
    reference_q25_burden = unname(quantile(
      reference_values,
      0.25,
      na.rm = TRUE)),
    reference_q75_burden = unname(quantile(
      reference_values,
      0.75,
      na.rm = TRUE)),
    n_case = length(case_values),
    case_median_burden = median(case_values, na.rm = TRUE),
    case_q25_burden = unname(quantile(case_values, 0.25, na.rm = TRUE)),
    case_q75_burden = unname(quantile(case_values, 0.75, na.rm = TRUE)),
    median_burden_difference = case_median_burden -
      reference_median_burden,
    p_value = if (is.null(wilcoxon_result)) {
      NA_real_
    } else {
      wilcoxon_result$p.value
    })
}

run_cbio_candidate_burden_associations <- function(
    burden_tbl,
    endpoint_rows) {

  burden_tbl %>%
    inner_join(
      endpoint_rows %>%
        filter(!is.na(endpoint_group)),
      by = c("study_id", "sample_id"),
      relationship = "many-to-many") %>%
    group_by(endpoint_id) %>%
    group_modify(~ summarise_cbio_burden_endpoint_test(.x)) %>%
    ungroup() %>%
    group_by(analysis_family) %>%
    mutate(q_value_analysis_family = p_adjust_with_na(p_value)) %>%
    ungroup() %>%
    arrange(endpoint_order)
}

prepare_cbio_candidate_tmb_data <- function(
    burden_tbl,
    endpoint_rows) {

  sample_type_rows <- endpoint_rows %>%
    filter(
      endpoint_id == "sample_type_pooled",
      !is.na(endpoint_group)) %>%
    select(study_id, sample_id, endpoint_group) %>%
    distinct()

  pooled_data <- burden_tbl %>%
    inner_join(
      sample_type_rows,
      by = c("study_id", "sample_id")) %>%
    mutate(
      analysis_scope = "Pooled",
      analysis_scope_order = 1L)

  bind_rows(
    pooled_data,
    pooled_data %>%
      filter(study_id == "cscc_hgsc_bcm_2014") %>%
      mutate(
        analysis_scope = "HGSC only",
        analysis_scope_order = 2L)) %>%
    mutate(
      candidate_status = if_else(
        n_mutated_candidate_genes > 0L,
        "Candidate-mutated",
        "No candidate mutation"),
      candidate_status = factor(
        candidate_status,
        levels = c(
          "No candidate mutation",
          "Candidate-mutated")),
      endpoint_group = factor(
        endpoint_group,
        levels = c("Primary", "Metastasis")),
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "HGSC only")),
      study_label = recode(
        study_id,
        cscc_dfarber_2015 = "DFCI",
        cscc_hgsc_bcm_2014 = "HGSC",
        cscc_ranson_2022 = "Ranson",
        cscc_ucsf_2021 = "UCSF"),
      study_label = factor(
        study_label,
        levels = c("DFCI", "HGSC", "Ranson", "UCSF")),
      log10_tmb_nonsynonymous = log10(tmb_nonsynonymous),
      tmb_definition =
        "cBioPortal clinical attribute TMB_NONSYNONYMOUS (mutations/Mb)",
      candidate_mutation_definition = paste(
        "At least one of the 30 final candidates has any nonsynonymous",
        "mutation in a mapped human orthologue; exact variants are not matched.")) %>%
    arrange(analysis_scope_order, study_id, sample_id)
}

summarise_cbio_candidate_tmb <- function(tmb_data) {
  tmb_data %>%
    group_by(
      analysis_scope,
      analysis_scope_order,
      endpoint_group,
      candidate_status) %>%
    summarise(
      n_samples = n_distinct(paste(study_id, sample_id, sep = "::")),
      tmb_min = min(tmb_nonsynonymous),
      tmb_q25 = unname(quantile(tmb_nonsynonymous, 0.25)),
      tmb_median = median(tmb_nonsynonymous),
      tmb_q75 = unname(quantile(tmb_nonsynonymous, 0.75)),
      tmb_max = max(tmb_nonsynonymous),
      .groups = "drop") %>%
    arrange(
      analysis_scope_order,
      endpoint_group,
      candidate_status)
}

run_cbio_candidate_tmb_wilcoxon <- function(
    tmb_data,
    analysis_scope,
    analysis_id,
    test_role) {

  test_data <- tmb_data %>%
    filter(
      as.character(.data$analysis_scope) == .env$analysis_scope,
      candidate_status == "Candidate-mutated")

  reference_values <- test_data$tmb_nonsynonymous[
    test_data$endpoint_group == "Primary"]
  case_values <- test_data$tmb_nonsynonymous[
    test_data$endpoint_group == "Metastasis"]
  test_status <- if (
      length(reference_values) > 0L &&
      length(case_values) > 0L) {
    "tested"
  } else {
    "not_testable_missing_disease_group"
  }
  test_result <- if (test_status == "tested") {
    suppressWarnings(wilcox.test(
      case_values,
      reference_values,
      alternative = "two.sided",
      exact = FALSE))
  } else {
    NULL
  }

  tibble(
    analysis_id = analysis_id,
    test_role = test_role,
    analysis_scope = analysis_scope,
    population = "Candidate-mutated tumours",
    test_method = "Two-sided Wilcoxon rank-sum test",
    model_formula = NA_character_,
    study_adjustment = "None",
    informative_studies = collapse_unique(
      test_data$study_id,
      n = Inf,
      sep = ";"),
    reference_group = "Primary",
    case_group = "Metastasis",
    n_reference = length(reference_values),
    n_case = length(case_values),
    reference_median_tmb = median(reference_values),
    reference_q25_tmb = unname(quantile(reference_values, 0.25)),
    reference_q75_tmb = unname(quantile(reference_values, 0.75)),
    case_median_tmb = median(case_values),
    case_q25_tmb = unname(quantile(case_values, 0.25)),
    case_q75_tmb = unname(quantile(case_values, 0.75)),
    estimate_type = "Median TMB difference (metastasis - primary)",
    estimate = median(case_values) - median(reference_values),
    confidence_interval_low = NA_real_,
    confidence_interval_high = NA_real_,
    p_value = if (is.null(test_result)) {
      NA_real_
    } else {
      test_result$p.value
    },
    test_status = test_status,
    interpretation_note = if_else(
      analysis_scope == "Pooled",
      paste(
        "Descriptive pooled comparison; disease status is strongly",
        "confounded with study."),
      "Within-HGSC sensitivity comparison."))
}

run_cbio_candidate_tmb_adjusted_model <- function(tmb_data) {
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("The sandwich package is required for HC3 robust standard errors.")
  }

  model_data <- tmb_data %>%
    filter(
      analysis_scope == "Pooled",
      candidate_status == "Candidate-mutated") %>%
    mutate(
      endpoint_group = factor(
        endpoint_group,
        levels = c("Primary", "Metastasis")),
      study_id = factor(
        study_id,
        levels = c(
          "cscc_hgsc_bcm_2014",
          "cscc_dfarber_2015",
          "cscc_ranson_2022",
          "cscc_ucsf_2021")))

  model <- lm(
    log10_tmb_nonsynonymous ~ endpoint_group + study_id,
    data = model_data)
  design_matrix <- model.matrix(model)
  model_is_full_rank <- qr(design_matrix)$rank == ncol(design_matrix)

  if (!model_is_full_rank) {
    stop("The study-adjusted candidate TMB model is rank deficient.")
  }

  robust_vcov <- sandwich::vcovHC(model, type = "HC3")
  coefficient_name <- "endpoint_groupMetastasis"
  coefficient_estimate <- unname(coef(model)[[coefficient_name]])
  coefficient_se <- sqrt(robust_vcov[
    coefficient_name,
    coefficient_name])
  model_df <- df.residual(model)
  critical_value <- qt(0.975, df = model_df)
  coefficient_t <- coefficient_estimate / coefficient_se
  p_value <- 2 * pt(abs(coefficient_t), df = model_df, lower.tail = FALSE)
  coefficient_ci <- coefficient_estimate +
    c(-1, 1) * critical_value * coefficient_se

  reference_values <- model_data$tmb_nonsynonymous[
    model_data$endpoint_group == "Primary"]
  case_values <- model_data$tmb_nonsynonymous[
    model_data$endpoint_group == "Metastasis"]

  tibble(
    analysis_id = "study_adjusted_log10_tmb",
    test_role = "Primary adjusted analysis",
    analysis_scope = "Pooled",
    population = "Candidate-mutated tumours",
    test_method = paste(
      "Linear regression of log10 TMB with study fixed effects",
      "and HC3 robust standard errors"),
    model_formula =
      "log10(TMB_NONSYNONYMOUS) ~ primary/metastatic status + study",
    study_adjustment = "Study fixed effects",
    informative_studies = "cscc_hgsc_bcm_2014",
    reference_group = "Primary",
    case_group = "Metastasis",
    n_reference = length(reference_values),
    n_case = length(case_values),
    reference_median_tmb = median(reference_values),
    reference_q25_tmb = unname(quantile(reference_values, 0.25)),
    reference_q75_tmb = unname(quantile(reference_values, 0.75)),
    case_median_tmb = median(case_values),
    case_q25_tmb = unname(quantile(case_values, 0.25)),
    case_q75_tmb = unname(quantile(case_values, 0.75)),
    estimate_type = "Metastatic / primary geometric-mean TMB ratio",
    estimate = 10^coefficient_estimate,
    confidence_interval_low = 10^coefficient_ci[[1]],
    confidence_interval_high = 10^coefficient_ci[[2]],
    p_value = p_value,
    test_status = "tested",
    interpretation_note = paste(
      "The disease coefficient is identified primarily by HGSC because",
      "DFCI and Ranson contain only metastatic tumours and UCSF contains",
      "only primary tumours."),
    model_n = nobs(model),
    model_residual_df = model_df,
    model_design_columns = ncol(design_matrix),
    model_design_rank = qr(design_matrix)$rank,
    model_is_full_rank = model_is_full_rank,
    log10_coefficient = coefficient_estimate,
    log10_coefficient_hc3_se = coefficient_se)
}

run_cbio_candidate_tmb_tests <- function(tmb_data) {
  adjusted_result <- run_cbio_candidate_tmb_adjusted_model(tmb_data)
  pooled_wilcoxon <- run_cbio_candidate_tmb_wilcoxon(
    tmb_data = tmb_data,
    analysis_scope = "Pooled",
    analysis_id = "pooled_candidate_mutated_wilcoxon",
    test_role = "Pooled sensitivity analysis")
  hgsc_wilcoxon <- run_cbio_candidate_tmb_wilcoxon(
    tmb_data = tmb_data,
    analysis_scope = "HGSC only",
    analysis_id = "hgsc_candidate_mutated_wilcoxon",
    test_role = "Within-study sensitivity analysis")

  bind_rows(
    adjusted_result,
    pooled_wilcoxon,
    hgsc_wilcoxon) %>%
    mutate(
      multiple_testing_adjustment =
        "None: one prespecified primary model and two labelled sensitivity tests") %>%
    arrange(match(
      analysis_id,
      c(
        "study_adjusted_log10_tmb",
        "pooled_candidate_mutated_wilcoxon",
        "hgsc_candidate_mutated_wilcoxon")))
}

summarise_cbio_candidate_mutation_counts <- function(tmb_data) {
  tmb_data %>%
    group_by(
      analysis_scope,
      analysis_scope_order,
      endpoint_group) %>%
    summarise(
      n_samples = n_distinct(paste(study_id, sample_id, sep = "::")),
      n_samples_with_candidate_mutation = sum(
        n_candidate_nonsynonymous_mutations > 0L),
      candidate_mutation_count_min =
        min(n_candidate_nonsynonymous_mutations),
      candidate_mutation_count_q25 = unname(quantile(
        n_candidate_nonsynonymous_mutations,
        0.25)),
      candidate_mutation_count_median =
        median(n_candidate_nonsynonymous_mutations),
      candidate_mutation_count_q75 = unname(quantile(
        n_candidate_nonsynonymous_mutations,
        0.75)),
      candidate_mutation_count_max =
        max(n_candidate_nonsynonymous_mutations),
      .groups = "drop") %>%
    arrange(analysis_scope_order, endpoint_group)
}

run_cbio_candidate_mutation_count_test <- function(
    tmb_data,
    analysis_scope,
    analysis_id) {

  test_data <- tmb_data %>%
    filter(as.character(.data$analysis_scope) == .env$analysis_scope)
  reference_values <- test_data$n_candidate_nonsynonymous_mutations[
    test_data$endpoint_group == "Primary"]
  case_values <- test_data$n_candidate_nonsynonymous_mutations[
    test_data$endpoint_group == "Metastasis"]
  test_result <- suppressWarnings(wilcox.test(
    case_values,
    reference_values,
    alternative = "two.sided",
    exact = FALSE))

  tibble(
    analysis_id = analysis_id,
    analysis_scope = analysis_scope,
    population = "All tumours, including zero candidate mutations",
    test_method = "Two-sided Wilcoxon rank-sum test",
    reference_group = "Primary",
    case_group = "Metastasis",
    n_reference = length(reference_values),
    n_case = length(case_values),
    reference_median_candidate_mutations = median(reference_values),
    reference_q25_candidate_mutations = unname(quantile(
      reference_values,
      0.25)),
    reference_q75_candidate_mutations = unname(quantile(
      reference_values,
      0.75)),
    case_median_candidate_mutations = median(case_values),
    case_q25_candidate_mutations = unname(quantile(case_values, 0.25)),
    case_q75_candidate_mutations = unname(quantile(case_values, 0.75)),
    estimate_type =
      "Median count difference (metastasis - primary)",
    estimate = median(case_values) - median(reference_values),
    p_value = test_result$p.value,
    multiple_testing_adjustment =
      "None: pooled descriptive analysis and HGSC sensitivity analysis",
    interpretation_note = if_else(
      analysis_scope == "Pooled",
      paste(
        "Descriptive pooled comparison; disease status is strongly",
        "confounded with study."),
      "Within-HGSC sensitivity comparison."),
    test_status = "tested")
}

run_cbio_candidate_mutation_count_tests <- function(tmb_data) {
  bind_rows(
    run_cbio_candidate_mutation_count_test(
      tmb_data,
      analysis_scope = "Pooled",
      analysis_id = "pooled_candidate_mutation_count_wilcoxon"),
    run_cbio_candidate_mutation_count_test(
      tmb_data,
      analysis_scope = "HGSC only",
      analysis_id = "hgsc_candidate_mutation_count_wilcoxon"))
}

prepare_cbio_candidate_gene_count_tmb_correlation_data <- function(tmb_data) {
  tmb_data %>%
    filter(analysis_scope == "Pooled") %>%
    transmute(
      study_id,
      sample_id,
      cbio_sample_uid,
      n_mutated_candidate_genes,
      tmb_nonsynonymous,
      tmb_definition,
      candidate_mutation_definition) %>%
    distinct()
}

run_cbio_candidate_gene_count_tmb_spearman_test <- function(
    correlation_data) {

  candidate_gene_counts <- correlation_data$n_mutated_candidate_genes
  tmb_values <- correlation_data$tmb_nonsynonymous
  test_status <- if (
      n_distinct(candidate_gene_counts) > 1L &&
      n_distinct(tmb_values) > 1L) {
    "tested"
  } else {
    "not_testable_no_variation"
  }
  test_result <- if (test_status == "tested") {
    suppressWarnings(cor.test(
      candidate_gene_counts,
      tmb_values,
      method = "spearman",
      alternative = "two.sided",
      exact = FALSE))
  } else {
    NULL
  }

  tibble(
    analysis_id = "pooled_candidate_gene_count_tmb_spearman",
    analysis_scope = "Pooled",
    test_role = "Primary analysis",
    population = paste(
      "All cBioPortal tumours, including samples with zero mutated",
      "candidate genes"),
    x_variable = "Number of distinct mutated candidate genes",
    y_variable = "TMB_NONSYNONYMOUS (mutations/Mb)",
    test_method = "Two-sided Spearman rank correlation",
    n_samples = nrow(correlation_data),
    n_zero_candidate_genes = sum(candidate_gene_counts == 0L),
    min_candidate_genes = min(candidate_gene_counts),
    max_candidate_genes = max(candidate_gene_counts),
    spearman_rho = if (is.null(test_result)) {
      NA_real_
    } else {
      unname(test_result$estimate)
    },
    p_value = if (is.null(test_result)) {
      NA_real_
    } else {
      pmin(1, pmax(0, test_result$p.value))
    },
    test_status,
    multiple_testing_adjustment =
      "None: one prespecified pooled correlation test",
    interpretation_note = paste(
      "Candidate mutations contribute to nonsynonymous TMB, so the",
      "correlation is partly mathematically coupled. Candidate mutation",
      "calls retain panel-aware profiling."))
}

plot_cbio_candidate_gene_count_tmb_correlation <- function(
    correlation_data,
    correlation_test) {

  annotation_label <- paste0(
    "Spearman rho = ",
    formatC(correlation_test$spearman_rho, format = "f", digits = 2),
    "\nP ",
    format_clinical_probability(correlation_test$p_value))
  max_candidate_genes <- max(
    correlation_data$n_mutated_candidate_genes)

  ggplot(
    correlation_data,
    aes(
      x = n_mutated_candidate_genes,
      y = tmb_nonsynonymous)) +
    geom_point(
      shape = 21,
      fill = "#BC667E",
      colour = "black",
      stroke = 0.2,
      size = 1.5,
      alpha = 0.65,
      position = position_jitter(
        width = 0.12,
        height = 0,
        seed = 20260729L)) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annotation_label,
      hjust = -0.08,
      vjust = 1.1,
      size = 7 / .pt) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      limits = c(-0.45, max_candidate_genes + 0.45),
      expand = expansion(mult = 0)) +
    scale_y_log10(
      breaks = scales::breaks_log(n = 6),
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.05, 0.2))) +
    labs(
      x = "Mutated candidate genes (n)",
      y = "Nonsynonymous TMB (mutations/Mb)") +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none")
}

validate_cbio_candidate_gene_count_tmb_correlation <- function(
    correlation_data,
    correlation_test) {

  if (nrow(correlation_data) != 176L ||
      n_distinct(correlation_data$cbio_sample_uid) != 176L ||
      sum(correlation_data$n_mutated_candidate_genes == 0L) != 32L ||
      any(
        correlation_data$n_mutated_candidate_genes < 0L |
        correlation_data$n_mutated_candidate_genes > 30L |
        correlation_data$n_mutated_candidate_genes !=
          as.integer(correlation_data$n_mutated_candidate_genes)) ||
      any(
        is.na(correlation_data$tmb_nonsynonymous) |
        !is.finite(correlation_data$tmb_nonsynonymous) |
        correlation_data$tmb_nonsynonymous <= 0)) {
    stop("cBioPortal candidate-gene/TMB correlation data failed validation.")
  }

  if (nrow(correlation_test) != 1L ||
      correlation_test$test_status != "tested" ||
      correlation_test$n_samples != 176L ||
      correlation_test$n_zero_candidate_genes != 32L ||
      !is.finite(correlation_test$spearman_rho) ||
      abs(correlation_test$spearman_rho) > 1 ||
      !is.finite(correlation_test$p_value) ||
      correlation_test$p_value < 0 |
      correlation_test$p_value > 1) {
    stop("cBioPortal candidate-gene/TMB Spearman test failed validation.")
  }

  invisible(TRUE)
}

prepare_cbio_candidate_recurrence_boxplot_data <- function(
    burden_tbl,
    endpoint_rows) {

  recurrence_rows <- endpoint_rows %>%
    filter(
      endpoint_id == "recurrence_pooled",
      !is.na(endpoint_group)) %>%
    select(study_id, sample_id, endpoint_group) %>%
    distinct()
  pooled_data <- burden_tbl %>%
    inner_join(
      recurrence_rows,
      by = c("study_id", "sample_id")) %>%
    mutate(
      analysis_scope = "Pooled",
      analysis_scope_order = 1L)

  bind_rows(
    pooled_data,
    pooled_data %>%
      filter(study_id == "cscc_dfarber_2015") %>%
      mutate(
        analysis_scope = "DFCI only",
        analysis_scope_order = 2L),
    pooled_data %>%
      filter(study_id == "cscc_hgsc_bcm_2014") %>%
      mutate(
        analysis_scope = "HGSC only",
        analysis_scope_order = 3L)) %>%
    mutate(
      candidate_status = if_else(
        n_mutated_candidate_genes > 0L,
        "Candidate-mutated",
        "No candidate mutation"),
      candidate_status = factor(
        candidate_status,
        levels = c(
          "No candidate mutation",
          "Candidate-mutated")),
      endpoint_group = factor(
        endpoint_group,
        levels = c("Disease free", "Recurred / progressed")),
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "DFCI only", "HGSC only")),
      log10_tmb_nonsynonymous = log10(tmb_nonsynonymous),
      recurrence_definition = paste(
        "DFCI DFS_STATUS and HGSC",
        "DISEASE_RECURRENCE_OR_PERSISTENCE, normalized to",
        "Disease free versus Recurred / progressed."),
      candidate_mutation_definition = paste(
        "Any nonsynonymous mutation in the 30 final candidates or mapped",
        "human orthologues; exact variants are not matched.")) %>%
    arrange(analysis_scope_order, study_id, sample_id)
}

summarise_cbio_candidate_recurrence_tmb <- function(recurrence_data) {
  recurrence_data %>%
    group_by(
      analysis_scope,
      analysis_scope_order,
      endpoint_group,
      candidate_status) %>%
    summarise(
      n_samples = n_distinct(paste(study_id, sample_id, sep = "::")),
      tmb_min = min(tmb_nonsynonymous),
      tmb_q25 = unname(quantile(tmb_nonsynonymous, 0.25)),
      tmb_median = median(tmb_nonsynonymous),
      tmb_q75 = unname(quantile(tmb_nonsynonymous, 0.75)),
      tmb_max = max(tmb_nonsynonymous),
      .groups = "drop") %>%
    arrange(
      analysis_scope_order,
      endpoint_group,
      candidate_status)
}

run_cbio_candidate_recurrence_tmb_wilcoxon <- function(
    recurrence_data,
    analysis_scope,
    analysis_id,
    test_role) {

  test_data <- recurrence_data %>%
    filter(
      as.character(.data$analysis_scope) == .env$analysis_scope,
      candidate_status == "Candidate-mutated")
  reference_values <- test_data$tmb_nonsynonymous[
    test_data$endpoint_group == "Disease free"]
  case_values <- test_data$tmb_nonsynonymous[
    test_data$endpoint_group == "Recurred / progressed"]
  test_result <- suppressWarnings(wilcox.test(
    case_values,
    reference_values,
    alternative = "two.sided",
    exact = FALSE))

  tibble(
    analysis_id = analysis_id,
    test_role = test_role,
    analysis_scope = analysis_scope,
    population = "Candidate-mutated tumours",
    test_method = "Two-sided Wilcoxon rank-sum test",
    model_formula = NA_character_,
    study_adjustment = "None",
    reference_group = "Disease free",
    case_group = "Recurred / progressed",
    n_reference = length(reference_values),
    n_case = length(case_values),
    reference_median_tmb = median(reference_values),
    reference_q25_tmb = unname(quantile(reference_values, 0.25)),
    reference_q75_tmb = unname(quantile(reference_values, 0.75)),
    case_median_tmb = median(case_values),
    case_q25_tmb = unname(quantile(case_values, 0.25)),
    case_q75_tmb = unname(quantile(case_values, 0.75)),
    estimate_type =
      "Median TMB difference (recurred/progressed - disease free)",
    estimate = median(case_values) - median(reference_values),
    confidence_interval_low = NA_real_,
    confidence_interval_high = NA_real_,
    p_value = test_result$p.value,
    test_status = "tested",
    interpretation_note = if_else(
      analysis_scope == "Pooled",
      "Pooled sensitivity comparison without study adjustment.",
      paste0("Within-", str_remove(analysis_scope, " only"), " comparison.")))
}

run_cbio_candidate_recurrence_tmb_adjusted_model <- function(
    recurrence_data) {

  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("The sandwich package is required for HC3 robust standard errors.")
  }
  model_data <- recurrence_data %>%
    filter(
      analysis_scope == "Pooled",
      candidate_status == "Candidate-mutated") %>%
    mutate(
      endpoint_group = factor(
        endpoint_group,
        levels = c("Disease free", "Recurred / progressed")),
      study_id = factor(
        study_id,
        levels = c("cscc_dfarber_2015", "cscc_hgsc_bcm_2014")))
  model <- lm(
    log10_tmb_nonsynonymous ~ endpoint_group + study_id,
    data = model_data)
  design_matrix <- model.matrix(model)
  model_is_full_rank <- qr(design_matrix)$rank == ncol(design_matrix)
  if (!model_is_full_rank) {
    stop("The study-adjusted recurrence TMB model is rank deficient.")
  }

  robust_vcov <- sandwich::vcovHC(model, type = "HC3")
  coefficient_name <- "endpoint_groupRecurred / progressed"
  coefficient_estimate <- unname(coef(model)[[coefficient_name]])
  coefficient_se <- sqrt(robust_vcov[
    coefficient_name,
    coefficient_name])
  model_df <- df.residual(model)
  critical_value <- qt(0.975, df = model_df)
  coefficient_t <- coefficient_estimate / coefficient_se
  p_value <- 2 * pt(abs(coefficient_t), df = model_df, lower.tail = FALSE)
  coefficient_ci <- coefficient_estimate +
    c(-1, 1) * critical_value * coefficient_se
  reference_values <- model_data$tmb_nonsynonymous[
    model_data$endpoint_group == "Disease free"]
  case_values <- model_data$tmb_nonsynonymous[
    model_data$endpoint_group == "Recurred / progressed"]

  tibble(
    analysis_id = "study_adjusted_log10_tmb_recurrence",
    test_role = "Primary adjusted analysis",
    analysis_scope = "Pooled",
    population = "Candidate-mutated tumours",
    test_method = paste(
      "Linear regression of log10 TMB with study fixed effects",
      "and HC3 robust standard errors"),
    model_formula = paste(
      "log10(TMB_NONSYNONYMOUS) ~ recurrence/progression status + study"),
    study_adjustment = "Study fixed effects",
    reference_group = "Disease free",
    case_group = "Recurred / progressed",
    n_reference = length(reference_values),
    n_case = length(case_values),
    reference_median_tmb = median(reference_values),
    reference_q25_tmb = unname(quantile(reference_values, 0.25)),
    reference_q75_tmb = unname(quantile(reference_values, 0.75)),
    case_median_tmb = median(case_values),
    case_q25_tmb = unname(quantile(case_values, 0.25)),
    case_q75_tmb = unname(quantile(case_values, 0.75)),
    estimate_type =
      "Recurred/progressed / disease-free geometric-mean TMB ratio",
    estimate = 10^coefficient_estimate,
    confidence_interval_low = 10^coefficient_ci[[1]],
    confidence_interval_high = 10^coefficient_ci[[2]],
    p_value = p_value,
    test_status = "tested",
    interpretation_note = paste(
      "Both DFCI and HGSC contain disease-free and",
      "recurrent/progressed tumours."),
    model_n = nobs(model),
    model_residual_df = model_df,
    model_design_columns = ncol(design_matrix),
    model_design_rank = qr(design_matrix)$rank,
    model_is_full_rank = model_is_full_rank,
    log10_coefficient = coefficient_estimate,
    log10_coefficient_hc3_se = coefficient_se)
}

run_cbio_candidate_recurrence_tmb_tests <- function(recurrence_data) {
  bind_rows(
    run_cbio_candidate_recurrence_tmb_adjusted_model(recurrence_data),
    run_cbio_candidate_recurrence_tmb_wilcoxon(
      recurrence_data,
      analysis_scope = "Pooled",
      analysis_id = "pooled_candidate_mutated_recurrence_wilcoxon",
      test_role = "Pooled sensitivity analysis"),
    run_cbio_candidate_recurrence_tmb_wilcoxon(
      recurrence_data,
      analysis_scope = "DFCI only",
      analysis_id = "dfci_candidate_mutated_recurrence_wilcoxon",
      test_role = "Within-study sensitivity analysis"),
    run_cbio_candidate_recurrence_tmb_wilcoxon(
      recurrence_data,
      analysis_scope = "HGSC only",
      analysis_id = "hgsc_candidate_mutated_recurrence_wilcoxon",
      test_role = "Within-study sensitivity analysis")) %>%
    mutate(
      multiple_testing_adjustment = paste(
        "None: one prespecified adjusted model and three labelled",
        "sensitivity analyses"))
}

summarise_cbio_candidate_recurrence_mutation_counts <- function(
    recurrence_data) {

  recurrence_data %>%
    group_by(
      analysis_scope,
      analysis_scope_order,
      endpoint_group) %>%
    summarise(
      n_samples = n_distinct(paste(study_id, sample_id, sep = "::")),
      n_samples_with_candidate_mutation = sum(
        n_candidate_nonsynonymous_mutations > 0L),
      candidate_mutation_count_min =
        min(n_candidate_nonsynonymous_mutations),
      candidate_mutation_count_q25 = unname(quantile(
        n_candidate_nonsynonymous_mutations,
        0.25)),
      candidate_mutation_count_median =
        median(n_candidate_nonsynonymous_mutations),
      candidate_mutation_count_q75 = unname(quantile(
        n_candidate_nonsynonymous_mutations,
        0.75)),
      candidate_mutation_count_max =
        max(n_candidate_nonsynonymous_mutations),
      .groups = "drop") %>%
    arrange(analysis_scope_order, endpoint_group)
}

run_cbio_candidate_recurrence_mutation_count_test <- function(
    recurrence_data,
    analysis_scope,
    analysis_id) {

  test_data <- recurrence_data %>%
    filter(as.character(.data$analysis_scope) == .env$analysis_scope)
  reference_values <- test_data$n_candidate_nonsynonymous_mutations[
    test_data$endpoint_group == "Disease free"]
  case_values <- test_data$n_candidate_nonsynonymous_mutations[
    test_data$endpoint_group == "Recurred / progressed"]
  test_result <- suppressWarnings(wilcox.test(
    case_values,
    reference_values,
    alternative = "two.sided",
    exact = FALSE))

  tibble(
    analysis_id = analysis_id,
    analysis_scope = analysis_scope,
    population = "All tumours, including zero candidate mutations",
    test_method = "Two-sided Wilcoxon rank-sum test",
    reference_group = "Disease free",
    case_group = "Recurred / progressed",
    n_reference = length(reference_values),
    n_case = length(case_values),
    reference_median_candidate_mutations = median(reference_values),
    reference_q25_candidate_mutations =
      unname(quantile(reference_values, 0.25)),
    reference_q75_candidate_mutations =
      unname(quantile(reference_values, 0.75)),
    case_median_candidate_mutations = median(case_values),
    case_q25_candidate_mutations = unname(quantile(case_values, 0.25)),
    case_q75_candidate_mutations = unname(quantile(case_values, 0.75)),
    estimate_type = paste(
      "Median count difference",
      "(recurred/progressed - disease free)"),
    estimate = median(case_values) - median(reference_values),
    p_value = test_result$p.value,
    multiple_testing_adjustment =
      "None: pooled and two labelled within-study comparisons",
    interpretation_note = if_else(
      analysis_scope == "Pooled",
      "Pooled comparison without study adjustment.",
      paste0("Within-", str_remove(analysis_scope, " only"), " comparison.")),
    test_status = "tested")
}

run_cbio_candidate_recurrence_mutation_count_tests <- function(
    recurrence_data) {

  bind_rows(
    run_cbio_candidate_recurrence_mutation_count_test(
      recurrence_data,
      analysis_scope = "Pooled",
      analysis_id = "pooled_recurrence_candidate_mutation_count_wilcoxon"),
    run_cbio_candidate_recurrence_mutation_count_test(
      recurrence_data,
      analysis_scope = "DFCI only",
      analysis_id = "dfci_recurrence_candidate_mutation_count_wilcoxon"),
    run_cbio_candidate_recurrence_mutation_count_test(
      recurrence_data,
      analysis_scope = "HGSC only",
      analysis_id = "hgsc_recurrence_candidate_mutation_count_wilcoxon"))
}

format_clinical_probability <- function(x) {
  case_when(
    is.na(x) ~ "not testable",
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3))
}

plot_cbio_candidate_clinical_associations <- function(
    association_tbl,
    gene_levels,
    endpoint_definitions) {

  effect_values <- association_tbl$corrected_log2_odds_ratio[
    association_tbl$test_status == "tested"]
  effect_values <- effect_values[is.finite(effect_values)]
  effect_limit <- if (length(effect_values) == 0) {
    2
  } else {
    max(2, min(6, ceiling(max(abs(effect_values)))))
  }

  plot_df <- association_tbl %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(gene_levels)),
      endpoint_label = factor(
        endpoint_label,
        levels = endpoint_definitions %>%
          arrange(endpoint_order) %>%
          pull(endpoint_label)),
      candidate_group_plot = factor(
        candidate_group_display(candidate_group),
        levels = CBIO_CANDIDATE_GROUP_LEVELS),
      evidence_strength = if_else(
        test_status == "tested",
        pmin(-log10(pmax(p_value, .Machine$double.xmin)), 5),
        NA_real_))

  ggplot(plot_df, aes(x = endpoint_label, y = gene_label)) +
    geom_tile(
      fill = "grey94",
      colour = "white",
      linewidth = 0.25) +
    geom_point(
      data = plot_df %>% filter(test_status == "tested"),
      aes(
        fill = corrected_log2_odds_ratio,
        size = evidence_strength),
      shape = 21,
      colour = "black",
      stroke = 0.25) +
    geom_text(
      data = plot_df %>%
        filter(
          test_status == "tested",
          q_value_analysis_family <= 0.05),
      label = "*",
      colour = "black",
      size = 7 / .pt,
      vjust = 0.35) +
    facet_grid(
      candidate_group_plot ~ .,
      scales = "free_y",
      space = "free_y") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-effect_limit, effect_limit),
      oob = scales::squish,
      name = expression(log[2]("odds ratio"))) +
    scale_size_continuous(
      range = c(1.2, 4),
      limits = c(0, 5),
      breaks = c(1, 3, 5),
      name = expression(-log[10](italic(P)))) +
    labs(x = "", y = "") +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      strip.text.y = element_text(angle = 0),
      panel.grid = element_blank(),
      legend.position = "right")
}

plot_cbio_primary_metastatic_prevalence <- function(
    prevalence_tbl,
    gene_levels) {

  plot_df <- prevalence_tbl %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(gene_levels)),
      candidate_group_plot = factor(
        candidate_group_display(candidate_group),
        levels = CBIO_CANDIDATE_GROUP_LEVELS))

  ggplot(
    plot_df,
    aes(
      x = prevalence_pct,
      y = gene_label,
      fill = endpoint_group)) +
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.25,
      size = 2.2,
      position = position_dodge(width = 0.55)) +
    facet_grid(
      candidate_group_plot ~ comparison_cohort,
      scales = "free_y",
      space = "free_y") +
    scale_fill_manual(
      values = c(
        Primary = "#4C78A8",
        Metastasis = "#D62728"),
      drop = FALSE) +
    scale_x_continuous(
      breaks = c(0, 50, 100),
      labels = function(x) paste0(x, "%"),
      limits = c(0, 100),
      expand = expansion(mult = c(0.02, 0.04))) +
    labs(
      x = "Samples with a nonsynonymous mutation",
      y = "",
      fill = "Sample type") +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 7),
      strip.text.y = element_text(angle = 0),
      legend.position = "bottom")
}

prepare_cbio_primary_metastatic_highlight_data <- function(association_tbl) {
  association_tbl %>%
    filter(endpoint_id == "sample_type_pooled") %>%
    mutate(
      association_direction = case_when(
        test_status != "tested" ~ "Not testable",
        q_value_within_endpoint <= 0.05 &
          corrected_log2_odds_ratio < 0 ~ "Primary-associated",
        q_value_within_endpoint <= 0.05 &
          corrected_log2_odds_ratio > 0 ~ "Metastatic-associated",
        TRUE ~ "Not significant"),
      association_direction = factor(
        association_direction,
        levels = c(
          "Primary-associated",
          "Metastatic-associated",
          "Not significant",
          "Not testable")),
      evidence_strength = if_else(
        is.na(q_value_within_endpoint),
        0,
        pmin(
          -log10(pmax(q_value_within_endpoint, .Machine$double.xmin)),
          6)),
      primary_mutation_label = paste0(
        n_reference_mutated,
        "/",
        n_reference_profiled),
      metastatic_mutation_label = paste0(
        n_case_mutated,
        "/",
        n_case_profiled),
      prevalence_label = paste0(
        formatC(reference_prevalence_pct, format = "f", digits = 1),
        "% vs ",
        formatC(case_prevalence_pct, format = "f", digits = 1),
        "%"),
      candidate_group_plot = factor(
        candidate_group_display(candidate_group),
        levels = CBIO_CANDIDATE_GROUP_LEVELS))
}

plot_cbio_primary_metastatic_candidate_volcano <- function(plot_data) {
  plot_df <- plot_data %>%
    mutate(
      volcano_y = if_else(
        is.na(q_value_within_endpoint),
        0,
        -log10(pmax(
          q_value_within_endpoint,
          .Machine$double.xmin))),
      volcano_label = sub("^.* / ", "", gene_label))

  association_colours <- c(
    "Primary-associated" = "#4C78A8",
    "Metastatic-associated" = "#D62728",
    "Not significant" = "grey72",
    "Not testable" = "white")
  association_breaks <- c(
    "Primary-associated",
    "Metastatic-associated",
    "Not significant",
    "Not testable")
  association_breaks <- association_breaks[
    association_breaks %in% as.character(plot_df$association_direction)]
  x_limit <- max(
    1,
    ceiling(max(
      abs(plot_df$corrected_log2_odds_ratio),
      na.rm = TRUE)))
  candidate_shapes <- setNames(
    c(21, 22, 24),
    CBIO_CANDIDATE_GROUP_LEVELS)
  n_not_testable <- sum(plot_df$association_direction == "Not testable")

  ggplot(
    plot_df,
    aes(
      x = corrected_log2_odds_ratio,
      y = volcano_y)) +
    geom_hline(
      yintercept = -log10(0.05),
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.3) +
    geom_vline(
      xintercept = 0,
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.3) +
    geom_point(
      aes(
        fill = association_direction,
        shape = candidate_group_plot),
      colour = "black",
      stroke = 0.35,
      size = 3) +
    ggrepel::geom_text_repel(
      data = plot_df %>%
        filter(association_direction %in% c(
          "Primary-associated",
          "Metastatic-associated")),
      aes(label = volcano_label),
      size = 7 / .pt,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE) +
    scale_fill_manual(
      values = association_colours,
      breaks = association_breaks,
      drop = TRUE,
      name = "FDR association") +
    scale_shape_manual(
      values = candidate_shapes,
      drop = FALSE,
      name = "Candidate group") +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.15))) +
    labs(
      x = expression(
        "Primary-associated" %<-% log[2]("corrected odds ratio") %->%
          "Metastatic-associated"),
      y = expression(-log[10]("BH Q")),
      caption = str_wrap(
        paste0(
          "Dashed horizontal line: BH Q = 0.05; labels mark ",
          "FDR-significant candidates. ",
          n_not_testable,
          "/",
          nrow(plot_df),
          " candidates are not testable and are shown as overlapping ",
          "white baseline points."),
        width = 75)) +
    my_theme +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          colour = "black",
          size = 3)),
      shape = guide_legend(
        override.aes = list(
          fill = "white",
          colour = "black",
          size = 3))) +
    theme(
      legend.position = "right")
}

plot_cbio_candidate_recurrence_volcano <- function(association_tbl) {
  label_genes <- association_tbl %>%
    filter(test_status == "tested") %>%
    arrange(p_value) %>%
    slice_head(n = 5L) %>%
    pull(gene_symbol) %>%
    union(
      association_tbl %>%
        filter(q_value <= 0.05) %>%
        pull(gene_symbol))
  plot_df <- association_tbl %>%
    mutate(
      volcano_y = if_else(
        is.na(p_value),
        NA_real_,
        -log10(pmax(p_value, .Machine$double.xmin))),
      volcano_label = if_else(
        gene_symbol %in% label_genes,
        sub("^.* / ", "", gene_label),
        NA_character_))
  association_colours <- c(
    "Disease-free-associated" = "#4C78A8",
    "Recurrence-associated" = "#D62728",
    "Not significant" = "grey72",
    "Not testable" = "white")
  candidate_shapes <- setNames(
    c(21, 22, 24),
    CBIO_CANDIDATE_GROUP_LEVELS)
  x_limit <- max(
    1,
    ceiling(max(
      abs(plot_df$corrected_log2_common_odds_ratio),
      na.rm = TRUE)))

  ggplot(
    plot_df,
    aes(
      x = corrected_log2_common_odds_ratio,
      y = volcano_y)) +
    geom_hline(
      yintercept = -log10(0.05),
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.3) +
    geom_vline(
      xintercept = 0,
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.3) +
    geom_point(
      data = plot_df %>%
        filter(
          !is.na(volcano_y),
          is.finite(corrected_log2_common_odds_ratio)),
      aes(
        fill = association_direction,
        shape = candidate_group_plot),
      colour = "black",
      stroke = 0.35,
      size = 3) +
    ggrepel::geom_text_repel(
      data = plot_df %>%
        filter(
          !is.na(volcano_label),
          !is.na(volcano_y),
          is.finite(corrected_log2_common_odds_ratio)),
      aes(label = volcano_label),
      size = 7 / .pt,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE) +
    scale_fill_manual(
      values = association_colours,
      drop = TRUE,
      name = "Study-adjusted FDR") +
    scale_shape_manual(
      values = candidate_shapes,
      drop = FALSE,
      name = "Candidate group") +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.18))) +
    labs(
      x = expression(
        "Disease-free-associated" %<-%
          log[2]("study-adjusted odds ratio") %->%
          "Recurrence-associated"),
      y = expression(-log[10]("P value")),
      caption = paste(
        "CMH tests control for DFCI/HGSC when both studies profile a gene;",
        "single-study candidates use within-study Fisher tests.",
        "\nFill indicates BH-FDR significance; the horizontal line marks",
        "nominal P = 0.05. Labels identify the five lowest P values.")) +
    my_theme +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          colour = "black",
          size = 3)),
      shape = guide_legend(
        override.aes = list(
          fill = "white",
          colour = "black",
          size = 3))) +
    theme(legend.position = "right")
}

plot_cbio_candidate_tmb_recurrence <- function(
    recurrence_data,
    recurrence_tests) {

  plot_df <- recurrence_data %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "DFCI only", "HGSC only")),
      endpoint_group = factor(
        endpoint_group,
        levels = c("Disease free", "Recurred / progressed")),
      candidate_status = factor(
        candidate_status,
        levels = c(
          "No candidate mutation",
          "Candidate-mutated")))
  annotation_df <- recurrence_tests %>%
    filter(analysis_id %in% c(
      "pooled_candidate_mutated_recurrence_wilcoxon",
      "dfci_candidate_mutated_recurrence_wilcoxon",
      "hgsc_candidate_mutated_recurrence_wilcoxon")) %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "DFCI only", "HGSC only")),
      annotation_label = paste0(
        "Candidate-mutated\nWilcoxon P ",
        format_clinical_probability(p_value)))
  adjusted_result <- recurrence_tests %>%
    filter(analysis_id == "study_adjusted_log10_tmb_recurrence")
  adjusted_caption <- paste0(
    "Study-adjusted recurrent/disease-free geometric-mean TMB ratio ",
    formatC(adjusted_result$estimate, format = "f", digits = 2),
    " (95% CI ",
    formatC(
      adjusted_result$confidence_interval_low,
      format = "f",
      digits = 2),
    "-",
    formatC(
      adjusted_result$confidence_interval_high,
      format = "f",
      digits = 2),
    "), HC3 P ",
    format_clinical_probability(adjusted_result$p_value),
    ".")

  ggplot(
    plot_df,
    aes(
      x = endpoint_group,
      y = tmb_nonsynonymous,
      fill = candidate_status)) +
    geom_boxplot(
      position = position_dodge(width = 0.72),
      width = 0.62,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.8) +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0.12,
        jitter.height = 0,
        dodge.width = 0.72,
        seed = 20260728L),
      shape = 21,
      colour = "black",
      stroke = 0.25,
      size = 1.6,
      alpha = 0.75) +
    geom_text(
      data = annotation_df,
      aes(
        x = 1.5,
        y = Inf,
        label = annotation_label),
      inherit.aes = FALSE,
      vjust = 1.15,
      size = 7 / .pt) +
    facet_wrap(vars(analysis_scope), nrow = 1) +
    scale_fill_manual(
      values = c(
        "No candidate mutation" = "grey78",
        "Candidate-mutated" = "#BC667E"),
      drop = FALSE,
      name = "Candidate status") +
    scale_x_discrete(
      labels = c(
        "Disease free" = "Disease free",
        "Recurred / progressed" = "Recurred /\nprogressed")) +
    scale_y_log10(
      breaks = scales::breaks_log(n = 5),
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.06, 0.30))) +
    labs(
      x = "",
      y = "Nonsynonymous TMB (mutations/Mb)",
      caption = adjusted_caption) +
    my_theme +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 21, size = 3))) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text = element_text(size = 8),
      legend.position = "bottom")
}

plot_cbio_candidate_mutation_count_recurrence <- function(
    recurrence_data,
    count_tests) {

  plot_df <- recurrence_data %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "DFCI only", "HGSC only")),
      endpoint_group = factor(
        endpoint_group,
        levels = c("Disease free", "Recurred / progressed")))
  annotation_df <- count_tests %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "DFCI only", "HGSC only")),
      annotation_label = paste0(
        "Wilcoxon P ",
        format_clinical_probability(p_value)))

  ggplot(
    plot_df,
    aes(
      x = endpoint_group,
      y = n_candidate_nonsynonymous_mutations,
      fill = endpoint_group)) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.8) +
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.25,
      size = 1.6,
      alpha = 0.75,
      position = position_jitter(
        width = 0.13,
        height = 0,
        seed = 20260728L)) +
    geom_text(
      data = annotation_df,
      aes(
        x = 1.5,
        y = Inf,
        label = annotation_label),
      inherit.aes = FALSE,
      vjust = 1.15,
      size = 7 / .pt) +
    facet_wrap(vars(analysis_scope), nrow = 1) +
    scale_fill_manual(
      values = c(
        "Disease free" = "#4C78A8",
        "Recurred / progressed" = "#D62728"),
      guide = "none") +
    scale_x_discrete(
      labels = c(
        "Disease free" = "Disease free",
        "Recurred / progressed" = "Recurred /\nprogressed")) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 6),
      expand = expansion(mult = c(0.03, 0.28))) +
    labs(
      x = "",
      y = "Nonsynonymous mutations\nin candidate genes",
      caption = paste(
        "Counts include all qualifying mutation records across the 30",
        "candidate genes; zero-count tumours are retained.")) +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text = element_text(size = 8),
      legend.position = "none")
}

plot_cbio_primary_metastatic_all_gene_volcano <- function(association_tbl) {
  plot_df <- association_tbl %>%
    mutate(
      volcano_y = if_else(
        is.na(q_value_all_mutated_genes),
        NA_real_,
        -log10(pmax(
          q_value_all_mutated_genes,
          .Machine$double.xmin))),
      volcano_category = factor(
        volcano_category,
        levels = c(
          "Candidate: primary-associated",
          "Candidate: metastatic-associated",
          "Other FDR-significant gene",
          "Not significant",
          "Not testable")))

  volcano_colours <- c(
    "Candidate: primary-associated" = "#4C78A8",
    "Candidate: metastatic-associated" = "#D62728",
    "Other FDR-significant gene" = "grey20",
    "Not significant" = "grey78",
    "Not testable" = "white")
  volcano_breaks <- names(volcano_colours)
  volcano_breaks <- volcano_breaks[
    volcano_breaks %in% as.character(
      plot_df$volcano_category[!is.na(plot_df$volcano_y)])]
  candidate_shapes <- setNames(
    c(21, 22, 24),
    CBIO_CANDIDATE_GROUP_LEVELS)
  x_limit <- max(
    1,
    ceiling(max(
      abs(plot_df$corrected_log2_odds_ratio),
      na.rm = TRUE)))
  n_tested <- sum(plot_df$test_status == "tested")
  n_significant <- sum(
    plot_df$association_direction %in% c(
      "Primary-associated",
      "Metastatic-associated"))
  n_significant_candidates <- sum(plot_df$is_significant_candidate)

  ggplot(
    plot_df,
    aes(
      x = corrected_log2_odds_ratio,
      y = volcano_y)) +
    geom_hline(
      yintercept = -log10(0.05),
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.3) +
    geom_vline(
      xintercept = 0,
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.3) +
    geom_point(
      data = plot_df %>%
        filter(!is_significant_candidate, !is.na(volcano_y)),
      aes(fill = volcano_category),
      shape = 21,
      colour = "grey35",
      stroke = 0.15,
      size = 1.2,
      alpha = 0.45) +
    geom_point(
      data = plot_df %>% filter(is_significant_candidate),
      aes(
        fill = volcano_category,
        shape = candidate_group_plot),
      colour = "black",
      stroke = 0.4,
      size = 3.2) +
    ggrepel::geom_text_repel(
      data = plot_df %>% filter(is_significant_candidate),
      aes(label = human_symbol),
      size = 7 / .pt,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE) +
    scale_fill_manual(
      values = volcano_colours,
      breaks = volcano_breaks,
      drop = TRUE,
      name = "Genome-wide FDR") +
    scale_shape_manual(
      values = candidate_shapes,
      drop = TRUE,
      name = "Candidate group") +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.14))) +
    labs(
      x = expression(
        "Primary-associated" %<-% log[2]("corrected odds ratio") %->%
          "Metastatic-associated"),
      y = expression(-log[10]("BH Q")),
      caption = paste0(
        format(n_tested, big.mark = ","),
        " testable genes from the ",
        format(nrow(plot_df), big.mark = ","),
        "-gene mutated universe; ",
        format(n_significant, big.mark = ","),
        " genes pass BH Q <= 0.05, including ",
        n_significant_candidates,
        " candidate gene",
        if_else(n_significant_candidates == 1L, "", "s"),
        ".")) +
    my_theme +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          colour = "black",
          alpha = 1,
          size = 3)),
      shape = guide_legend(
        override.aes = list(
          fill = "white",
          colour = "black",
          alpha = 1,
          size = 3))) +
    theme(
      legend.position = "right")
}

plot_cbio_candidate_burden_associations <- function(
    burden_endpoint_tbl,
    burden_association_tbl,
    endpoint_definitions) {

  endpoint_levels <- endpoint_definitions %>%
    arrange(endpoint_order) %>%
    pull(endpoint_label)
  endpoint_group_levels <- unique(c(
    endpoint_definitions %>%
      arrange(endpoint_order) %>%
      pull(reference_group),
    endpoint_definitions %>%
      arrange(endpoint_order) %>%
      pull(case_group)))

  plot_df <- burden_endpoint_tbl %>%
    filter(!is.na(endpoint_group), !is.na(candidate_burden_fraction)) %>%
    mutate(
      endpoint_label = factor(endpoint_label, levels = endpoint_levels),
      endpoint_group = factor(
        endpoint_group,
        levels = endpoint_group_levels),
      endpoint_role = factor(
        endpoint_role,
        levels = c("Reference", "Higher-risk")))

  annotation_df <- burden_association_tbl %>%
    mutate(
      endpoint_label = factor(endpoint_label, levels = endpoint_levels),
      label = paste0(
        "P ",
        format_clinical_probability(p_value),
        "\nQ ",
        format_clinical_probability(q_value_analysis_family)),
      annotation_x = 1.5)

  ggplot(
    plot_df,
    aes(
      x = endpoint_group,
      y = candidate_burden_fraction,
      fill = endpoint_role)) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      linewidth = 0.3) +
    geom_jitter(
      width = 0.12,
      height = 0,
      shape = 21,
      colour = "black",
      stroke = 0.2,
      size = 1.5,
      alpha = 0.75) +
    geom_text(
      data = annotation_df,
      aes(x = annotation_x, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.1,
      size = 7 / .pt) +
    facet_wrap(
      vars(endpoint_label),
      ncol = 3,
      scales = "free_x") +
    scale_fill_manual(
      values = c(
        Reference = "grey75",
        `Higher-risk` = "#BC667E"),
      drop = FALSE) +
    scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      expand = expansion(mult = c(0.05, 0.23))) +
    labs(
      x = "",
      y = "Mutated / profiled candidate genes",
      fill = "") +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
      strip.text = element_text(size = 7),
      legend.position = "bottom")
}

plot_cbio_candidate_tmb_primary_metastatic <- function(
    tmb_data,
    tmb_tests) {

  plot_df <- tmb_data %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "HGSC only")),
      endpoint_group = factor(
        endpoint_group,
        levels = c("Primary", "Metastasis")),
      candidate_status = factor(
        candidate_status,
        levels = c(
          "No candidate mutation",
          "Candidate-mutated")))

  annotation_df <- tmb_tests %>%
    filter(analysis_id %in% c(
      "pooled_candidate_mutated_wilcoxon",
      "hgsc_candidate_mutated_wilcoxon")) %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "HGSC only")),
      annotation_x = 1.5,
      annotation_label = paste0(
        "Candidate-mutated: Primary vs metastatic\nWilcoxon P ",
        format_clinical_probability(p_value)))

  adjusted_result <- tmb_tests %>%
    filter(analysis_id == "study_adjusted_log10_tmb")
  adjusted_caption <- paste0(
    "Study-adjusted metastatic/primary geometric-mean TMB ratio ",
    formatC(adjusted_result$estimate, format = "f", digits = 2),
    " (95% CI ",
    formatC(
      adjusted_result$confidence_interval_low,
      format = "f",
      digits = 2),
    "-",
    formatC(
      adjusted_result$confidence_interval_high,
      format = "f",
      digits = 2),
    "), HC3 P ",
    format_clinical_probability(adjusted_result$p_value),
    ". Pooled boxes are descriptive; HGSC is the only cohort with both disease groups.")

  ggplot(
    plot_df,
    aes(
      x = endpoint_group,
      y = tmb_nonsynonymous,
      fill = candidate_status)) +
    geom_boxplot(
      position = position_dodge(width = 0.72),
      width = 0.62,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.8) +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0.12,
        jitter.height = 0,
        dodge.width = 0.72,
        seed = 20260728L),
      shape = 21,
      colour = "black",
      stroke = 0.25,
      size = 1.7,
      alpha = 0.75) +
    geom_text(
      data = annotation_df,
      aes(
        x = annotation_x,
        y = Inf,
        label = annotation_label),
      inherit.aes = FALSE,
      vjust = 1.15,
      size = 7 / .pt) +
    facet_wrap(
      vars(analysis_scope),
      nrow = 1) +
    scale_fill_manual(
      values = c(
        "No candidate mutation" = "grey78",
        "Candidate-mutated" = "#BC667E"),
      drop = FALSE,
      name = "Candidate status") +
    scale_y_log10(
      breaks = scales::breaks_log(n = 5),
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.06, 0.28))) +
    labs(
      x = "",
      y = "Nonsynonymous TMB (mutations/Mb)",
      caption = adjusted_caption) +
    my_theme +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(shape = 21, size = 3))) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text = element_text(size = 8),
      legend.position = "bottom")
}

plot_cbio_candidate_mutation_count_primary_metastatic <- function(
    tmb_data,
    count_tests) {

  plot_df <- tmb_data %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "HGSC only")),
      endpoint_group = factor(
        endpoint_group,
        levels = c("Primary", "Metastasis")))

  annotation_df <- count_tests %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "HGSC only")),
      annotation_label = paste0(
        "Primary vs metastatic\nWilcoxon P ",
        format_clinical_probability(p_value)))

  ggplot(
    plot_df,
    aes(
      x = endpoint_group,
      y = n_candidate_nonsynonymous_mutations,
      fill = endpoint_group)) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.8) +
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.25,
      size = 1.7,
      alpha = 0.75,
      position = position_jitter(
        width = 0.13,
        height = 0,
        seed = 20260728L)) +
    geom_text(
      data = annotation_df,
      aes(
        x = 1.5,
        y = Inf,
        label = annotation_label),
      inherit.aes = FALSE,
      vjust = 1.15,
      size = 7 / .pt) +
    facet_wrap(vars(analysis_scope), nrow = 1) +
    scale_fill_manual(
      values = c(
        Primary = "#4C78A8",
        Metastasis = "#D62728"),
      guide = "none") +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 6),
      expand = expansion(mult = c(0.03, 0.28))) +
    labs(
      x = "",
      y = "Nonsynonymous mutations\nin candidate genes",
      caption = paste(
        "Counts include all qualifying mutation records across the 30",
        "candidate genes; zero-count tumours are retained.",
        "The pooled comparison is descriptive because disease status is",
        "strongly confounded with study.")) +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text = element_text(size = 8),
      legend.position = "none")
}

validate_cbio_candidate_tmb_analysis <- function(
    tmb_data,
    tmb_summary,
    tmb_tests) {

  scope_counts <- tmb_data %>%
    distinct(analysis_scope, study_id, sample_id) %>%
    count(analysis_scope, name = "n_samples") %>%
    mutate(analysis_scope = as.character(analysis_scope))
  expected_scope_counts <- tibble(
    analysis_scope = c("Pooled", "HGSC only"),
    n_samples = c(176L, 39L))

  if (!identical(
      arrange(scope_counts, analysis_scope),
      arrange(expected_scope_counts, analysis_scope))) {
    stop("Unexpected sample counts in the candidate-associated TMB analysis.")
  }

  duplicate_samples <- tmb_data %>%
    count(analysis_scope, study_id, sample_id) %>%
    filter(n != 1L)

  if (nrow(duplicate_samples) > 0L) {
    stop("Candidate-associated TMB data contain duplicated samples.")
  }

  if (any(
      is.na(tmb_data$tmb_nonsynonymous) |
      !is.finite(tmb_data$tmb_nonsynonymous) |
      tmb_data$tmb_nonsynonymous <= 0)) {
    stop("Candidate-associated TMB values must be complete, finite and positive.")
  }

  observed_group_counts <- tmb_data %>%
    count(
      analysis_scope,
      endpoint_group,
      candidate_status,
      name = "n_samples") %>%
    mutate(
      across(
        c(analysis_scope, endpoint_group, candidate_status),
        as.character)) %>%
    arrange(analysis_scope, endpoint_group, candidate_status)
  expected_group_counts <- tribble(
    ~analysis_scope, ~endpoint_group, ~candidate_status, ~n_samples,
    "Pooled", "Primary", "Candidate-mutated", 94L,
    "Pooled", "Primary", "No candidate mutation", 21L,
    "Pooled", "Metastasis", "Candidate-mutated", 50L,
    "Pooled", "Metastasis", "No candidate mutation", 11L,
    "HGSC only", "Primary", "Candidate-mutated", 29L,
    "HGSC only", "Primary", "No candidate mutation", 3L,
    "HGSC only", "Metastasis", "Candidate-mutated", 7L) %>%
    arrange(analysis_scope, endpoint_group, candidate_status)

  if (!identical(observed_group_counts, expected_group_counts)) {
    stop("Unexpected disease/candidate-status counts in the TMB analysis.")
  }

  if (nrow(tmb_summary) != nrow(expected_group_counts) ||
      any(tmb_summary$n_samples <= 0L)) {
    stop("The candidate-associated TMB summary is incomplete.")
  }

  if (!all(str_detect(
      tmb_data$candidate_mutation_definition,
      "exact variants are not matched"))) {
    stop("Candidate TMB mutation definitions do not exclude exact-variant matching.")
  }

  expected_test_ids <- c(
    "study_adjusted_log10_tmb",
    "pooled_candidate_mutated_wilcoxon",
    "hgsc_candidate_mutated_wilcoxon")

  if (nrow(tmb_tests) != 3L ||
      !setequal(tmb_tests$analysis_id, expected_test_ids) ||
      any(tmb_tests$test_status != "tested")) {
    stop("The candidate-associated TMB statistical tests are incomplete.")
  }

  adjusted_result <- tmb_tests %>%
    filter(analysis_id == "study_adjusted_log10_tmb")
  pooled_result <- tmb_tests %>%
    filter(analysis_id == "pooled_candidate_mutated_wilcoxon")
  hgsc_result <- tmb_tests %>%
    filter(analysis_id == "hgsc_candidate_mutated_wilcoxon")

  if (adjusted_result$n_reference != 94L ||
      adjusted_result$n_case != 50L ||
      adjusted_result$model_n != 144L ||
      !isTRUE(adjusted_result$model_is_full_rank) ||
      adjusted_result$model_design_rank !=
        adjusted_result$model_design_columns) {
    stop("The study-adjusted candidate TMB model failed validation.")
  }

  if (pooled_result$n_reference != 94L ||
      pooled_result$n_case != 50L ||
      hgsc_result$n_reference != 29L ||
      hgsc_result$n_case != 7L) {
    stop("Candidate TMB Wilcoxon group sizes are unexpected.")
  }

  probability_values <- tmb_tests$p_value[!is.na(tmb_tests$p_value)]
  if (any(probability_values < 0 | probability_values > 1)) {
    stop("Candidate-associated TMB P values fall outside [0, 1].")
  }

  if (adjusted_result$confidence_interval_low <= 0 ||
      adjusted_result$confidence_interval_high <=
        adjusted_result$confidence_interval_low ||
      adjusted_result$estimate < adjusted_result$confidence_interval_low ||
      adjusted_result$estimate > adjusted_result$confidence_interval_high) {
    stop("The adjusted TMB ratio or confidence interval is invalid.")
  }

  invisible(TRUE)
}

validate_cbio_candidate_mutation_count_analysis <- function(
    tmb_data,
    count_summary,
    count_tests) {

  if (any(
      is.na(tmb_data$n_candidate_nonsynonymous_mutations) |
      tmb_data$n_candidate_nonsynonymous_mutations < 0 |
      tmb_data$n_candidate_nonsynonymous_mutations !=
        as.integer(tmb_data$n_candidate_nonsynonymous_mutations))) {
    stop("Candidate mutation counts must be complete non-negative integers.")
  }

  if (any(
      (tmb_data$n_candidate_nonsynonymous_mutations > 0L) !=
        (tmb_data$candidate_status == "Candidate-mutated"))) {
    stop("Candidate mutation counts disagree with candidate mutation status.")
  }

  expected_summary_counts <- tribble(
    ~analysis_scope, ~endpoint_group, ~n_samples,
    "Pooled", "Primary", 115L,
    "Pooled", "Metastasis", 61L,
    "HGSC only", "Primary", 32L,
    "HGSC only", "Metastasis", 7L) %>%
    mutate(
      analysis_scope = factor(
        analysis_scope,
        levels = c("Pooled", "HGSC only")),
      endpoint_group = factor(
        endpoint_group,
        levels = c("Primary", "Metastasis"))) %>%
    arrange(analysis_scope, endpoint_group)
  observed_summary_counts <- count_summary %>%
    select(analysis_scope, endpoint_group, n_samples) %>%
    arrange(analysis_scope, endpoint_group)

  if (!identical(observed_summary_counts, expected_summary_counts)) {
    stop("Candidate mutation-count sample sizes are unexpected.")
  }

  expected_test_sizes <- tribble(
    ~analysis_id, ~n_reference, ~n_case,
    "pooled_candidate_mutation_count_wilcoxon", 115L, 61L,
    "hgsc_candidate_mutation_count_wilcoxon", 32L, 7L)
  observed_test_sizes <- count_tests %>%
    select(analysis_id, n_reference, n_case)

  if (!identical(observed_test_sizes, expected_test_sizes) ||
      any(count_tests$test_status != "tested") ||
      any(count_tests$p_value < 0 | count_tests$p_value > 1)) {
    stop("Candidate mutation-count tests failed validation.")
  }

  invisible(TRUE)
}

validate_cbio_clinical_analysis <- function(
    endpoint_rows,
    mutation_matrix,
    association_tbl,
    burden_association_tbl,
    focused_association_tbl) {

  endpoint_counts <- endpoint_rows %>%
    filter(!is.na(endpoint_group)) %>%
    distinct(endpoint_id, study_id, sample_id, endpoint_group) %>%
    count(endpoint_id, endpoint_group, name = "n_samples")

  assert_endpoint_count <- function(endpoint_id, endpoint_group, expected_n) {
    observed_n <- endpoint_counts %>%
      filter(
        .data$endpoint_id == .env$endpoint_id,
        .data$endpoint_group == .env$endpoint_group) %>%
      pull(n_samples)

    if (length(observed_n) != 1L || observed_n != expected_n) {
      stop(
        "Expected ",
        expected_n,
        " ",
        endpoint_group,
        " samples for ",
        endpoint_id,
        "; observed ",
        if_else(length(observed_n) == 0L, "0", paste(observed_n, collapse = ",")),
        ".")
    }
  }

  assert_endpoint_count("sample_type_pooled", "Primary", 115L)
  assert_endpoint_count("sample_type_pooled", "Metastasis", 61L)
  assert_endpoint_count("sample_type_explicit", "Primary", 32L)
  assert_endpoint_count("sample_type_explicit", "Metastasis", 32L)
  assert_endpoint_count("sample_type_hgsc", "Primary", 32L)
  assert_endpoint_count("sample_type_hgsc", "Metastasis", 7L)
  assert_endpoint_count("recurrence_pooled", "Disease free", 41L)
  assert_endpoint_count(
    "recurrence_pooled",
    "Recurred / progressed",
    25L)

  pooled_study_counts <- endpoint_rows %>%
    filter(
      endpoint_id == "sample_type_pooled",
      !is.na(endpoint_group)) %>%
    distinct(study_id, sample_id, endpoint_group) %>%
    count(study_id, endpoint_group, name = "n_samples")
  expected_pooled_study_counts <- tribble(
    ~study_id, ~endpoint_group, ~n_samples,
    "cscc_dfarber_2015", "Metastasis", 29L,
    "cscc_hgsc_bcm_2014", "Primary", 32L,
    "cscc_hgsc_bcm_2014", "Metastasis", 7L,
    "cscc_ranson_2022", "Metastasis", 25L,
    "cscc_ucsf_2021", "Primary", 83L)

  if (!identical(
      arrange(pooled_study_counts, study_id, endpoint_group),
      arrange(expected_pooled_study_counts, study_id, endpoint_group))) {
    stop("The all-cohort primary/metastatic study assignments are unexpected.")
  }

  recurrence_study_counts <- endpoint_rows %>%
    filter(
      endpoint_id == "recurrence_pooled",
      !is.na(endpoint_group)) %>%
    distinct(study_id, sample_id, endpoint_group) %>%
    count(study_id, endpoint_group, name = "n_samples")
  expected_recurrence_study_counts <- tribble(
    ~study_id, ~endpoint_group, ~n_samples,
    "cscc_dfarber_2015", "Disease free", 18L,
    "cscc_dfarber_2015", "Recurred / progressed", 11L,
    "cscc_hgsc_bcm_2014", "Disease free", 23L,
    "cscc_hgsc_bcm_2014", "Recurred / progressed", 14L)

  if (!identical(
      arrange(recurrence_study_counts, study_id, endpoint_group),
      arrange(expected_recurrence_study_counts, study_id, endpoint_group))) {
    stop("The pooled recurrence study assignments are unexpected.")
  }

  expected_candidate_genes <- sort(unique(candidate_genes$gene_symbol))
  observed_candidate_genes <- sort(unique(mutation_matrix$gene_symbol))

  if (!identical(expected_candidate_genes, observed_candidate_genes) ||
      length(observed_candidate_genes) != 30L) {
    stop("The clinical mutation matrix does not contain the 30 final candidates.")
  }

  if (any(mutation_matrix$gene_symbol %in% CBIO_EXTRA_ONCOPLOT_GENES)) {
    stop("Supplemental cSCC driver genes entered the clinical analysis.")
  }

  stage_endpoint_ids <- c(
    "hgsc_t_stage",
    "hgsc_n_stage",
    "hgsc_m_stage",
    "ranson_nodal_stage")
  stage_groups <- endpoint_rows %>%
    filter(endpoint_id %in% stage_endpoint_ids, !is.na(endpoint_group)) %>%
    distinct(endpoint_id, endpoint_group) %>%
    count(endpoint_id, name = "n_groups")

  if (nrow(stage_groups) != length(stage_endpoint_ids) ||
      any(stage_groups$n_groups != 2L)) {
    stop("One or more normalized staging endpoints lacks two analysis groups.")
  }

  probability_columns <- c(
    association_tbl$p_value,
    association_tbl$q_value_within_endpoint,
    association_tbl$q_value_analysis_family,
    burden_association_tbl$p_value,
    burden_association_tbl$q_value_analysis_family)
  probability_columns <- probability_columns[!is.na(probability_columns)]

  if (any(probability_columns < 0 | probability_columns > 1)) {
    stop("Clinical association P or Q values fall outside [0, 1].")
  }

  expected_focused_rows <- length(CBIO_FOCUSED_PRIMARY_METASTATIC_GENES) *
    nrow(CBIO_FOCUSED_SENSITIVITY_DEFINITIONS)
  if (nrow(focused_association_tbl) != expected_focused_rows ||
      !setequal(
        focused_association_tbl$gene_symbol,
        CBIO_FOCUSED_PRIMARY_METASTATIC_GENES)) {
    stop("The focused KNDC1/FOXP4 sensitivity analysis is incomplete.")
  }

  invisible(TRUE)
}

validate_cbio_candidate_recurrence_analysis <- function(
    adjusted_associations,
    candidate_tbl) {

  expected_genes <- sort(unique(candidate_tbl$gene_symbol))
  observed_genes <- sort(unique(adjusted_associations$gene_symbol))

  if (nrow(adjusted_associations) != 30L ||
      !identical(expected_genes, observed_genes)) {
    stop("The adjusted recurrence analysis does not contain 30 candidates.")
  }

  if (any(
      adjusted_associations$n_reference_profiled > 41L |
      adjusted_associations$n_case_profiled > 25L |
      adjusted_associations$n_reference_mutated >
        adjusted_associations$n_reference_profiled |
      adjusted_associations$n_case_mutated >
        adjusted_associations$n_case_profiled)) {
    stop("Candidate recurrence counts exceed endpoint denominators.")
  }

  probability_values <- c(
    adjusted_associations$p_value,
    adjusted_associations$q_value)
  probability_values <- probability_values[!is.na(probability_values)]

  if (any(probability_values < 0 | probability_values > 1)) {
    stop("Candidate recurrence P or Q values fall outside [0, 1].")
  }

  tested_rows <- adjusted_associations %>%
    filter(test_status == "tested")
  if (nrow(tested_rows) == 0L ||
      any(
        tested_rows$common_odds_ratio <= 0 |
        tested_rows$common_odds_ratio_conf_low <= 0 |
        tested_rows$common_odds_ratio_conf_high <=
          tested_rows$common_odds_ratio_conf_low |
        tested_rows$common_odds_ratio <
          tested_rows$common_odds_ratio_conf_low |
        tested_rows$common_odds_ratio >
          tested_rows$common_odds_ratio_conf_high)) {
    stop("Study-adjusted recurrence odds ratios or intervals are invalid.")
  }

  if (!all(str_detect(
      adjusted_associations$mutation_definition,
      "exact mouse variants are not matched"))) {
    stop("Recurrence analysis unexpectedly uses exact-variant matching.")
  }

  invisible(TRUE)
}

validate_cbio_candidate_recurrence_boxplots <- function(
    recurrence_data,
    tmb_summary,
    tmb_tests,
    count_summary,
    count_tests) {

  scope_counts <- recurrence_data %>%
    distinct(analysis_scope, study_id, sample_id) %>%
    count(analysis_scope, name = "n_samples") %>%
    mutate(analysis_scope = as.character(analysis_scope))
  expected_scope_counts <- tibble(
    analysis_scope = c("Pooled", "DFCI only", "HGSC only"),
    n_samples = c(66L, 29L, 37L))
  if (!identical(
      arrange(scope_counts, analysis_scope),
      arrange(expected_scope_counts, analysis_scope))) {
    stop("Unexpected sample counts in the recurrence boxplot data.")
  }

  duplicate_samples <- recurrence_data %>%
    count(analysis_scope, study_id, sample_id) %>%
    filter(n != 1L)
  if (nrow(duplicate_samples) > 0L) {
    stop("Recurrence boxplot data contain duplicated samples.")
  }
  if (any(
      is.na(recurrence_data$tmb_nonsynonymous) |
      !is.finite(recurrence_data$tmb_nonsynonymous) |
      recurrence_data$tmb_nonsynonymous <= 0 |
      is.na(recurrence_data$n_candidate_nonsynonymous_mutations) |
      recurrence_data$n_candidate_nonsynonymous_mutations < 0)) {
    stop("Recurrence TMB or candidate mutation counts are invalid.")
  }

  observed_group_counts <- recurrence_data %>%
    count(
      analysis_scope,
      endpoint_group,
      candidate_status,
      name = "n_samples") %>%
    mutate(
      across(
        c(analysis_scope, endpoint_group, candidate_status),
        as.character)) %>%
    arrange(analysis_scope, endpoint_group, candidate_status)
  expected_group_counts <- tribble(
    ~analysis_scope, ~endpoint_group, ~candidate_status, ~n_samples,
    "Pooled", "Disease free", "Candidate-mutated", 31L,
    "Pooled", "Disease free", "No candidate mutation", 10L,
    "Pooled", "Recurred / progressed", "Candidate-mutated", 21L,
    "Pooled", "Recurred / progressed", "No candidate mutation", 4L,
    "DFCI only", "Disease free", "Candidate-mutated", 11L,
    "DFCI only", "Disease free", "No candidate mutation", 7L,
    "DFCI only", "Recurred / progressed", "Candidate-mutated", 7L,
    "DFCI only", "Recurred / progressed", "No candidate mutation", 4L,
    "HGSC only", "Disease free", "Candidate-mutated", 20L,
    "HGSC only", "Disease free", "No candidate mutation", 3L,
    "HGSC only", "Recurred / progressed", "Candidate-mutated", 14L) %>%
    arrange(analysis_scope, endpoint_group, candidate_status)
  if (!identical(observed_group_counts, expected_group_counts) ||
      nrow(tmb_summary) != nrow(expected_group_counts)) {
    stop("Unexpected recurrence disease/candidate-status counts.")
  }

  adjusted_result <- tmb_tests %>%
    filter(analysis_id == "study_adjusted_log10_tmb_recurrence")
  expected_wilcoxon_sizes <- tribble(
    ~analysis_id, ~n_reference, ~n_case,
    "pooled_candidate_mutated_recurrence_wilcoxon", 31L, 21L,
    "dfci_candidate_mutated_recurrence_wilcoxon", 11L, 7L,
    "hgsc_candidate_mutated_recurrence_wilcoxon", 20L, 14L)
  observed_wilcoxon_sizes <- tmb_tests %>%
    filter(analysis_id != "study_adjusted_log10_tmb_recurrence") %>%
    select(analysis_id, n_reference, n_case)
  if (nrow(tmb_tests) != 4L ||
      nrow(adjusted_result) != 1L ||
      adjusted_result$n_reference != 31L ||
      adjusted_result$n_case != 21L ||
      adjusted_result$model_n != 52L ||
      !isTRUE(adjusted_result$model_is_full_rank) ||
      !identical(observed_wilcoxon_sizes, expected_wilcoxon_sizes)) {
    stop("Recurrence TMB tests or model dimensions are unexpected.")
  }

  expected_count_sizes <- tribble(
    ~analysis_id, ~n_reference, ~n_case,
    "pooled_recurrence_candidate_mutation_count_wilcoxon", 41L, 25L,
    "dfci_recurrence_candidate_mutation_count_wilcoxon", 18L, 11L,
    "hgsc_recurrence_candidate_mutation_count_wilcoxon", 23L, 14L)
  observed_count_sizes <- count_tests %>%
    select(analysis_id, n_reference, n_case)
  if (nrow(count_summary) != 6L ||
      !identical(observed_count_sizes, expected_count_sizes)) {
    stop("Recurrence candidate mutation-count outputs are incomplete.")
  }

  probability_values <- c(tmb_tests$p_value, count_tests$p_value)
  probability_values <- probability_values[!is.na(probability_values)]
  if (any(probability_values < 0 | probability_values > 1) ||
      adjusted_result$confidence_interval_low <= 0 ||
      adjusted_result$confidence_interval_high <=
        adjusted_result$confidence_interval_low) {
    stop("Recurrence boxplot test probabilities or confidence interval are invalid.")
  }

  invisible(TRUE)
}

validate_cbio_all_gene_primary_metastatic_analysis <- function(
    association_tbl,
    incidence_tbl,
    candidate_tbl) {

  expected_gene_symbols <- sort(unique(incidence_tbl$human_symbol))
  observed_gene_symbols <- sort(unique(association_tbl$human_symbol))

  if (length(expected_gene_symbols) != 18517L ||
      nrow(association_tbl) != 18517L ||
      !identical(expected_gene_symbols, observed_gene_symbols)) {
    stop(
      "The all-gene primary/metastatic analysis does not contain the ",
      "18,517-gene cBioPortal universe.")
  }

  if (any(!association_tbl$profile_denominator_matches_incidence)) {
    mismatched_genes <- association_tbl %>%
      filter(!profile_denominator_matches_incidence) %>%
      pull(human_symbol)
    stop(
      "Panel-aware denominators do not match cBioPortal_Mutated_Genes.txt for: ",
      paste(head(mismatched_genes, 20), collapse = ", "),
      if_else(length(mismatched_genes) > 20L, ", ...", ""),
      ".")
  }

  if (!setequal(
      unique(association_tbl$n_reference_profiled),
      115L) ||
      !setequal(
        unique(association_tbl$n_case_profiled),
        c(32L, 61L))) {
    stop("Unexpected primary/metastatic genome-wide profiling denominators.")
  }

  if (any(
      association_tbl$n_reference_mutated >
        association_tbl$n_reference_profiled |
      association_tbl$n_case_mutated >
        association_tbl$n_case_profiled)) {
    stop("Genome-wide mutated counts exceed profiled denominators.")
  }

  probability_values <- c(
    association_tbl$p_value,
    association_tbl$q_value_all_mutated_genes)
  probability_values <- probability_values[!is.na(probability_values)]

  if (any(probability_values < 0 | probability_values > 1)) {
    stop("Genome-wide primary/metastatic P or Q values fall outside [0, 1].")
  }

  expected_candidate_human_symbols <- build_cbio_candidate_human_metadata(
    candidate_tbl) %>%
    semi_join(incidence_tbl, by = "human_symbol") %>%
    pull(human_symbol) %>%
    sort()
  observed_candidate_human_symbols <- association_tbl %>%
    filter(is_candidate) %>%
    pull(human_symbol) %>%
    sort()

  if (!identical(
      expected_candidate_human_symbols,
      observed_candidate_human_symbols)) {
    stop("Candidate genes were not mapped correctly onto the all-gene analysis.")
  }

  if (any(
      association_tbl$is_significant_candidate &
        (
          is.na(association_tbl$q_value_all_mutated_genes) |
          association_tbl$q_value_all_mutated_genes > 0.05))) {
    stop("A candidate was highlighted without genome-wide FDR significance.")
  }

  invisible(TRUE)
}

cbio_clinical_long <- load_cbio_clinical_data(CBIO_STUDY_IDS)
cbio_sample_clinical <- build_cbio_sample_clinical_data(cbio_clinical_long)
cbio_clinical_endpoint_rows <- build_cbio_clinical_endpoint_rows(
  cbio_sample_clinical)
cbio_clinical_endpoint_availability <- summarise_cbio_endpoint_availability(
  cbio_clinical_endpoint_rows)

cbio_all_mutated_human_symbols <- cbio_mutated_gene_incidence %>%
  pull(human_symbol) %>%
  unique() %>%
  sort()
cbio_all_mutated_gene_map <- tibble(
  gene_symbol = cbio_all_mutated_human_symbols,
  human_symbol = cbio_all_mutated_human_symbols)
cbio_all_mutated_panel_gene_members <- load_cbio_panel_gene_members(
  panel_ids = cbio_sample_gene_panels$gene_panel,
  query_human_symbols = cbio_all_mutated_human_symbols)
cbio_all_mutated_profiled_sample_map <- build_cbio_profiled_sample_map(
  mouse_gene_map = cbio_all_mutated_gene_map,
  sample_gene_panels = cbio_sample_gene_panels,
  panel_gene_members = cbio_all_mutated_panel_gene_members)
cbio_all_study_mutations <- load_cbio_all_study_mutations(
  profiled_samples = cbio_profiled_samples,
  query_human_symbols = cbio_all_mutated_human_symbols)
cbio_all_gene_primary_metastatic_associations <-
  run_cbio_all_gene_primary_metastatic_associations(
    incidence_tbl = cbio_mutated_gene_incidence,
    profiled_map = cbio_all_mutated_profiled_sample_map,
    cbio_mutations = cbio_all_study_mutations,
    endpoint_rows = cbio_clinical_endpoint_rows,
    candidate_tbl = candidate_genes)

cbio_candidate_clinical_matrix <- prepare_cbio_candidate_clinical_matrix(
  oncoplot_df = cbio_candidate_oncoplot_df,
  candidate_tbl = candidate_genes)
cbio_candidate_clinical_associations <- run_cbio_candidate_clinical_associations(
  mutation_matrix = cbio_candidate_clinical_matrix,
  endpoint_rows = cbio_clinical_endpoint_rows)
cbio_candidate_recurrence_unadjusted_associations <-
  cbio_candidate_clinical_associations %>%
  filter(endpoint_id %in% c(
    "recurrence_pooled",
    "dfci_dfs_status",
    "hgsc_recurrence_persistence"))
cbio_candidate_recurrence_adjusted_associations <-
  run_cbio_candidate_recurrence_adjusted_associations(
    mutation_matrix = cbio_candidate_clinical_matrix,
    endpoint_rows = cbio_clinical_endpoint_rows)
cbio_primary_metastatic_prevalence <- summarise_cbio_primary_metastatic_prevalence(
  mutation_matrix = cbio_candidate_clinical_matrix,
  endpoint_rows = cbio_clinical_endpoint_rows)
cbio_primary_metastatic_highlight_data <-
  prepare_cbio_primary_metastatic_highlight_data(
    cbio_candidate_clinical_associations)
cbio_focused_sensitivity_endpoint_rows <-
  build_cbio_focused_sensitivity_endpoint_rows(cbio_clinical_endpoint_rows)
cbio_focused_primary_metastatic_associations <-
  run_cbio_focused_primary_metastatic_associations(
    mutation_matrix = cbio_candidate_clinical_matrix,
    focused_endpoint_rows = cbio_focused_sensitivity_endpoint_rows)
cbio_focused_all_clinical_associations <-
  cbio_candidate_clinical_associations %>%
  filter(gene_symbol %in% CBIO_FOCUSED_PRIMARY_METASTATIC_GENES)
cbio_focused_primary_metastatic_prevalence <-
  candidate_genes %>%
  filter(gene_symbol %in% CBIO_FOCUSED_PRIMARY_METASTATIC_GENES) %>%
  transmute(
    gene_symbol,
    gene_label = as.character(gene_label),
    candidate_group = as.character(candidate_group),
    human_symbols) %>%
  distinct() %>%
  crossing(CBIO_PRIMARY_METASTATIC_COHORT_GROUPS) %>%
  left_join(
    cbio_primary_metastatic_prevalence %>%
      filter(gene_symbol %in% CBIO_FOCUSED_PRIMARY_METASTATIC_GENES) %>%
      mutate(
        comparison_cohort = as.character(comparison_cohort),
        endpoint_group = as.character(endpoint_group)) %>%
      select(
        gene_symbol,
        comparison_cohort,
        endpoint_group,
        n_profiled,
        n_mutated,
        prevalence_pct),
    by = c(
      "gene_symbol",
      "comparison_cohort",
      "endpoint_group")) %>%
  mutate(
    n_profiled = replace_na(n_profiled, 0L),
    n_mutated = replace_na(n_mutated, 0L),
    comparison_cohort = factor(
      comparison_cohort,
      levels = c("DFCI", "HGSC", "Ranson", "UCSF", "Pooled")),
    endpoint_group = factor(
      endpoint_group,
      levels = c("Primary", "Metastasis"))) %>%
  arrange(gene_symbol, comparison_cohort, endpoint_group)

cbio_candidate_burden <- build_cbio_candidate_burden(
  mutation_matrix = cbio_candidate_clinical_matrix,
  sample_clinical = cbio_sample_clinical)
cbio_candidate_burden_endpoint_data <- cbio_candidate_burden %>%
  inner_join(
    cbio_clinical_endpoint_rows,
    by = c("study_id", "sample_id"),
    relationship = "many-to-many")
cbio_candidate_burden_associations <- run_cbio_candidate_burden_associations(
  burden_tbl = cbio_candidate_burden,
  endpoint_rows = cbio_clinical_endpoint_rows)
cbio_candidate_tmb_primary_metastatic_data <-
  prepare_cbio_candidate_tmb_data(
    burden_tbl = cbio_candidate_burden,
    endpoint_rows = cbio_clinical_endpoint_rows)
cbio_candidate_tmb_primary_metastatic_summary <-
  summarise_cbio_candidate_tmb(
    cbio_candidate_tmb_primary_metastatic_data)
cbio_candidate_tmb_primary_metastatic_tests <-
  run_cbio_candidate_tmb_tests(
    cbio_candidate_tmb_primary_metastatic_data)
cbio_candidate_mutation_count_primary_metastatic_summary <-
  summarise_cbio_candidate_mutation_counts(
    cbio_candidate_tmb_primary_metastatic_data)
cbio_candidate_mutation_count_primary_metastatic_tests <-
  run_cbio_candidate_mutation_count_tests(
    cbio_candidate_tmb_primary_metastatic_data)
cbio_candidate_gene_count_tmb_correlation_data <-
  prepare_cbio_candidate_gene_count_tmb_correlation_data(
    cbio_candidate_tmb_primary_metastatic_data)
cbio_candidate_gene_count_tmb_correlation_test <-
  run_cbio_candidate_gene_count_tmb_spearman_test(
    cbio_candidate_gene_count_tmb_correlation_data)
cbio_candidate_recurrence_boxplot_data <-
  prepare_cbio_candidate_recurrence_boxplot_data(
    burden_tbl = cbio_candidate_burden,
    endpoint_rows = cbio_clinical_endpoint_rows)
cbio_candidate_tmb_recurrence_summary <-
  summarise_cbio_candidate_recurrence_tmb(
    cbio_candidate_recurrence_boxplot_data)
cbio_candidate_tmb_recurrence_tests <-
  run_cbio_candidate_recurrence_tmb_tests(
    cbio_candidate_recurrence_boxplot_data)
cbio_candidate_mutation_count_recurrence_summary <-
  summarise_cbio_candidate_recurrence_mutation_counts(
    cbio_candidate_recurrence_boxplot_data)
cbio_candidate_mutation_count_recurrence_tests <-
  run_cbio_candidate_recurrence_mutation_count_tests(
    cbio_candidate_recurrence_boxplot_data)

validate_cbio_clinical_analysis(
  endpoint_rows = cbio_clinical_endpoint_rows,
  mutation_matrix = cbio_candidate_clinical_matrix,
  association_tbl = cbio_candidate_clinical_associations,
  burden_association_tbl = cbio_candidate_burden_associations,
  focused_association_tbl = cbio_focused_primary_metastatic_associations)
validate_cbio_all_gene_primary_metastatic_analysis(
  association_tbl = cbio_all_gene_primary_metastatic_associations,
  incidence_tbl = cbio_mutated_gene_incidence,
  candidate_tbl = candidate_genes)
validate_cbio_candidate_tmb_analysis(
  tmb_data = cbio_candidate_tmb_primary_metastatic_data,
  tmb_summary = cbio_candidate_tmb_primary_metastatic_summary,
  tmb_tests = cbio_candidate_tmb_primary_metastatic_tests)
validate_cbio_candidate_mutation_count_analysis(
  tmb_data = cbio_candidate_tmb_primary_metastatic_data,
  count_summary = cbio_candidate_mutation_count_primary_metastatic_summary,
  count_tests = cbio_candidate_mutation_count_primary_metastatic_tests)
validate_cbio_candidate_gene_count_tmb_correlation(
  correlation_data = cbio_candidate_gene_count_tmb_correlation_data,
  correlation_test = cbio_candidate_gene_count_tmb_correlation_test)
validate_cbio_candidate_recurrence_analysis(
  adjusted_associations =
    cbio_candidate_recurrence_adjusted_associations,
  candidate_tbl = candidate_genes)
validate_cbio_candidate_recurrence_boxplots(
  recurrence_data = cbio_candidate_recurrence_boxplot_data,
  tmb_summary = cbio_candidate_tmb_recurrence_summary,
  tmb_tests = cbio_candidate_tmb_recurrence_tests,
  count_summary = cbio_candidate_mutation_count_recurrence_summary,
  count_tests = cbio_candidate_mutation_count_recurrence_tests)

clinical_gene_levels <- levels(skin_gene_vaf$gene_label)
if (is.null(clinical_gene_levels)) {
  clinical_gene_levels <- candidate_genes %>%
    arrange(candidate_group, gene_label) %>%
    pull(gene_label) %>%
    as.character()
}

cbio_candidate_clinical_association_plot <- plot_cbio_candidate_clinical_associations(
  association_tbl = cbio_candidate_clinical_associations,
  gene_levels = clinical_gene_levels,
  endpoint_definitions = CBIO_CLINICAL_ENDPOINT_DEFINITIONS)
cbio_primary_metastatic_prevalence_plot <- plot_cbio_primary_metastatic_prevalence(
  prevalence_tbl = cbio_primary_metastatic_prevalence,
  gene_levels = clinical_gene_levels)
cbio_primary_metastatic_candidate_volcano_plot <-
  plot_cbio_primary_metastatic_candidate_volcano(
    cbio_primary_metastatic_highlight_data)
cbio_candidate_recurrence_volcano_plot <-
  plot_cbio_candidate_recurrence_volcano(
    cbio_candidate_recurrence_adjusted_associations)
cbio_candidate_tmb_recurrence_plot <-
  plot_cbio_candidate_tmb_recurrence(
    recurrence_data = cbio_candidate_recurrence_boxplot_data,
    recurrence_tests = cbio_candidate_tmb_recurrence_tests)
cbio_candidate_mutation_count_recurrence_plot <-
  plot_cbio_candidate_mutation_count_recurrence(
    recurrence_data = cbio_candidate_recurrence_boxplot_data,
    count_tests = cbio_candidate_mutation_count_recurrence_tests)
cbio_primary_metastatic_all_gene_volcano_plot <-
  plot_cbio_primary_metastatic_all_gene_volcano(
    cbio_all_gene_primary_metastatic_associations)
cbio_focused_primary_metastatic_prevalence_plot <-
  plot_cbio_focused_primary_metastatic_prevalence(
    prevalence_tbl = cbio_focused_primary_metastatic_prevalence)
cbio_focused_primary_metastatic_sensitivity_plot <-
  plot_cbio_focused_primary_metastatic_sensitivity(
    focused_results = cbio_focused_primary_metastatic_associations)
cbio_candidate_burden_association_plot <- plot_cbio_candidate_burden_associations(
  burden_endpoint_tbl = cbio_candidate_burden_endpoint_data,
  burden_association_tbl = cbio_candidate_burden_associations,
  endpoint_definitions = CBIO_CLINICAL_ENDPOINT_DEFINITIONS)
cbio_candidate_tmb_primary_metastatic_plot <-
  plot_cbio_candidate_tmb_primary_metastatic(
    tmb_data = cbio_candidate_tmb_primary_metastatic_data,
    tmb_tests = cbio_candidate_tmb_primary_metastatic_tests)
cbio_candidate_mutation_count_primary_metastatic_plot <-
  plot_cbio_candidate_mutation_count_primary_metastatic(
    tmb_data = cbio_candidate_tmb_primary_metastatic_data,
    count_tests = cbio_candidate_mutation_count_primary_metastatic_tests)
cbio_candidate_gene_count_tmb_correlation_plot <-
  plot_cbio_candidate_gene_count_tmb_correlation(
    correlation_data = cbio_candidate_gene_count_tmb_correlation_data,
    correlation_test = cbio_candidate_gene_count_tmb_correlation_test)
write_csv(
  cbio_clinical_long,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_clinical_metadata_long_data.csv"))
write_csv(
  cbio_sample_clinical,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_clinical_sample_metadata_data.csv"))
write_csv(
  CBIO_CLINICAL_ENDPOINT_DEFINITIONS,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_clinical_endpoint_definitions_data.csv"))
write_csv(
  cbio_clinical_endpoint_availability,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_clinical_endpoint_availability_data.csv"))
write_csv(
  cbio_clinical_endpoint_rows,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_clinical_endpoint_sample_assignments_data.csv"))
write_csv(
  CBIO_FOCUSED_SENSITIVITY_DEFINITIONS,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_KNDC1_FOXP4_sensitivity_definitions_data.csv"))
write_csv(
  cbio_focused_all_clinical_associations,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_KNDC1_FOXP4_all_clinical_associations_data.csv"))
write_csv(
  CBIO_TREATMENT_RESPONSE_AUDIT,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_treatment_response_availability_data.csv"))
write_csv(
  cbio_candidate_burden,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_burden_data.csv"))
write_csv(
  cbio_candidate_burden_associations,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_burden_associations_data.csv"))
write_csv(
  cbio_candidate_recurrence_unadjusted_associations,
  file.path(
    OUTPUT_DIR,
    paste0(
      "Figure4_cbioportal_candidate_recurrence_unadjusted",
      "_associations_data.csv")))
write_csv(
  cbio_candidate_recurrence_adjusted_associations,
  file.path(
    OUTPUT_DIR,
    paste0(
      "Figure4_cbioportal_candidate_recurrence_study_adjusted",
      "_associations_data.csv")))
write_csv(
  cbio_candidate_tmb_recurrence_summary,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_TMB_recurrence_summary_data.csv"))
write_csv(
  cbio_candidate_tmb_recurrence_tests,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_TMB_recurrence_tests_data.csv"))
write_csv(
  cbio_candidate_mutation_count_recurrence_summary,
  file.path(
    OUTPUT_DIR,
    paste0(
      "Figure4_cbioportal_candidate_mutation_count_recurrence",
      "_summary_data.csv")))
write_csv(
  cbio_candidate_mutation_count_recurrence_tests,
  file.path(
    OUTPUT_DIR,
    paste0(
      "Figure4_cbioportal_candidate_mutation_count_recurrence",
      "_tests_data.csv")))
write_csv(
  cbio_candidate_tmb_primary_metastatic_summary,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_TMB_primary_metastatic_summary_data.csv"))
write_csv(
  cbio_candidate_tmb_primary_metastatic_tests,
  file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_TMB_primary_metastatic_tests_data.csv"))
write_csv(
  cbio_candidate_mutation_count_primary_metastatic_summary,
  file.path(
    OUTPUT_DIR,
    paste0(
      "Figure4_cbioportal_candidate_mutation_count_primary_metastatic",
      "_summary_data.csv")))
write_csv(
  cbio_candidate_mutation_count_primary_metastatic_tests,
  file.path(
    OUTPUT_DIR,
    paste0(
      "Figure4_cbioportal_candidate_mutation_count_primary_metastatic",
      "_tests_data.csv")))
write_csv(
  cbio_candidate_gene_count_tmb_correlation_test,
  file.path(
    OUTPUT_DIR,
    paste0(
      "Figure4_cbioportal_candidate_gene_count",
      "_TMB_spearman_tests_data.csv")))

export_ggplot_data(
  plot = cbio_candidate_clinical_association_plot,
  data = cbio_candidate_clinical_associations,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_clinical_associations"),
  width = 250,
  height = max(210, 6.5 * length(clinical_gene_levels) + 40))
export_ggplot_data(
  plot = cbio_primary_metastatic_prevalence_plot,
  data = cbio_primary_metastatic_prevalence,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_primary_metastatic_prevalence"),
  width = 350,
  height = max(210, 6.5 * length(clinical_gene_levels) + 40))
export_ggplot_data(
  plot = cbio_primary_metastatic_candidate_volcano_plot,
  data = cbio_primary_metastatic_highlight_data,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_primary_metastatic_candidate_effects_volcano"),
  width = 180,
  height = 145)
export_ggplot_data(
  plot = cbio_candidate_recurrence_volcano_plot,
  data = cbio_candidate_recurrence_adjusted_associations,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_recurrence_effects_volcano"),
  width = 185,
  height = 150)
export_ggplot_data(
  plot = cbio_candidate_tmb_recurrence_plot,
  data = cbio_candidate_recurrence_boxplot_data,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_TMB_recurrence"),
  width = 285,
  height = 155)
export_ggplot_data(
  plot = cbio_candidate_mutation_count_recurrence_plot,
  data = cbio_candidate_recurrence_boxplot_data,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_mutation_count_recurrence"),
  width = 285,
  height = 155)
export_ggplot_data(
  plot = cbio_primary_metastatic_all_gene_volcano_plot,
  data = cbio_all_gene_primary_metastatic_associations,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_primary_metastatic_all_mutated_genes_volcano"),
  width = 190,
  height = 150)
export_ggplot_data(
  plot = cbio_focused_primary_metastatic_prevalence_plot,
  data = cbio_focused_primary_metastatic_prevalence,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_KNDC1_FOXP4_primary_metastatic_prevalence"),
  width = 210,
  height = 115)
export_ggplot_data(
  plot = cbio_focused_primary_metastatic_sensitivity_plot,
  data = cbio_focused_primary_metastatic_associations,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_KNDC1_FOXP4_primary_metastatic_sensitivity"),
  width = 220,
  height = 135)
export_ggplot_data(
  plot = cbio_candidate_burden_association_plot,
  data = cbio_candidate_burden_endpoint_data,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_burden_distributions"),
  width = 240,
  height = 235)
export_ggplot_data(
  plot = cbio_candidate_tmb_primary_metastatic_plot,
  data = cbio_candidate_tmb_primary_metastatic_data,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_TMB_primary_metastatic"),
  width = 225,
  height = 155)
export_ggplot_data(
  plot = cbio_candidate_mutation_count_primary_metastatic_plot,
  data = cbio_candidate_tmb_primary_metastatic_data,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_mutation_count_primary_metastatic"),
  width = 225,
  height = 155)
ggplot2::set_last_plot(cbio_candidate_gene_count_tmb_correlation_plot)
export_plot_data(
  data = cbio_candidate_gene_count_tmb_correlation_data,
  file_name = file.path(
    OUTPUT_DIR,
    "Figure4_cbioportal_candidate_gene_count_vs_TMB_spearman"))
}

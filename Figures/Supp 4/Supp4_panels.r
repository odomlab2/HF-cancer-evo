# 02/06/2026
# Author: Yoav Avi-Guy
# Purpose: Supplementary Figure 4 panels - Landscape of initiation and 
# progression mutated genes

#-------------------------------------------------------------------------------
# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(UpSetR)
library(gghalves)
library(ggpattern)
repo_root <- normalizePath(Sys.getenv("HF_SCC_ROOT", unset = "."), mustWork = TRUE)
output_root <- Sys.getenv("HF_SCC_OUTPUT_ROOT")
if (!nzchar(output_root)) stop("Set HF_SCC_OUTPUT_ROOT to an external output directory.")

# Colours
source(file.path(repo_root, "TES", "01_source", "01_plotting.r"))

# Loading files
tes_metadata <- readRDS(file.path(repo_root, "TES/02_data/01_metadata/sample_metadata.rds"))
wes_metadata <- readRDS(file.path(repo_root, "WES/02_data/01_metadata/sample_metadata.rds"))
tes_unique <- readRDS(file.path(repo_root, "TES/02_data/02_processed/mutations_unique_dp.rds"))
wes_unique <- readRDS(file.path(repo_root, "WES/02_data/02_processed/mutations_unique_dp.rds"))
wes_dnds <- readRDS(file.path(repo_root, "WES/02_data/02_processed/dnds_category_results.rds"))
figure4_candidate_genes <- readr::read_csv(
  file.path(repo_root, "Figures/Figure 4/Figure4_candidate_genes_data.csv"),
  show_col_types = FALSE)
tes_outliers <- readRDS(file.path(repo_root, "TES/02_data/02_processed/technical_outliers.rds"))
wes_outliers <- readRDS(file.path(repo_root, "WES/02_data/02_processed/technical_outliers.rds"))

#-------------------------------------------------------------------------------
### Parameters ###
#-------------------------------------------------------------------------------
supp4_dir <- file.path(output_root, "Figures", "Supp 4")
dir.create(supp4_dir, recursive = TRUE, showWarnings = FALSE)
refcds_path <- Sys.getenv("HF_SCC_REFCDS")
if (!nzchar(refcds_path)) stop("Set HF_SCC_REFCDS to the reference CDS RDA file.")
cbio_mutated_genes_file <- file.path(repo_root, "Figures/Figure 4/cBioPortal_Mutated_Genes.txt")
min_gt_af <- 0.01

manual_human_orthologues <- tibble::tibble(
  gene_symbol = c(
    "C130026I21Rik",
    "Fntb",
    "Hras",
    "Notch1",
    "Tp53",
    "Ttll3",
    "Zrsr2"),
  human_symbols = c(
    "SP140",
    "FNTB",
    "HRAS",
    "NOTCH1",
    "TP53",
    "TTLL3",
    "ZRSR2"))

shared_gene_group_levels <- c(
  "Morphologically normal",
  "HF papilloma",
  "HF SCC",
  "Skin papilloma",
  "Skin SCC")

ras_gene_levels <- c("Hras", "Kras", "Nras")
ras_category_levels <- c(
  "Morphologically normal",
  "Papilloma",
  "SCC")
ras_category_colours <- c(
  SCC = unname(col_palette$category["SCC"]),
  "Morphologically normal" = unname(col_palette$category["Visually normal"]),
  Papilloma = unname(col_palette$category["Papilloma"]),
  Acetone = unname(col_palette$category["Acetone"]))
ras_tissue_levels <- c("Hair follicle", "Skin")

dnds_category_levels <- c("Visually normal", "Papilloma")
dnds_class_levels <- c("Missense", "Nonsense", "Indel/DBS")
dnds_positive_omega_cutoff <- 1.6
dnds_q_cutoff <- 0.05

candidate_group_levels <- c("Initiation (I)", "I+P", "Predisposition")
candidate_mutation_type_levels <- c(
  "Missense",
  "Insertion",
  "Deletion",
  "Stop gain/loss",
  "DBS",
  "Splice")
candidate_mutation_type_colours <- c(
  Missense = "#3D1B1B",
  Insertion = "#7EA3CF",
  Deletion = "#436503",
  "Stop gain/loss" = "#9E2F2F",
  DBS = "#2F9A8B",
  "Splice" = "#6F5A8C")
candidate_mutation_bar_height_mm <- 4
candidate_mutation_random_iterations <- 10000L
candidate_mutation_random_n_genes <- 30L
candidate_mutation_random_seed <- 20260602L

#-------------------------------------------------------------------------------
### Small helpers ###
#-------------------------------------------------------------------------------
export_pdf_plot_data <- function(
  data,
  file_name,
  draw,
  width,
  height,
  cols = NULL) {

  out <- if (is.null(cols)) data else dplyr::select(data, any_of(cols))
  readr::write_csv(out, paste0(file_name, "_data.csv"))

  pdf(paste0(file_name, ".pdf"), width = width, height = height)
  draw()
  dev.off()

  invisible(out)
}

filter_last_cycle <- function(df) {
  if ("last_cycle" %in% names(df)) {
    df <- df %>% filter(last_cycle == TRUE)
  }

  df
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

is_priority_nonsynonymous <- function(consequence, impact) {
  consequence <- as.character(consequence)
  impact <- as.character(impact)

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

prepare_mutations <- function(
  raw_tbl,
  outliers,
  dataset_label,
  use_last_cycle = FALSE,
  min_gt_af_filter = min_gt_af) {

  tbl <- raw_tbl %>%
    mutate(gt_AF = suppressWarnings(as.numeric(gt_AF)))

  if (!"time" %in% names(tbl)) {
    tbl$time <- NA_character_
  }

  if (use_last_cycle) {
    tbl <- filter_last_cycle(tbl)
  }

  tbl <- tbl %>%
    filter(
      !is.na(sample_name),
      !is.na(gt_AF),
      !sample_name %in% outliers$sample_name)

  if (!is.null(min_gt_af_filter)) {
    tbl <- tbl %>% filter(gt_AF >= min_gt_af_filter)
  }

  tbl %>%
    mutate(
      dataset = dataset_label,
      sample_uid = paste(dataset_label, sample_name, sep = "::"),
      mutation_id = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
    filter(
      !(
        (
          treatment == "Acetone" |
            (treatment == "DT" & category == "Visually normal")
        ) &
          IMPACT %in% c("LOW", "MODIFIER")
      ))
}

build_nonsynonymous_base <- function(presence_tbl) {
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

assign_shared_gene_group <- function(tissue, treatment, category) {
  case_when(
    treatment == "DT" & category == "Visually normal" ~
      "Morphologically normal",
    treatment == "DT" & tissue == "Hair follicle" & category == "Papilloma" ~
      "HF papilloma",
    treatment == "DT" & tissue == "Hair follicle" & category == "SCC" ~
      "HF SCC",
    treatment == "DT" & tissue == "Skin" & category == "Papilloma" ~
      "Skin papilloma",
    treatment == "DT" & tissue == "Skin" & category == "SCC" ~
      "Skin SCC",
    TRUE ~ NA_character_)
}

make_shared_gene_groups <- function(mutations, group_levels) {
  shared_gene_long <- mutations %>%
    mutate(
      gene_symbol = as.character(SYMBOL),
      group = assign_shared_gene_group(
        as.character(tissue),
        as.character(treatment),
        as.character(category))) %>%
    filter(!is.na(group), !is.na(gene_symbol), gene_symbol != "") %>%
    distinct(
      gene_symbol,
      group,
      dataset,
      sample_uid,
      sample_name,
      tissue,
      treatment,
      category,
      mutation_id,
      CHROM,
      POS,
      REF,
      ALT,
      Consequence,
      IMPACT,
      gt_AF)

  shared_gene_matrix <- shared_gene_long %>%
    distinct(gene_symbol, group) %>%
    mutate(present = 1L) %>%
    pivot_wider(
      names_from = group,
      values_from = present,
      values_fill = 0L)

  shared_gene_matrix[setdiff(group_levels, names(shared_gene_matrix))] <- 0L

  shared_gene_export <- shared_gene_matrix %>%
    dplyr::select(gene_symbol, all_of(group_levels)) %>%
    as.data.frame()
  rownames(shared_gene_export) <- shared_gene_export$gene_symbol

  shared_gene_summary <- shared_gene_long %>%
    group_by(gene_symbol, group) %>%
    summarise(
      n_samples = n_distinct(sample_uid),
      n_mutations = n_distinct(mutation_id),
      datasets = paste(sort(unique(dataset)), collapse = "; "),
      samples = paste(sort(unique(sample_name)), collapse = "; "),
      max_gt_AF = max(gt_AF, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    arrange(gene_symbol, group)

  list(
    matrix = shared_gene_export %>%
      dplyr::select(all_of(group_levels)) %>%
      as.data.frame(),
    export = shared_gene_export,
    long = shared_gene_long,
    summary = shared_gene_summary)
}

save_shared_gene_upset <- function(shared_groups, group_levels, file_stub) {
  readr::write_csv(
    shared_groups$summary,
    paste0(file_stub, "_long_data.csv"))

  upset_plot <- UpSetR::upset(
    shared_groups$matrix,
    sets = group_levels,
    nsets = length(group_levels),
    nintersects = 2 ^ length(group_levels),
    keep.order = TRUE,
    order.by = c("degree", "freq"),
    mb.ratio = c(0.65, 0.35),
    mainbar.y.label = "# mutated genes across group intersections",
    sets.x.label = "# mutated genes in group",
    main.bar.color = "black",
    sets.bar.color = "black",
    matrix.color = "black",
    shade.color = "white",
    shade.alpha = 0,
    text.scale = c(1.5, 1.4, 1.2, 1.2, 1.4, 1.2),
    point.size = 3,
    line.size = 1)

  export_pdf_plot_data(
    data = shared_groups$export,
    file_name = file_stub,
    draw = function() print(upset_plot),
    width = 6.5,
    height = 5)

  if (interactive()) print(upset_plot)

  invisible(upset_plot)
}

sel_long_df <- function(dndscv_res) {
  sel <- dndscv_res$sel_cv

  tibble::tibble(
    gene = rep(sel$gene_name, 3),
    dnds_class = rep(dnds_class_levels, each = nrow(sel)),
    omega = c(sel$wmis_cv, sel$wnon_cv, sel$wind_cv),
    q = c(sel$qmis_cv, sel$qtrunc_cv, sel$qind_cv))
}

prep_positive_dnds <- function(
  dnds_list,
  category_levels = dnds_category_levels,
  w_pos = dnds_positive_omega_cutoff,
  q_cut = dnds_q_cutoff) {

  dnds_long <- lapply(category_levels, function(category_name) {
    if (!category_name %in% names(dnds_list)) {
      return(tibble::tibble())
    }

    sel_long_df(dnds_list[[category_name]]) %>%
      mutate(
        category = category_name,
        sig = q < q_cut,
        positive_selection = omega > w_pos & sig)
  }) %>%
    bind_rows()

  selected_genes <- dnds_long %>%
    filter(positive_selection) %>%
    distinct(category, gene)

  dnds_long %>%
    inner_join(selected_genes, by = c("category", "gene")) %>%
    group_by(category, gene) %>%
    mutate(gene_score = max(omega, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      category = factor(category, levels = category_levels),
      dnds_class = factor(dnds_class, levels = dnds_class_levels)) %>%
    arrange(category, desc(gene_score), gene) %>%
    mutate(
      gene_label = factor(
        paste(category, gene, sep = "::"),
        levels = unique(paste(category, gene, sep = "::"))))
}

split_amino_acids <- function(x) {
  if (is.na(x) || x == "") {
    return(character(0))
  }

  strsplit(x, "", fixed = FALSE)[[1]]
}

ras_amino_acid_label <- function(protein_position, amino_acids) {
  protein_position <- coalesce(as.character(protein_position), "")
  amino_acids <- coalesce(as.character(amino_acids), "")

  start_position <- suppressWarnings(as.numeric(
    str_match(protein_position, "^(\\d+)")[, 2]))
  if (is.na(start_position) || amino_acids == "") {
    return(tibble::tibble(
      codon_position = NA_real_,
      amino_acid_label = NA_character_))
  }

  amino_acids <- str_replace_all(amino_acids, "\\*", "X")
  aa_parts <- str_split_fixed(amino_acids, "/", 2)
  ref_chars <- split_amino_acids(aa_parts[, 1])
  alt_chars <- split_amino_acids(aa_parts[, 2])
  max_len <- max(length(ref_chars), length(alt_chars))

  if (max_len == 0) {
    return(tibble::tibble(
      codon_position = NA_real_,
      amino_acid_label = NA_character_))
  }

  length(ref_chars) <- max_len
  length(alt_chars) <- max_len
  changed_index <- which(
    (is.na(ref_chars) != is.na(alt_chars)) |
      (!is.na(ref_chars) & !is.na(alt_chars) & ref_chars != alt_chars))[1]

  if (is.na(changed_index)) {
    changed_index <- 1L
  }

  codon_position <- start_position + changed_index - 1
  ref_aa <- ref_chars[changed_index]
  alt_aa <- alt_chars[changed_index]
  ref_label <- ifelse(is.na(ref_aa) || ref_aa == "-", "", ref_aa)
  alt_label <- ifelse(is.na(alt_aa) || alt_aa %in% c("", "-"), "X", alt_aa)

  tibble::tibble(
    codon_position = codon_position,
    amino_acid_label = paste0(ref_label, codon_position, alt_label))
}

build_ras_mutation_map_data <- function(mutation_tbl) {
  if (nrow(mutation_tbl) == 0) {
    return(tibble::tibble(
      gene_symbol = factor(levels = ras_gene_levels),
      codon_position = numeric(),
      amino_acid_label = character(),
      mutation = factor(),
      morphology = factor(levels = ras_category_levels),
      n_samples = integer(),
      n_mutations = integer(),
      tissues = character(),
      datasets = character(),
      sample_names = character(),
      mutation_ids = character(),
      protein_positions = character(),
      amino_acids = character(),
      consequences = character(),
      mean_gt_AF = numeric(),
      max_gt_AF = numeric()))
  }

  label_df <- dplyr::bind_rows(lapply(seq_len(nrow(mutation_tbl)), function(i) {
    ras_amino_acid_label(
      mutation_tbl$Protein_position[[i]],
      mutation_tbl$Amino_acids[[i]])
  }))

  mutation_tbl %>%
    dplyr::bind_cols(label_df) %>%
    mutate(
      gene_symbol = as.character(SYMBOL),
      morphology = recode(
        as.character(category),
        "Visually normal" = "Morphologically normal",
        .default = as.character(category))) %>%
    filter(
      gene_symbol %in% ras_gene_levels,
      morphology %in% ras_category_levels,
      !is.na(codon_position),
      !is.na(amino_acid_label)) %>%
    mutate(
      gene_symbol = factor(gene_symbol, levels = ras_gene_levels),
      morphology = factor(morphology, levels = ras_category_levels)) %>%
    group_by(
      gene_symbol,
      codon_position,
      amino_acid_label,
      morphology) %>%
    summarise(
      n_samples = n_distinct(sample_uid),
      n_mutations = n_distinct(mutation_id),
      tissues = paste(sort(unique(tissue)), collapse = "; "),
      datasets = paste(sort(unique(dataset)), collapse = "; "),
      sample_names = paste(sort(unique(sample_name)), collapse = "; "),
      mutation_ids = paste(sort(unique(mutation_id)), collapse = "; "),
      protein_positions = paste(sort(unique(Protein_position)), collapse = "; "),
      amino_acids = paste(sort(unique(Amino_acids)), collapse = "; "),
      consequences = paste(sort(unique(Consequence)), collapse = "; "),
      mean_gt_AF = mean(gt_AF, na.rm = TRUE),
      max_gt_AF = max(gt_AF, na.rm = TRUE),
      .groups = "drop") %>%
    arrange(gene_symbol, codon_position, amino_acid_label, morphology) %>%
    mutate(
      mutation = paste(gene_symbol, amino_acid_label, sep = "::"),
      mutation = factor(mutation, levels = unique(mutation))) %>%
    dplyr::select(
      gene_symbol,
      codon_position,
      amino_acid_label,
      mutation,
      morphology,
      everything())
}

plot_ras_mutation_map <- function(plot_df) {
  observed_morphologies <- ras_category_levels[
    ras_category_levels %in% as.character(plot_df$morphology)]
  max_samples <- if (nrow(plot_df) == 0) {
    0
  } else {
    plot_df %>%
      group_by(gene_symbol, mutation) %>%
      summarise(n_samples = sum(n_samples), .groups = "drop") %>%
      summarise(max_samples = max(n_samples)) %>%
      pull(max_samples)
  }

  ggplot(
    plot_df,
    aes(x = mutation, y = n_samples, fill = morphology)) +
    geom_col(
      width = 0.75,
      colour = "black",
      linewidth = 0.2,
      position = position_stack(reverse = TRUE)) +
    facet_grid(
      . ~ gene_symbol,
      scales = "free_x",
      space = "free_x",
      drop = FALSE) +
    scale_fill_manual(
      values = ras_category_colours,
      breaks = observed_morphologies,
      drop = TRUE) +
    scale_x_discrete(
      labels = function(x) sub("^.*::", "", x),
      expand = expansion(add = 0.45)) +
    scale_y_continuous(
      breaks = unique(c(0, max_samples)),
      labels = function(x) format(x, scientific = FALSE, trim = TRUE),
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "Mutation",
      y = "# samples mutated",
      fill = "Morphology") +
    my_theme +
    theme(
      text = element_text(size = 8),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      strip.text = element_text(face = "italic", size = 8))
}

plot_positive_dnds <- function(dnds_df) {
  ggplot(
    dnds_df,
    aes(x = gene_label, y = omega, fill = dnds_class, alpha = sig)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    # geom_hline(
    #   yintercept = 1,
    #   linetype = "dashed",
    #   linewidth = 0.4,
    #   alpha = 0.6) +
    facet_grid(. ~ category, scales = "free_x", space = "free_x") +
    scale_x_discrete(
      labels = function(x) str_replace(x, ".*::", "")) +
    scale_y_continuous(
      labels = my_label,
      expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = col_palette$dnds_class) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25)) +
    labs(
      x = "",
      y = "dN/dS ratio",
      fill = "") +
    guides(alpha = "none") +
    my_theme +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(hjust = 1))
}

candidate_mutation_type <- function(consequence, variant_class, ref, alt) {
  consequence <- coalesce(as.character(consequence), "")
  variant_class <- coalesce(as.character(variant_class), "")
  ref <- coalesce(as.character(ref), "")
  alt <- coalesce(as.character(alt), "")

  case_when(
    str_detect(consequence, "stop_gained|stop_lost") ~ "Stop gain/loss",
    variant_class == "insertion" | nchar(ref) < nchar(alt) ~ "Insertion",
    variant_class == "deletion" | nchar(ref) > nchar(alt) ~ "Deletion",
    variant_class %in% c("substitution", "MNV", "DNV") |
      (nchar(ref) == 2 & nchar(alt) == 2) ~ "DBS",
    str_detect(consequence, "missense_variant") ~ "Missense",
    TRUE ~ "Splice")
}

candidate_gene_order <- function(candidate_matrix_tbl) {
  candidate_matrix_tbl %>%
    distinct(
      gene_symbol,
      row_label,
      gene_label,
      candidate_group,
      lollipop_order) %>%
    arrange(lollipop_order) %>%
    mutate(
      row_label = factor(row_label, levels = rev(unique(row_label))),
      candidate_group = factor(candidate_group, levels = candidate_group_levels))
}

build_candidate_mutation_type_counts <- function(
  mutation_tbl,
  candidate_matrix_tbl) {

  gene_order <- candidate_gene_order(candidate_matrix_tbl)
  mutation_type_rank <- setNames(
    seq_along(candidate_mutation_type_levels),
    candidate_mutation_type_levels)

  candidate_mutations <- mutation_tbl %>%
    filter(
      treatment == "DT",
      SYMBOL %in% gene_order$gene_symbol) %>%
    mutate(
      gene_symbol = as.character(SYMBOL),
      mutation_type = candidate_mutation_type(
        Consequence,
        VARIANT_CLASS,
        REF,
        ALT),
      mutation_type_rank = mutation_type_rank[mutation_type]) %>%
    group_by(sample_name, gene_symbol, mutation_id) %>%
    arrange(mutation_type_rank, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()

  observed_mutation_type_levels <- candidate_mutation_type_levels[
    candidate_mutation_type_levels %in% candidate_mutations$mutation_type]

  count_df <- candidate_mutations %>%
    count(gene_symbol, mutation_type, name = "n_mutations") %>%
    right_join(gene_order, by = "gene_symbol") %>%
    mutate(mutation_type = factor(
      mutation_type,
      levels = observed_mutation_type_levels)) %>%
    tidyr::complete(
      nesting(
        gene_symbol,
        row_label,
        gene_label,
        candidate_group,
        lollipop_order),
      mutation_type = factor(
        observed_mutation_type_levels,
        levels = observed_mutation_type_levels),
      fill = list(n_mutations = 0L)) %>%
    group_by(gene_symbol) %>%
    mutate(total_mutations = sum(n_mutations, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      row_label = factor(row_label, levels = levels(gene_order$row_label)),
      candidate_group = factor(candidate_group, levels = candidate_group_levels),
      mutation_type = factor(
        mutation_type,
        levels = observed_mutation_type_levels)) %>%
    arrange(lollipop_order, mutation_type)

  list(
    counts = count_df,
    mutations = candidate_mutations,
    mutation_type_levels = observed_mutation_type_levels)
}

plot_candidate_mutation_type_counts <- function(count_obj) {
  plot_df <- count_obj$counts
  mutation_type_levels <- count_obj$mutation_type_levels

  ggplot(
    plot_df,
    aes(x = n_mutations, y = row_label, fill = mutation_type)) +
    geom_col(width = 0.75) +
    facet_grid(
      candidate_group ~ .,
      scales = "free_y",
      space = "free_y") +
    scale_x_continuous(
      breaks = scales::breaks_width(2),
      expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(
      values = candidate_mutation_type_colours[mutation_type_levels],
      limits = mutation_type_levels,
      breaks = mutation_type_levels,
      drop = TRUE) +
    labs(
      x = "# mutations",
      y = "",
      fill = "") +
    my_theme +
    theme(
      legend.position = "bottom",
      strip.text.y = element_text(angle = 0),
      axis.text.y = element_text(size = 8))
}

export_candidate_mutation_type_plot <- function(
  plot,
  count_obj,
  file_name,
  bar_height_mm = candidate_mutation_bar_height_mm,
  width = 95,
  extra_height = 50,
  cols = NULL) {

  out <- if (is.null(cols)) {
    count_obj$counts
  } else {
    dplyr::select(count_obj$counts, any_of(cols))
  }

  facet_heights <- count_obj$counts %>%
    distinct(candidate_group, row_label) %>%
    count(candidate_group, name = "n_genes") %>%
    mutate(candidate_group = factor(
      candidate_group,
      levels = candidate_group_levels)) %>%
    arrange(candidate_group) %>%
    pull(n_genes) * bar_height_mm

  plot <- plot +
    ggh4x::force_panelsizes(rows = grid::unit(facet_heights, "mm"))

  readr::write_csv(out, paste0(file_name, "_data.csv"))
  ggsave(
    filename = paste0(file_name, ".pdf"),
    plot = plot,
    width = width,
    height = sum(facet_heights) + extra_height,
    units = "mm")

  invisible(out)
}

build_mutation_type_gene_matrix <- function(
  mutation_tbl,
  mutation_type_levels = candidate_mutation_type_levels) {

  mutation_type_rank <- setNames(
    seq_along(mutation_type_levels),
    mutation_type_levels)

  mutation_tbl %>%
    filter(
      treatment == "DT",
      !is.na(SYMBOL),
      SYMBOL != "") %>%
    mutate(
      gene_symbol = as.character(SYMBOL),
      mutation_type = candidate_mutation_type(
        Consequence,
        VARIANT_CLASS,
        REF,
        ALT),
      mutation_type_rank = mutation_type_rank[mutation_type]) %>%
    group_by(sample_name, gene_symbol, mutation_id) %>%
    arrange(mutation_type_rank, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    count(gene_symbol, mutation_type, name = "n_mutations") %>%
    complete(
      gene_symbol,
      mutation_type = mutation_type_levels,
      fill = list(n_mutations = 0L)) %>%
    pivot_wider(
      names_from = mutation_type,
      values_from = n_mutations,
      values_fill = 0L) %>%
    arrange(gene_symbol)
}

run_candidate_mutation_type_null <- function(
  mutation_tbl,
  observed_counts,
  n_iterations = candidate_mutation_random_iterations,
  n_genes = candidate_mutation_random_n_genes,
  seed = candidate_mutation_random_seed,
  mutation_type_levels = candidate_mutation_type_levels) {

  gene_matrix <- build_mutation_type_gene_matrix(
    mutation_tbl = mutation_tbl,
    mutation_type_levels = mutation_type_levels)
  gene_symbols <- gene_matrix$gene_symbol

  if (length(gene_symbols) < n_genes) {
    stop("Random gene universe contains fewer genes than requested.")
  }

  background_present <- mutation_type_levels[
    colSums(as.matrix(gene_matrix[, mutation_type_levels, drop = FALSE])) > 0]
  observed_present <- observed_counts %>%
    group_by(mutation_type) %>%
    summarise(
      observed_count = sum(n_mutations, na.rm = TRUE),
      .groups = "drop") %>%
    filter(observed_count > 0) %>%
    pull(mutation_type) %>%
    as.character()
  mutation_type_levels <- mutation_type_levels[
    mutation_type_levels %in% union(background_present, observed_present)]

  if (length(mutation_type_levels) == 0) {
    stop("No mutation types were present in the random gene universe.")
  }

  count_matrix <- as.matrix(gene_matrix[, mutation_type_levels, drop = FALSE])
  rownames(count_matrix) <- gene_symbols

  set.seed(seed)

  sampled_indices <- replicate(
    n_iterations,
    sample.int(length(gene_symbols), n_genes, replace = FALSE))
  n_sampled_genes <- apply(sampled_indices, 2, function(idx) length(unique(idx)))

  if (any(n_sampled_genes != n_genes)) {
    stop("At least one random iteration did not sample distinct genes.")
  }

  null_matrix <- t(vapply(seq_len(n_iterations), function(iteration) {
    sampled_genes <- sampled_indices[, iteration]
    colSums(count_matrix[sampled_genes, , drop = FALSE], na.rm = TRUE)
  }, numeric(length(mutation_type_levels))))
  colnames(null_matrix) <- mutation_type_levels

  null_long <- tibble::as_tibble(null_matrix) %>%
    mutate(
      iteration = row_number(),
      n_sampled_genes = n_sampled_genes) %>%
    pivot_longer(
      cols = all_of(mutation_type_levels),
      names_to = "mutation_type",
      values_to = "null_count") %>%
    mutate(mutation_type = factor(mutation_type, levels = mutation_type_levels))

  observed_long <- observed_counts %>%
    group_by(mutation_type) %>%
    summarise(
      observed_count = sum(n_mutations, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(mutation_type = as.character(mutation_type)) %>%
    right_join(
      tibble::tibble(mutation_type = mutation_type_levels),
      by = "mutation_type") %>%
    mutate(
      observed_count = replace_na(observed_count, 0L),
      mutation_type = factor(mutation_type, levels = mutation_type_levels))

  summary_tbl <- null_long %>%
    left_join(observed_long, by = "mutation_type") %>%
    group_by(mutation_type) %>%
    summarise(
      n_iterations = n_distinct(iteration),
      n_sampled_genes = dplyr::first(n_sampled_genes),
      n_random_genes = n_genes,
      n_universe_genes = length(gene_symbols),
      observed_count = dplyr::first(observed_count),
      null_mean = mean(null_count, na.rm = TRUE),
      null_median = median(null_count, na.rm = TRUE),
      `null_q2.5` = unname(quantile(null_count, 0.025, na.rm = TRUE)),
      `null_q97.5` = unname(quantile(null_count, 0.975, na.rm = TRUE)),
      empirical_p_upper = (sum(null_count >= observed_count, na.rm = TRUE) + 1) /
        (n_iterations + 1),
      observed_percentile = 100 * mean(null_count <= observed_count, na.rm = TRUE),
      .groups = "drop")

  list(
    long = null_long,
    observed = observed_long,
    summary = summary_tbl,
    gene_matrix = gene_matrix,
    mutation_type_levels = mutation_type_levels)
}

plot_candidate_mutation_type_null <- function(null_obj) {
  mutation_type_levels <- null_obj$mutation_type_levels

  x_lookup <- tibble::tibble(
    mutation_type = factor(mutation_type_levels, levels = mutation_type_levels),
    mutation_type_index = seq_along(mutation_type_levels))

  null_df <- null_obj$long %>%
    mutate(mutation_type = factor(mutation_type, levels = mutation_type_levels)) %>%
    left_join(x_lookup, by = "mutation_type") %>%
    filter(null_count <= 100)

  summary_df <- null_obj$summary %>%
    mutate(
      mutation_type = factor(mutation_type, levels = mutation_type_levels)) %>%
    left_join(x_lookup, by = "mutation_type")

  significant_df <- summary_df %>%
    filter(!is.na(empirical_p_upper), empirical_p_upper < 0.05) %>%
    mutate(
      p_label = if_else(
        empirical_p_upper < 0.001,
        "italic(P) < 0.001",
        paste0(
          "italic(P) == ",
          formatC(empirical_p_upper, format = "f", digits = 3))),
      p_y = pmax(`null_q97.5`, observed_count, na.rm = TRUE) + 5,
      p_y = if_else(p_y > 10 & p_y < 20, 21, p_y),
      p_hjust = case_when(
        mutation_type_index == min(x_lookup$mutation_type_index) ~ 0,
        mutation_type_index == max(x_lookup$mutation_type_index) ~ 1,
        .default = 0.5))

  base_plot <- ggplot(
    null_df,
    aes(
      x = mutation_type_index,
      y = null_count,
      group = mutation_type)) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      fill = "grey70",
      colour = "grey40",
      linewidth = 0.25) +
    geom_segment(
      data = summary_df,
      aes(
        x = mutation_type_index - 0.3,
        xend = mutation_type_index + 0.3,
        y = observed_count,
        yend = observed_count),
      inherit.aes = FALSE,
      colour = "red",
      linewidth = 0.5) +
    scale_x_continuous(
      breaks = x_lookup$mutation_type_index,
      labels = mutation_type_levels,
      limits = c(0.2, max(x_lookup$mutation_type_index) + 0.8),
      expand = expansion(mult = c(0, 0))) +
    labs(
      x = "",
      y = "# mutations per 30 genes") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5))

  upper_plot <- base_plot +
    geom_text(
      data = significant_df,
      aes(
        x = mutation_type_index,
        y = p_y,
        label = p_label,
        hjust = p_hjust),
      inherit.aes = FALSE,
      size = 8 / .pt,
      vjust = 0,
      parse = TRUE) +
    annotate(
      "segment",
      x = 0.2,
      xend = 0.45,
      y = 20,
      yend = 30,
      linewidth = 0.4) +
    scale_y_continuous(
      breaks = c(20, 50, 100),
      expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(20, 100), clip = "on") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(15, 5.5, 1.5, 5.5, unit = "pt"))

  lower_plot <- base_plot +
    annotate(
      "segment",
      x = 0.2,
      xend = 0.45,
      y = 8.5,
      yend = 10,
      linewidth = 0.4) +
    scale_y_continuous(
      breaks = c(0, 5, 10),
      expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, 10), clip = "on") +
    theme(plot.margin = margin(1.5, 5.5, 15, 5.5, unit = "pt"))

  patchwork::wrap_plots(
    upper_plot,
    lower_plot,
    ncol = 1,
    heights = c(1, 1),
    axis_titles = "collect_y")
}

build_candidate_gene_sample_counts <- function(
  mutation_tbl,
  candidate_matrix_tbl) {

  candidate_symbols <- candidate_matrix_tbl %>%
    distinct(gene_symbol) %>%
    pull(gene_symbol)

  sample_universe <- mutation_tbl %>%
    filter(
      treatment == "DT",
      tissue %in% c("Hair follicle", "Skin"),
      !is.na(sample_uid),
      !is.na(sample_name)) %>%
    distinct(sample_uid, sample_name, dataset, tissue, category) %>%
    mutate(
      tissue = factor(tissue, levels = c("Hair follicle", "Skin")),
      category = factor(category, levels = names(col_palette$category)))

  mutated_candidate_genes <- mutation_tbl %>%
    filter(
      treatment == "DT",
      tissue %in% c("Hair follicle", "Skin"),
      SYMBOL %in% candidate_symbols) %>%
    mutate(gene_symbol = as.character(SYMBOL)) %>%
    distinct(sample_uid, gene_symbol)

  sample_universe %>%
    left_join(
      mutated_candidate_genes %>%
        count(sample_uid, name = "n_candidate_genes_mutated"),
      by = "sample_uid") %>%
    mutate(
      n_candidate_genes_mutated = replace_na(n_candidate_genes_mutated, 0L)) %>%
    arrange(tissue, n_candidate_genes_mutated, sample_name) %>%
    mutate(
      sample_order = row_number(),
      sample_uid_ordered = factor(sample_uid, levels = sample_uid))
}

plot_candidate_gene_sample_counts <- function(
  sample_counts,
  x_title = "Samples",
  y_title = "# candidate genes mutated",
  y_breaks = scales::breaks_width(2),
  legend_key_size = grid::unit(1, "lines")) {

  ggplot(
    sample_counts,
    aes(
      x = sample_uid_ordered,
      y = n_candidate_genes_mutated,
      fill = category)) +
    geom_col(width = 0.8) +
    scale_fill_manual(
      values = c(
        col_palette$category,
        "Macr. normal" = unname(col_palette$category["Visually normal"]))) +
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = x_title,
      y = y_title,
      fill = NULL) +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = legend_key_size)
}

plot_candidate_gene_sample_counts_violin <- function(sample_counts) {
  category_colours <- c(
    "Macr. normal" = unname(col_palette$category["Visually normal"]),
    "SCC" = unname(col_palette$category["SCC"]),
    "Papilloma" = unname(col_palette$category["Papilloma"]))

  ggplot(
    sample_counts,
    aes(
      x = category,
      y = n_candidate_genes_mutated,
      fill = category)) +
    geom_violin(
      trim = TRUE,
      scale = "width",
      colour = "black",
      linewidth = 0.25) +
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.2,
      size = 0.8,
      alpha = 0.8,
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 20260717)) +
    scale_fill_manual(values = category_colours, guide = "none") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 20),
      breaks = c(0, 10, 20),
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "",
      y = "# candidate genes mutated") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5))
}

load_refcds_gene_lengths <- function(refcds_file = refcds_path) {
  ref_env <- new.env(parent = emptyenv())
  load(refcds_file, envir = ref_env)

  if (!"RefCDS" %in% ls(ref_env)) {
    stop("RefCDS was not found in ", refcds_file)
  }

  tibble::tibble(
    gene_symbol = vapply(
      ref_env$RefCDS,
      function(gene_ref) as.character(gene_ref$gene_name),
      character(1)),
    gene_id = vapply(
      ref_env$RefCDS,
      function(gene_ref) as.character(gene_ref$gene_id),
      character(1)),
    protein_id = vapply(
      ref_env$RefCDS,
      function(gene_ref) as.character(gene_ref$protein_id),
      character(1)),
    gene_length_bp = vapply(
      ref_env$RefCDS,
      function(gene_ref) as.numeric(gene_ref$CDS_length),
      numeric(1))) %>%
    filter(
      !is.na(gene_symbol),
      gene_symbol != "",
      !is.na(gene_length_bp),
      gene_length_bp > 0) %>%
    distinct(gene_symbol, .keep_all = TRUE)
}

build_gene_length_mutation_counts <- function(
  mutation_tbl,
  candidate_matrix_tbl,
  gene_lengths) {

  candidate_groups <- candidate_matrix_tbl %>%
    distinct(gene_symbol, candidate_group)

  mutation_tbl %>%
    filter(
      treatment == "DT",
      !is.na(SYMBOL),
      SYMBOL != "") %>%
    mutate(gene_symbol = as.character(SYMBOL)) %>%
    distinct(sample_uid, gene_symbol, mutation_id) %>%
    count(gene_symbol, name = "n_mutations") %>%
    left_join(gene_lengths, by = "gene_symbol") %>%
    left_join(candidate_groups, by = "gene_symbol") %>%
    filter(!is.na(gene_length_bp)) %>%
    mutate(
      is_candidate = !is.na(candidate_group),
      gene_class = factor(
        if_else(is_candidate, as.character(candidate_group), "Other genes"),
        levels = c("Other genes", names(col_palette$candidate)))) %>%
    arrange(is_candidate)
}

plot_gene_length_mutation_counts <- function(gene_counts) {
  top_candidate_genes <- gene_counts %>%
    filter(is_candidate) %>%
    arrange(desc(n_mutations), gene_symbol) %>%
    slice_head(n = 5)

  label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = top_candidate_genes,
      aes(label = gene_symbol),
      colour = "black",
      size = 6 / .pt,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.colour = "grey30",
      segment.size = 0.3,
      force = 10,
      force_pull = 0.05,
      max.time = 3,
      max.iter = 20000,
      max.overlaps = Inf,
      seed = 1)
  } else {
    geom_text(
      data = top_candidate_genes,
      aes(label = gene_symbol),
      colour = "black",
      size = 6 / .pt,
      hjust = -0.1,
      vjust = 0.5)
  }

  spearman_test <- suppressWarnings(cor.test(
    gene_counts$gene_length_bp,
    gene_counts$n_mutations,
    method = "spearman",
    exact = FALSE))
  spearman_rho_label <- paste0(
    "rho = ",
    formatC(unname(spearman_test$estimate), format = "f", digits = 2))

  ggplot(
    gene_counts,
    aes(
      x = gene_length_bp / 1000,
      y = n_mutations,
      colour = gene_class)) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      formula = y ~ x,
      fullrange = TRUE,
      se = FALSE,
      colour = "grey30",
      linewidth = 0.35) +
    geom_point(
      data = filter(gene_counts, !is_candidate),
      size = 0.35,
      alpha = 0.4) +
    geom_point(
      data = filter(gene_counts, is_candidate),
      size = 0.5,
      alpha = 1) +
    label_layer +
    annotate(
      "text",
      x = 19,
      y = 9.5,
      label = spearman_rho_label,
      hjust = 1,
      vjust = 1,
      size = 8 / .pt,
      colour = "black") +
    annotate(
      "text",
      x = 19,
      y = 8.5,
      label = "P < 2.2e-16",
      hjust = 1,
      vjust = 1,
      size = 8 / .pt,
      colour = "black") +
    scale_colour_manual(
      values = c(
        "Other genes" = "grey70",
        col_palette$candidate),
      name = NULL) +
    scale_x_continuous(
      breaks = seq(0, 20, by = 10),
      labels = scales::label_number(accuracy = 1, big.mark = ""),
      expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(
      breaks = seq(0, 15, by = 2),
      labels = scales::label_number(accuracy = 1, big.mark = ""),
      expand = expansion(mult = c(0, 0))) +
    coord_cartesian(
      xlim = c(0, 20),
      ylim = c(0, 15)) +
    labs(
      x = "Gene length (Kb)",
      y = "# mutations per gene") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.position = "none")
}

split_human_symbols <- function(orthologue_tbl) {
  orthologue_tbl %>%
    dplyr::select(gene_symbol, human_symbols) %>%
    mutate(human_symbols = coalesce(human_symbols, "")) %>%
    separate_rows(human_symbols, sep = ";") %>%
    transmute(
      gene_symbol,
      human_symbol = trimws(human_symbols)) %>%
    filter(human_symbol != "")
}

query_babelgene <- function(gene_symbols, min_support) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    return(tibble::tibble(gene_symbol = character(), human_symbol = character()))
  }

  gene_symbols <- sort(unique(na.omit(gene_symbols)))
  if (length(gene_symbols) == 0) {
    return(tibble::tibble(gene_symbol = character(), human_symbol = character()))
  }

  ortholog_args <- list(genes = gene_symbols, species = "mouse", human = FALSE)
  ortholog_formals <- names(formals(babelgene::orthologs))

  if ("min_support" %in% ortholog_formals) {
    ortholog_args$min_support <- min_support
  }

  if ("top" %in% ortholog_formals) {
    ortholog_args$top <- FALSE
  }

  tryCatch(
    do.call(babelgene::orthologs, ortholog_args) %>%
      as_tibble() %>%
      filter(!is.na(symbol), !is.na(human_symbol), human_symbol != "") %>%
      transmute(gene_symbol = symbol, human_symbol) %>%
      distinct(),
    error = function(err) {
      warning("babelgene orthologue lookup failed: ", conditionMessage(err))
      tibble::tibble(gene_symbol = character(), human_symbol = character())
    })
}

build_human_orthologues <- function(gene_symbols) {
  manual_hits <- manual_human_orthologues %>%
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
      tibble::tibble(gene_symbol = sort(unique(na.omit(gene_symbols)))),
      by = "gene_symbol") %>%
    group_by(gene_symbol) %>%
    summarise(
      human_symbols = paste(sort(unique(na.omit(human_symbol))), collapse = ";"),
      .groups = "drop") %>%
    mutate(
      human_symbols = coalesce(human_symbols, ""),
      has_human_orthologue = human_symbols != "")
}

load_cbio_mutated_gene_universe <- function(file_name = cbio_mutated_genes_file) {
  readr::read_tsv(file_name, show_col_types = FALSE) %>%
    transmute(human_symbol = as.character(Gene)) %>%
    filter(!is.na(human_symbol), human_symbol != "") %>%
    distinct()
}

build_mouse_cbio_mutated_gene_overlap <- function(
  mutation_tbl,
  cbio_gene_universe,
  treatment_filter = "DT") {

  mouse_mutated_genes <- mutation_tbl %>%
    filter(
      treatment == treatment_filter,
      !is.na(SYMBOL),
      SYMBOL != "") %>%
    transmute(gene_symbol = as.character(SYMBOL)) %>%
    distinct() %>%
    arrange(gene_symbol)

  orthologues <- build_human_orthologues(mouse_mutated_genes$gene_symbol) %>%
    mutate(
      human_ortholog = if_else(
        has_human_orthologue,
        str_split_fixed(human_symbols, ";", 2)[, 1],
        NA_character_))

  mouse_mutated_genes %>%
    left_join(
      orthologues %>%
        dplyr::select(gene_symbol, human_ortholog),
      by = "gene_symbol") %>%
    mutate(
      appears_in_cbioportal =
        !is.na(human_ortholog) &
        human_ortholog %in% cbio_gene_universe$human_symbol,
      appears_in_cbioportal_mutated_gene_universe = appears_in_cbioportal)
}

build_cbio_overlap_bar_summary <- function(overlap_tbl) {
  overlap_levels <- c(
    "Mutated in human cohorts",
    "Non-mutated in human cohorts")

  summary_tbl <- overlap_tbl %>%
    mutate(
      cbio_overlap_group = if_else(
        appears_in_cbioportal,
        "Mutated in human cohorts",
        "Non-mutated in human cohorts"),
      has_human_ortholog = !is.na(human_ortholog)) %>%
    group_by(cbio_overlap_group) %>%
    summarise(
      n_genes = n_distinct(gene_symbol),
      n_genes_without_human_ortholog = n_distinct(
        gene_symbol[!has_human_ortholog]),
      .groups = "drop") %>%
    right_join(
      tibble::tibble(cbio_overlap_group = overlap_levels),
      by = "cbio_overlap_group") %>%
    mutate(
      n_genes = replace_na(n_genes, 0L),
      n_genes_without_human_ortholog =
        replace_na(n_genes_without_human_ortholog, 0L),
      cbio_overlap_group = factor(
        cbio_overlap_group,
        levels = overlap_levels),
      x_position = as.numeric(cbio_overlap_group)) %>%
    arrange(cbio_overlap_group)

  summary_tbl %>%
    transmute(
      cbio_overlap_group,
      overlap_segment = "Genes with human ortholog",
      n_genes = n_genes - n_genes_without_human_ortholog,
      n_genes_total = n_genes,
      n_genes_without_human_ortholog,
      x_position) %>%
    bind_rows(
      summary_tbl %>%
        transmute(
          cbio_overlap_group,
          overlap_segment = "Genes without human ortholog",
          n_genes = n_genes_without_human_ortholog,
          n_genes_total = n_genes,
          n_genes_without_human_ortholog,
          x_position)) %>%
    filter(n_genes > 0) %>%
    mutate(
      overlap_segment = factor(
        overlap_segment,
        levels = c(
          "Genes with human ortholog",
          "Genes without human ortholog")))
}

plot_cbio_overlap_bar_summary <- function(summary_tbl) {
  if (!requireNamespace("ggpattern", quietly = TRUE)) {
    stop(
      "ggpattern is required for striped stacked bars. ",
      "Install ggpattern before plotting Supp4_mouse_cbio_mutated_gene_overlap_bar.",
      call. = FALSE)
  }

  x_labels <- summary_tbl %>%
    distinct(cbio_overlap_group, x_position) %>%
    arrange(x_position)

  ggplot(
    summary_tbl,
    aes(
      x = x_position,
      y = n_genes,
      fill = cbio_overlap_group,
      pattern = overlap_segment)) +
    ggpattern::geom_col_pattern(
      width = 0.7,
      colour = "black",
      linewidth = 0.15,
      pattern_fill = "black",
      pattern_colour = "black",
      pattern_angle = 45,
      pattern_density = 0.25,
      pattern_spacing = 0.05,
      pattern_key_scale_factor = 0.6) +
    scale_fill_manual(
      values = c(
        "Mutated in human cohorts" = "grey70",
        "Non-mutated in human cohorts" = "grey70")) +
    ggpattern::scale_pattern_manual(
      values = c(
        "Genes with human ortholog" = "none",
        "Genes without human ortholog" = "stripe"),
      guide = "none") +
    scale_x_continuous(
      breaks = x_labels$x_position,
      labels = as.character(x_labels$cbio_overlap_group),
      expand = expansion(mult = c(0.2, 0.2))) +
    scale_y_continuous(
      labels = my_label,
      expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "",
      y = "# genes") +
    my_theme +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(hjust = 0.5))
}

#-------------------------------------------------------------------------------
### Shared mutated gene UpSet plot ###
#-------------------------------------------------------------------------------
tes_presence <- prepare_mutations(
  raw_tbl = tes_unique,
  outliers = tes_outliers,
  dataset_label = "TES")

wes_presence <- prepare_mutations(
  raw_tbl = wes_unique,
  outliers = wes_outliers,
  dataset_label = "WES",
  use_last_cycle = TRUE)

shared_gene_base <- bind_rows(tes_presence, wes_presence) %>%
  build_nonsynonymous_base()

validate_outlier_exclusion(
  shared_gene_base,
  tes_outliers,
  "TES",
  "shared_gene_base")
validate_outlier_exclusion(
  shared_gene_base,
  wes_outliers,
  "WES",
  "shared_gene_base")

#-------------------------------------------------------------------------------
### Ras mutation counts ###
#-------------------------------------------------------------------------------
ras_mutation_base <- bind_rows(
  prepare_mutations(
    raw_tbl = tes_unique,
    outliers = tes_outliers,
    dataset_label = "TES",
    min_gt_af_filter = NULL),
  prepare_mutations(
    raw_tbl = wes_unique,
    outliers = wes_outliers,
    dataset_label = "WES",
    use_last_cycle = TRUE,
    min_gt_af_filter = NULL)) %>%
  build_nonsynonymous_base() %>%
  filter(
    SYMBOL %in% ras_gene_levels,
    tissue %in% ras_tissue_levels)

validate_outlier_exclusion(
  ras_mutation_base,
  tes_outliers,
  "TES",
  "ras_mutation_base")
validate_outlier_exclusion(
  ras_mutation_base,
  wes_outliers,
  "WES",
  "ras_mutation_base")

# Rows without parseable protein changes cannot be assigned a mutation label.
ras_mutation_map_data <- build_ras_mutation_map_data(ras_mutation_base)

Supp4_Ras_mutation_map <- plot_ras_mutation_map(ras_mutation_map_data)

export_plot_data(
  data = ras_mutation_map_data,
  file_name = file.path(supp4_dir, "Supp4_Ras_mutation_map"),
  width = 30,
  height = 30,
  panel_width = c(16.2, 5.4, 5.4),
  panel_height = 29,
  tick_size = 1,
  cols = c(
    "gene_symbol",
    "codon_position",
    "amino_acid_label",
    "mutation",
    "morphology",
    "n_samples",
    "n_mutations",
    "tissues",
    "datasets",
    "sample_names",
    "mutation_ids",
    "protein_positions",
    "amino_acids",
    "consequences",
    "mean_gt_AF",
    "max_gt_AF"))
Supp4_Ras_mutation_map

shared_gene_groups <- make_shared_gene_groups(
  mutations = shared_gene_base,
  group_levels = shared_gene_group_levels)

Supp4_shared_mutated_genes_upset <- save_shared_gene_upset(
  shared_groups = shared_gene_groups,
  group_levels = shared_gene_group_levels,
  file_stub = file.path(supp4_dir, "Supp4_shared_mutated_genes_upset"))

Supp4_shared_mutated_genes_upset

#-------------------------------------------------------------------------------
### WES dN/dS positive selection ###
#-------------------------------------------------------------------------------
wes_positive_dnds <- prep_positive_dnds(wes_dnds)

Supp4_WES_positive_dnds <- plot_positive_dnds(wes_positive_dnds)
export_plot_data(
  data = wes_positive_dnds,
  file_name = file.path(supp4_dir, "Supp4_WES_positive_dnds"),
  cols = c(
    "category",
    "gene",
    "dnds_class",
    "omega",
    "q",
    "sig",
    "positive_selection",
    "gene_score"))
Supp4_WES_positive_dnds

#-------------------------------------------------------------------------------
### Candidate gene mutation type composition ###
#-------------------------------------------------------------------------------
candidate_mutation_type_counts <- build_candidate_mutation_type_counts(
  mutation_tbl = shared_gene_base,
  candidate_matrix_tbl = figure4_candidate_genes)

Supp4_candidate_mutation_types <- plot_candidate_mutation_type_counts(
  candidate_mutation_type_counts)
export_candidate_mutation_type_plot(
  plot = Supp4_candidate_mutation_types,
  count_obj = candidate_mutation_type_counts,
  file_name = file.path(supp4_dir, "Supp4_candidate_mutation_types"),
  cols = c(
    "gene_symbol",
    "row_label",
    "gene_label",
    "candidate_group",
    "lollipop_order",
    "mutation_type",
    "n_mutations",
    "total_mutations"))
readr::write_csv(
  candidate_mutation_type_counts$mutations,
  file.path(supp4_dir, "Supp4_candidate_mutation_types_mutations_data.csv"))
Supp4_candidate_mutation_types

#-------------------------------------------------------------------------------
### Candidate gene mutation type random-gene null ###
#-------------------------------------------------------------------------------
candidate_mutation_type_null <- run_candidate_mutation_type_null(
  mutation_tbl = shared_gene_base,
  observed_counts = candidate_mutation_type_counts$counts)

Supp4_candidate_mutation_type_random_null <- plot_candidate_mutation_type_null(
  candidate_mutation_type_null)

readr::write_csv(
  candidate_mutation_type_null$long %>%
    select(any_of(c(
    "iteration",
    "n_sampled_genes",
    "mutation_type",
    "null_count"))),
  file.path(
    supp4_dir,
    "Supp4_candidate_mutation_type_random_null_data.csv"))
ggsave(
  filename = file.path(
    supp4_dir,
    "Supp4_candidate_mutation_type_random_null.pdf"),
  plot = Supp4_candidate_mutation_type_random_null,
  width = 60,
  height = 60,
  units = "mm")
readr::write_csv(
  candidate_mutation_type_null$summary,
  file.path(
    supp4_dir,
    "Supp4_candidate_mutation_type_random_null_summary_data.csv"))
Supp4_candidate_mutation_type_random_null

#-------------------------------------------------------------------------------
### Candidate genes mutated per sample ###
#-------------------------------------------------------------------------------
candidate_gene_sample_counts <- build_candidate_gene_sample_counts(
  mutation_tbl = shared_gene_base,
  candidate_matrix_tbl = figure4_candidate_genes)

candidate_gene_sample_counts_nonzero <- candidate_gene_sample_counts %>%
  filter(n_candidate_genes_mutated > 0) %>%
  mutate(
    category = recode(
      as.character(category),
      "Visually normal" = "Macr. normal"),
    category = factor(
      category,
      levels = c("Macr. normal", "SCC", "Papilloma"))) %>%
  arrange(category, n_candidate_genes_mutated, sample_name) %>%
  mutate(
    sample_order = row_number(),
    sample_uid_ordered = factor(
      as.character(sample_uid),
      levels = as.character(sample_uid)))

Supp4_candidate_genes_per_sample_nonzero <- plot_candidate_gene_sample_counts(
  candidate_gene_sample_counts_nonzero,
  x_title = "Samples with at least one candidate gene mutated",
  y_title = "# of unique candidate genes mutated",
  y_breaks = c(0, 9, 18),
  legend_key_size = grid::unit(0.5, "lines"))

export_plot_data(
  data = candidate_gene_sample_counts_nonzero,
  file_name = file.path(supp4_dir, "Supp4_candidate_genes_per_sample_nonzero"),
  cols = c(
    "sample_order",
    "sample_uid",
    "sample_name",
    "dataset",
    "tissue",
    "category",
    "n_candidate_genes_mutated"))
Supp4_candidate_genes_per_sample_nonzero

Supp4_candidate_genes_per_sample_nonzero_violin <-
  plot_candidate_gene_sample_counts_violin(
    candidate_gene_sample_counts_nonzero)

export_plot_data(
  data = candidate_gene_sample_counts_nonzero,
  file_name = file.path(
    supp4_dir,
    "Supp4_candidate_genes_per_sample_nonzero_violin"),
  panel_width = 29,
  panel_height = 29,
  tick_size = 1,
  cols = c(
    "sample_order",
    "sample_uid",
    "sample_name",
    "dataset",
    "tissue",
    "category",
    "n_candidate_genes_mutated"))
Supp4_candidate_genes_per_sample_nonzero_violin

#-------------------------------------------------------------------------------
### Gene length versus mutation count ###
#-------------------------------------------------------------------------------
refcds_gene_lengths <- load_refcds_gene_lengths(refcds_path)

gene_length_mutation_counts <- build_gene_length_mutation_counts(
  mutation_tbl = shared_gene_base,
  candidate_matrix_tbl = figure4_candidate_genes,
  gene_lengths = refcds_gene_lengths)

Supp4_gene_length_mutation_counts <- plot_gene_length_mutation_counts(
  gene_length_mutation_counts)

export_plot_data(
  data = gene_length_mutation_counts,
  file_name = file.path(supp4_dir, "Supp4_gene_length_mutation_counts"),
  cols = c(
    "gene_symbol",
    "gene_id",
    "protein_id",
    "gene_length_bp",
    "n_mutations",
    "is_candidate",
    "candidate_group",
    "gene_class"))
Supp4_gene_length_mutation_counts

#-------------------------------------------------------------------------------
### Mouse mutated genes in cBioPortal mutated-gene universe ###
#-------------------------------------------------------------------------------
cbio_mutated_gene_universe <- load_cbio_mutated_gene_universe(
  cbio_mutated_genes_file)

Supp4_mouse_cbio_mutated_gene_overlap <- build_mouse_cbio_mutated_gene_overlap(
  mutation_tbl = shared_gene_base,
  cbio_gene_universe = cbio_mutated_gene_universe)

readr::write_csv(
  Supp4_mouse_cbio_mutated_gene_overlap,
  file.path(supp4_dir, "Supp4_mouse_cbio_mutated_gene_overlap_data.csv"))
Supp4_mouse_cbio_mutated_gene_overlap

Supp4_mouse_cbio_mutated_gene_overlap_summary <-
  build_cbio_overlap_bar_summary(Supp4_mouse_cbio_mutated_gene_overlap)

Supp4_mouse_cbio_mutated_gene_overlap_bar <-
  plot_cbio_overlap_bar_summary(Supp4_mouse_cbio_mutated_gene_overlap_summary)

export_plot_data(
  data = Supp4_mouse_cbio_mutated_gene_overlap_summary,
  file_name = file.path(supp4_dir, "Supp4_mouse_cbio_mutated_gene_overlap_bar"),
  cols = c(
    "cbio_overlap_group",
    "overlap_segment",
    "n_genes",
    "n_genes_total",
    "n_genes_without_human_ortholog"))
Supp4_mouse_cbio_mutated_gene_overlap_bar

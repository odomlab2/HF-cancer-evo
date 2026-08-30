################################################################################
# @Project - TES Hair Follicle
# @Date - 05/01/2026
# @Author - Yoav Avi-Guy
# @Description - This file provides the colour palettes and plotting parameters
################################################################################

#-------------------------------------------------------------------------------
# Libraries
library(ggplot2)
library(readr)

# Publication writers are redirected to an explicit isolated output root.
if (!exists("publication_output_path", mode = "function")) {
  plotting_root <- Sys.getenv("HF_SCC_ROOT", unset = "")
  if (!nzchar(plotting_root)) {
    plotting_root <- normalizePath(getwd(), mustWork = TRUE)
    while (!dir.exists(file.path(plotting_root, ".git"))) {
      parent <- dirname(plotting_root)
      if (identical(parent, plotting_root)) {
        stop("Set HF_SCC_ROOT to the repository root.")
      }
      plotting_root <- parent
    }
  }
  source(file.path(plotting_root, "scripts", "publication_output.R"))
}

#-------------------------------------------------------------------------------
# ggplot theme
my_theme <- theme_minimal() +
    theme(
    text = element_text(family = "sans", colour = "black", size = 8),
    axis.text.x = element_text(color = "black", margin = margin(t = 2)),
    axis.text.y = element_text(color = "black", margin = margin(r = 2)),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", element_line(1)),
    axis.ticks.length = grid::unit(1, "mm"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5))

#-------------------------------------------------------------------------------
# Exponential axis labelling 
my_label <- function(x) {
  vapply(x, function(y) {
    if (is.na(y)) return(NA_character_)
    if (y == 0)  return("0")

    e <- floor(log10(abs(y)))
    m <- y / (10^e)

    m_str <- formatC(m, format = "g", digits = 3)   # mantissa precision
    e_str <- sub("^0+", "", as.character(abs(e)))   # drop leading zeros

    paste0(m_str, "e", ifelse(e >= 0, "+", "-"), e_str)

  }, character(1))
}

#-------------------------------------------------------------------------------
# Plot knitting parameters
# knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)

#-------------------------------------------------------------------------------
# Colour palettes
col_palette <- list(
  core = c(
    "#de4243",
    "#bf2525",
    "#f1631b",
    "#ff8848",
    "#feb334",
    "#fed933",
    "#ffe757",
    "#75d478",
    "#1a8551",
    "#025245"
  ),
  treatment = c(
    Acetone = "#F28B50",
    DT = "#732642"
  ),
  tissue = c(
    Skin = "#518E95",
    "Hair follicle" = "#BF6B63" 
  ),
  category = c(
    Acetone = "#cbcecd",
    "Visually normal" = "#F2B27A",
    Papilloma = "#b5ae7a",
    SCC = "#84b6e3"
  ),
  mutation_type = c(
    SNV = "#293f63", 
    deletion = "#95bfb8", 
    insertion = "#027333", 
    substitution = "#d9cd29"
  ),
  mutation_signature = c(
    # "C>A" = "#E8CCCCD8",
    # "C>G" = "#EBC0CC",
    # "C>T" = "#CE7A92",
    # "T>A" = "#143259",
    # "T>C" = "#8FB6D9",
    # "T>G" = "#CEC3B8FF"
    "C>A" = "#E5C8C8",
    "C>G" = "#DDA9B8",
    "C>T" = "#BC667E",
    "T>A" = "#1A3A63",
    "T>C" = "#7FA9D0",
    "T>G" = "#a9b5c4ff"
  ),
  dnds_class = c(
    "Missense" = "#10755F",
    "Nonsense" = "#699F4F",
    "Indel/DBS" = "#ADCE71"
  ),
  adjacency_class = c(
    "Adjacent" = "#BB377D",
    "Distal" = "#F084BD"
  ),
  mouse = c(
    "21" = "#5A878C",
    "19" = "#848C46"
  ),
  time = c(
    "Week 8" = "#bee5ee",
    "Week 12" = "#85d3e2",
    "Week 14" = "#01b8d8",
    "Week 17" = "#0e7cb8",
    "Week 19" = "#1b2454"
  ),
  candidate = c(
    "Initiation (I)" = "#2ca091",
    "I+P" = "#81ae81",
    "Predisposition" = "#e4bf6e"
  )
)

#-------------------------------------------------------------------------------
# Ordering variables for plotting
plot_levels <- function(df) {
  lvl <- list(
    tissue    = c("Hair follicle", "Skin"),
    treatment = c("Acetone", "DT"),
    category  = c("Acetone", "SCC", "Papilloma", "Visually normal"),
    time      = c("Week 8", "Week 14", "Week 17", "Week 19"),
    mutation  = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  )

  for (nm in intersect(names(lvl), names(df))) {
    df[[nm]] <- factor(df[[nm]], levels = lvl[[nm]], ordered = TRUE)
  }
  df
}

#-------------------------------------------------------------------------------
# Dataframe export
export_plot_data <- function(
  file_name,
  data,
  cols = NULL,
  width = NULL,
  height = NULL,
  panel_width = 29,
  panel_height = 29,
  tick_size = 1) {

  out <- if (is.null(cols)) data else select(data, any_of(cols))
  plot <- last_plot()
  plot_layout <- ggplot_build(plot)$layout$layout
  n_panel_cols <- max(plot_layout$COL)
  n_panel_rows <- max(plot_layout$ROW)

  plot <- plot +
    ggh4x::force_panelsizes(
      cols = grid::unit(panel_width, "mm"),
      rows = grid::unit(panel_height, "mm")
    )

  if (is.null(width)) width <- n_panel_cols * (panel_width + tick_size) + 30
  if (is.null(height)) height <- n_panel_rows * (panel_height + tick_size) + 30

  write_csv(out, paste0(file_name, "_data.csv"))

  ggsave(
    filename = paste0(file_name, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = "mm"
  )

  invisible(out)
}

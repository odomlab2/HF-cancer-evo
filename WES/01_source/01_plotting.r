################################################################################
# @Project - TES Hair Follicle
# @Date - 05/01/2026
# @Author - Yoav Avi-Guy
# @Description - This file provides the colour palettes and plotting parameters
################################################################################

#-------------------------------------------------------------------------------
# Libraries
library(ggplot2)

#-------------------------------------------------------------------------------
# ggplot theme
my_theme <- theme_minimal() +
    theme(
    text = element_text(family = "sans", colour = "black", size = 18),
    axis.text.x = element_text( color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", element_line(1)),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1))

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
    Acetone = "#90e0ef",
    "Visually normal" = "#00b4d8",
    Papilloma = "#0077b6",
    SCC = "#03045e"
  ),
  mutation_type = c(
    SNV = "#293F63", 
    deletion = "#95BFB8", 
    insertion = "#027333", 
    substitution = "#D9CD29"
  ),
  mutation_signature = c(
    # "C>A" = "#F7E3E3",
    # "C>G" = "#EBC0CC",
    # "C>T" = "#CE7A92",
    # "T>A" = "#143259",
    # "T>C" = "#8FB6D9",
    # "T>G" = "#F2E7DC"
    "C>A" = "#E5C8C8",
    "C>G" = "#DDA9B8",
    "C>T" = "#BC667E",
    "T>A" = "#1A3A63",
    "T>C" = "#7FA9D0",
    "T>G" = "#D1C2B4"
  ),
  dnds_class = c(
    "Missense" = "#10755F",
    "Nonsense" = "#699F4F",
    "Indel/DBS" = "#ADCE71"
  ),
  adjacency_class = c(
    "Adjacent" = "#BF300F",
    "Distal" = "#A49A35"
  ),
  mouse = c(
    "3" = "#A8A384",
    "4" = "#5A878C",
    "5" = "#848C46",
    "9" = "#87624C"
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


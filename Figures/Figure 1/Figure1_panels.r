# 10/04/2026
# Author: Yoav Avi-Guy
# Purpose: Figure 1 panels - Longitudinal and spatial sampling enables high-
# resolution mutational profiling of hair follicles during skin carcinogenesis
#-------------------------------------------------------------------------------
# Libraries
library(dplyr)
repo_root <- normalizePath(Sys.getenv("HF_SCC_ROOT", unset = "."), mustWork = TRUE)
output_dir <- Sys.getenv("HF_SCC_OUTPUT_ROOT")
if (!nzchar(output_dir)) stop("Set HF_SCC_OUTPUT_ROOT to an external output directory.")
output_dir <- file.path(output_dir, "Figures", "Figure 1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Colours and theme
source(file.path(repo_root, "TES", "01_source", "01_plotting.r"))

# Loading files
TES_DNA_conc <- read.csv(file.path(repo_root, "TES", "02_data", "01_metadata", "TES_dna_concentrations.csv"))
TES_DNA_conc <- TES_DNA_conc %>% filter(tissue != "liver")

#-------------------------------------------------------------------------------
### DNA yield per biopsy ###
#-------------------------------------------------------------------------------
# Capping plot at 50 ng
TES_DNA_df <- TES_DNA_conc %>%
  mutate(
    y_plot_DNA = if_else(
      tissue == "Hair follicle",
      pmin(input_conc, 50),
      input_conc),
    above_cap_reads = tissue == 
      "Hair follicle" & input_conc > 50)

TES_DNA_df <- plot_levels(TES_DNA_df)

TES_DNA_yield <- ggplot(
  TES_DNA_df %>% filter(tissue == "Hair follicle"),
  aes(x = time, y = y_plot_DNA, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) +
  my_theme +
  scale_fill_manual(values = col_palette$tissue) +
  labs(y = "Input DNA (ng)", x = "") +
  theme(legend.position = "") +
  geom_jitter(
    data = ~ dplyr::filter(.x, !above_cap_reads),
    position = position_jitterdodge(
      jitter.width = 0.35,
      dodge.width = 0.75),
    size = 0.5,
    alpha = 0.6,
    shape = 16) +
  geom_jitter(
    data = ~ dplyr::filter(.x, above_cap_reads),
    position = position_jitterdodge(
      jitter.width = 0.35,
      dodge.width = 0.75),
    size = 1.8,
    alpha = 0.9,
    shape = 17)
export_plot_data(
  data = TES_DNA_df %>% filter(tissue == "Hair follicle"),
  file_name = file.path(output_dir, "Figure1_TES_DNA_yield"),
  cols = c(
    "sample_name",
    "tissue",
    "time",
    "input_conc",
    "y_plot_DNA",
    "above_cap_reads"))
TES_DNA_yield

# =========================================================
# Papilloma over time + Kaplan-Meier plot in R
# Reproduces what I did:
# 1) line plot of papilloma count per mouse
# 2) Kaplan-Meier plot of tumor-free mice over time
# 3) saves both plots into one PDF
# 4) first plot x-axis shown from week 2 to 20
# =========================================================

# Install packages once if needed:
# install.packages(c("ggplot2", "tidyr", "dplyr", "survival", "survminer"))

library(ggplot2)
library(tidyr)
library(dplyr)
library(survival)
library(survminer)
library(gridExtra)

# ---------------------------
# 1. Enter your data
# ---------------------------
df <- data.frame(
  Weeks = 1:20,
  `14M` = c(0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,4),
  `16M` = c(0,0,0,0,0,0,0,0,0,0,1,1,2,2,4,6,6,6,7,6),
  `23M` = c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,2,2,5,4,4),
  `24M` = c(0,0,0,0,0,0,0,0,0,0,1,2,2,3,3,3,4,5,4,4),
  `21M` = c(0,0,0,0,0,0,0,1,1,3,3,4,6,6,7,7,6,6,6,6),
  check.names = FALSE
)

# ---------------------------
# 2. Convert to long format for plotting
# ---------------------------
df_long <- df %>%
  pivot_longer(
    cols = -Weeks,
    names_to = "Mouse",
    values_to = "Papilloma_Count"
  )

# ---------------------------
# 3. Plot papilloma count over time
# ---------------------------
p1 <- ggplot(df_long, aes(x = Weeks, y = Papilloma_Count, color = Mouse)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # scale_x_continuous(limits = c(2, 20), breaks = 2:20) +
  labs(
    title = "Papilloma growth per mouse",
    x = "Weeks",
    y = "Papilloma count"
  ) +
  my_theme

# ---------------------------
# 4. Build Kaplan-Meier input
# Event = first week with papilloma count > 0
# ---------------------------
km_df <- df_long %>%
  group_by(Mouse) %>%
  summarise(
    event_time = ifelse(
      any(Papilloma_Count > 0),
      min(Weeks[Papilloma_Count > 0]),
      max(Weeks)),
    status = ifelse(any(Papilloma_Count > 0), 1, 0),
    .groups = "drop"
  )

print(km_df)

# ---------------------------
# 5. Kaplan-Meier fit
# ---------------------------
fit <- survfit(Surv(event_time, status) ~ 1, data = km_df)

# KM plot
p2 <- ggsurvplot(
  fit,
  data = km_df,
  conf.int = FALSE,
  censor = FALSE,
  risk.table = FALSE,
  xlab = "Weeks",
  ylab = "Tumor-free survival",
  title = "Kaplan-Meier curve",
  ggtheme = my_theme
)

# ---------------------------
# 6. Save both plots into one PDF
# ---------------------------
pdf(file.path(output_dir, "Figure1_tumour_plots.pdf"))

grid.arrange(
  p1,
  p2$plot,
  ncol = 1
)

dev.off()

# =========================================================
# Optional: save each plot separately too
# =========================================================

# ggsave("papilloma_growth_updated.pdf", plot = p1, width = 8, height = 5)

# ggsave("kaplan_meier.pdf", plot = p2$plot, width = 8, height = 5)

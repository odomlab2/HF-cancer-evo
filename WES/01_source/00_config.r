################################################################################
# @Project - TES Hair Follicle
# @Description - This file is responsible for the loading the commonly used libraries, paths, and functions
################################################################################

#-------------------------------------------------------------------------------
# repository and publication-output contract
find_repo_root <- function(path = getwd()) {
  configured <- Sys.getenv("HF_SCC_ROOT", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, ".git"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Set HF_SCC_ROOT to the repository root.")
    path <- parent
  }
}
HF_SCC_ROOT <- find_repo_root()
source(file.path(HF_SCC_ROOT, "scripts", "publication_output.R"))
setwd(file.path(HF_SCC_ROOT, "WES"))

#-------------------------------------------------------------------------------
# libraries
library(dplyr)
library(tidyr)

#-------------------------------------------------------------------------------
# PATHS
BASE_PATH <- file.path(HF_SCC_ROOT, "WES")
SOURCE_PATH <- file.path(BASE_PATH, "01_source")
INPUT_PATH <- file.path(BASE_PATH, "02_data/00_raw")
METADATA_PATH <- file.path(BASE_PATH, "02_data/01_metadata")
PROCESSED_PATH <- file.path(BASE_PATH, "02_data/02_processed")
SESSION_INFO_PATH <- file.path(BASE_PATH, "03_analysis/session_info")
PLOTS_PATH <- file.path(BASE_PATH, "04_results/01_plots")

#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, warn = 1, digits = 17)
Sys.setenv(TZ = "UTC")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript scripts/public/reproduce_submission_snapshot.R public-output/NAME",
    call. = FALSE)
}
if (!grepl("^public-output/[^/]+", args[[1]]) || grepl("(^|/)\\.\\.(/|$)", args[[1]])) {
  stop("Output must be a new directory below public-output/.", call. = FALSE)
}
root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "scripts/public/common.R"))) {
  stop("Run this command from the repository root.", call. = FALSE)
}
output <- file.path(root, args[[1]])
if (dir.exists(output)) stop("Refusing to overwrite: ", args[[1]], call. = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

source(file.path(root, "scripts/public/common.R"))
source(file.path(root, "scripts/public/figure2a.R"))
source(file.path(root, "scripts/public/frozen_statistics.R"))
source(file.path(root, "scripts/public/figure5g.R"))
source(file.path(root, "scripts/public/submitted_panels.R"))

figure2a <- reproduce_figure2a(root)
frozen <- validate_frozen_statistics(root)
figure5g <- reproduce_figure5g(root)
panels <- validate_submitted_panels(root)
aliases <- read_publication_csv(root,
  "publication_inputs/figure5e_display_aliases.csv")

write_result(figure2a, output, "figure2a_kruskal_wallis.csv")
write_result(frozen$sharing, output, "shared_mutation_statistics.csv")
write_result(frozen$supp3c, output, "suppfigure3c_beta_binomial.csv")
write_result(figure5g, output, "figure5g_offline_permutation.csv")
write_result(panels$figure4f, output, "figure4f_submitted_values.csv")
write_result(panels$figure5c, output, "figure5c_exact_equality_proof.csv")
write_result(aliases, output, "figure5e_display_aliases.csv")

decisions <- data.frame(invariant = c("Figure 2A TES", "Figure 2A WES",
  "Supplementary Figure 3C", "overall HF/skin sharing",
  "tumour-skin sharing", "Figure 5G summary"),
  authoritative = c(format(figure2a$p_value[[1]], digits = 17),
    format(figure2a$p_value[[2]], digits = 17),
    "OR=1.8121921461362716; BH P=0.00063761766723724066",
    "214/8133; display=2.63%", "115/953; display=12.1%",
    "median=13.95; iterations=10000; P=9.999000099990002e-5"),
  status = "VALIDATED")
write_result(decisions, output, "publication_discrepancy_table.csv")
capture.output(sessionInfo(), file = file.path(output, "r_session_info.txt"))
writeLines(c("PASS: all submitted-manuscript snapshot assertions succeeded.",
  "No live external service was used.",
  "No generated RDS or PDF is required by this public validation path."),
  file.path(output, "execution_assertions.txt"))
message("PASS: manuscript snapshot reproduced in ", args[[1]])

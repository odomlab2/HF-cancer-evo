################################################################################
# @Project - TES Hair Follicle
# @Date - 05/01/2026
# @Author - Yoav Avi-Guy
# @Description - This file makes the targeted gene panel compatible with Mutect2
# gene list for further analysis
################################################################################

#-------------------------------------------------------------------------------
# SOURCE
source(file.path(Sys.getenv("HF_SCC_ROOT", unset = "."), "TES/01_source/00_config.r"))

#-------------------------------------------------------------------------------
# Loading work files
mutations_annotated <- readRDS(
  file.path(PROCESSED_PATH, "mutations_annotated.rds"))
gene_panel <- read.table(
  file.path(INPUT_PATH, "targeted_gene_panel.bed"))

#-------------------------------------------------------------------------------
# Correcting capitalisation of gene panel genes according to Mutect2
lookup <- mutations_annotated %>%
  distinct(SYMBOL) %>%
  mutate(V4_upper = toupper(SYMBOL))

gene_panel_caps <- gene_panel %>%
  mutate(V4_upper  = V4) %>%
  left_join(lookup, by = "V4_upper") %>%
  mutate(V4 = coalesce(SYMBOL, V4_upper)) %>% 
  select(-V4_upper, -SYMBOL)

#-------------------------------------------------------------------------------
# Finding genes that are in the panel but are not found in Mutect2
only_in_panel_genes <- setdiff(
  unique(gene_panel_caps$V4),
  unique(mutations_annotated$SYMBOL)
)

#-------------------------------------------------------------------------------
# Changing panel gene names to match the most recently assigned gene names and
# assigning column names
gene_panel_corrected <- gene_panel_caps %>%
  mutate(
    V4 = case_when(
      V4 == "VPRBP" ~ "Dcaf1",
      V4 == "GCN1L1" ~ "Gcn1",
      V4 == "GLTSCR1" ~ "Bicra",
      V4 == "HIST1H2BM" ~ "H2bc14",
      TRUE ~ V4
    )
  )

colnames(gene_panel_corrected) <- c("CHROM", "start", "end", "gene")

#-------------------------------------------------------------------------------
# saving annotated dataframe
saveRDS(gene_panel_corrected,
file.path(PROCESSED_PATH, "/gene_panel.rds"))

#-------------------------------------------------------------------------------
# session info
capture.output(sessionInfo(),
  file = file.path(SESSION_INFO_PATH, "06_gene_panel_session.txt"))

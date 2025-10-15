# ---- Setup ----
source("renv/activate.R")
source(snakemake@input[["functions"]])
source(snakemake@input[["settings"]])

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

# Load data --------------------------------------------------------------
df_long <- readRDS(snakemake@input[["df_long"]])
df_long_l2fc <- readRDS(snakemake@input[["df_long_l2fc"]])

# Cytosolic ribosomal proteins -------------------------------------------
print_section_header("Creating Ribosomal cytosolic proteins plots (Fig. 6)")

# Terms of interest
terms <- c("GO:0006413", "GO:0006414", "GO:0022625", "GO:0022627")
term_labels <- c(
  "GO:0006413" = "Init.\nfactors",
  "GO:0006414" = "Elong.\nfactors",
  "GO:0022627" = "Small ribo.\nsubunits",
  "GO:0022625" = "Large ribo.\nsubunits"
)

# Set order
term_order <- c(
  "GO:0006413",
  "GO:0006414",
  "GO:0022627",
  "GO:0022625"
)

# Retrieve plot and results
ribosomal <- plot_terms(terms, term_labels, term_order)

plot_ribo <- ribosomal[["plot"]] +
  coord_cartesian(ylim = c(-0.15, 0.6))

# Save plot
ggplot2::ggsave(snakemake@output[["fig_ribosomal_proteins"]], width = 200, height = 50, units = "mm", plot = plot_ribo)


# Mitochondrial proteins -------------------------------------------------
print_section_header("Creating Mitochondrial complex proteins plots (Fig. 6)")

# Complexes to iterate over
complexes <- c("CI", "CII", "CIII", "CIV", "CV")

# Retrieve plot and results
mito <- plot_complexes(complexes)

plot_mito <- mito[["plot"]] +
  coord_cartesian(ylim = c(-0.6, 0.3))

# Save plot
ggplot2::ggsave(snakemake@output[["fig_mitochondrial_proteins"]], width = 200, height = 50, units = "mm", plot = plot_mito)


# Structural GO-terms ----------------------------------------------------
print_section_header("Creating structural GO-term plots (Fig. 6)")

# Choose GO:terms to analyze
terms <- c("GO:0031012", "GO:0006936", "GO:0030239", "GO:0030036")
term_labels <- c(
  "GO:0031012" = "ECM",
  "GO:0006936" = "Muscle\ncontraction",
  "GO:0030239" = "Myofibril\nassembly",
  "GO:0030036" = "Actin\ncytosk.org."
)

# Set order
term_order <- c(
  "GO:0030036",
  "GO:0006936",
  "GO:0031012",
  "GO:0030239"
)

# Retrieve plot and results
structural <- plot_terms(terms, term_labels, term_order)

plot_structural <- structural[["plot"]] +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

# Save plot
ggplot2::ggsave(snakemake@output[["fig_structural_go_terms"]], width = 100, height = 50, units = "mm", plot = plot_structural)


# Metabolic GO-terms -----------------------------------------------------
print_section_header("Creating metabolic GO-term plots (Fig. 6)")

# Choose GO:terms to analyze
terms <- c("GO:0006099", "GO:0006635", "GO:0061621")
term_labels <- c(
  "GO:0006099" = "TCA",
  "GO:0006635" = "FA ox.",
  "GO:0061621" = "Gly-\ncolysis"
)

# Set order
term_order <- c(
  "GO:0006635",
  "GO:0061621",
  "GO:0006099"
)

# Retrieve plot and results
metabolic <- plot_terms(terms, term_labels, term_order)

plot_metabolic <- metabolic[["plot"]] +
  coord_cartesian(ylim = c(-0.35, 0.3))

# Save plot
ggsave(snakemake@output[["fig_metabolic_go_terms"]], width = 100, height = 50, units = "mm", plot = plot_metabolic)

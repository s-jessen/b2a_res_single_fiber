# ---- Setup ----
source("renv/activate.R")
source(snakemake@input[["functions"]])
source(snakemake@input[["settings"]])

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)})

# Load results  ----------------------------------------------------------
# Intervention group independent
results_main <- readRDS(snakemake@input[["results_main"]])
results_i <- readRDS(snakemake@input[["results_i"]])
results_ii <- readRDS(snakemake@input[["results_ii"]])
results_interaction <- readRDS(snakemake@input[["results_interaction"]])

# Terbutaline
results_ter_main <- readRDS(snakemake@input[["results_ter_main"]])
results_ter_i <- readRDS(snakemake@input[["results_ter_i"]])
results_ter_ii <- readRDS(snakemake@input[["results_ter_ii"]])
results_ter_interaction <- readRDS(snakemake@input[["results_ter_interaction"]])

# Resistance training
results_res_main <- readRDS(snakemake@input[["results_res_main"]])
results_res_i <- readRDS(snakemake@input[["results_res_i"]])
results_res_ii <- readRDS(snakemake@input[["results_res_ii"]])
results_res_interaction <- readRDS(snakemake@input[["results_res_interaction"]])


# Volcano plots ----------------------------------------------------------
print_section_header("Creating volcano plots (Fig. 3)")

# Independent, main
labels_main <- c(
  "S100A13",
  "MYBPH",
  "THY1",
  "RRAD",
  "AKR1C3",
  "ANKRD2",
  "VIM",
  "YBX1",
  "SERPINH1",
  "ART3",
  "FAM185A",
  "DHRSC7"
)

volcano_main <- volcano_plot(results_main, labels_main) +
  labs(title = "Main effect")

# Independent, type 1
labels_i <- c(
  "S100A13",
  "THY1",
  "TUBB2B",
  "VIM",
  "SULT1A3",
  "CYC1",
  "DNAJC11"
)

volcano_i <- volcano_plot(results_i, labels_i) +
  labs(title = "Type I")

# Independent, type 2
labels_ii <- c(
  "MYBPH",
  "S100A13",
  "ACTC1",
  "MYL2",
  "THY1",
  "AKR1C3",
  "AKR1C2",
  "SERPINH1",
  "YBX1",
  "SERPINB1",
  "ALDH1A2",
  "ART3"
)

volcano_ii <- volcano_plot(results_ii, labels_ii) +
  labs(title = "Type II")

# Independent, interaction
labels_interaction <- c(
  "MYH6",
  "MYL2",
  "ATP2A2",
  "PPP1R3D",
  "CREG1",
  "NIBAN1",
  "SERPINB3",
  "MAPRE3"
)

volcano_interaction <- volcano_plot(results_interaction, labels_interaction) +
  labs(title = "Interaction effect (Type I vs. Type II)")

# Res, main
labels_res_main <- c(
  "MUSTN1",
  "ANKRD2",
  "MYBPH",
  "THY1",
  "S100A13",
  "TMSB4X",
  "CSRP3",
  "KLHL40",
  "VIM"
)

volcano_res_main <- volcano_plot(results_res_main, labels_res_main) +
  labs(title = "Main effect")

# Res, type 1
labels_res_i <- c(
  "MUSTN1",
  "ANKRD2",
  "MYBPH",
  "S100A13"
)

volcano_res_i <- volcano_plot(results_res_i, labels_res_i) +
  labs(title = "Type I")

# Res, type 2
labels_res_ii <- c(
  "MYH6",
  "KLHL40",
  "DNAJA4",
  "MUSTN1",
  "S100A13",
  "THY1",
  "CSRP3",
  "KLHL41"
)

volcano_res_ii <- volcano_plot(results_res_ii, labels_res_ii) +
  labs(title = "Type II")

# Res, interaction
labels_res_interaction <- c(
  "SERPINB3",
  "PGM5",
  "MYH6",
  "MYH4",
  "TNNC1",
  "TNNI1",
  "MYL3",
  "MYL2"
)

volcano_res_interaction <- volcano_plot(results_res_interaction, labels_res_interaction) +
  labs(title = "Interaction effect (Type I vs. Type II)")

# Ter main
labels_ter_main <- c(
  "S100A13",
  "FTL",
  "PRKAA1",
  "AAMDC",
  "SORB51"
)

volcano_ter_main <- volcano_plot(results_ter_main, labels_ter_main) +
  labs(title = "Main effect")

# Ter type I
labels_ter_i <- c(
  "S100A13",
  "UNC45B",
  "FTL",
  "ATP2B1",
  "RPS19",
  "RPL3L",
  "DSP"
)

volcano_ter_i <- volcano_plot(results_ter_i, labels_ter_i) +
  labs(title = "Type I")

# Ter type II
labels_ter_ii <- c(
  "S100A13",
  "LCP1",
  "SOD3",
  "HIBADH",
  "PRKAA1",
  "OSTF1"
)

volcano_ter_ii <- volcano_plot(results_ter_ii, labels_ter_ii) +
  labs(title = "Type II")

# Ter interaction
labels_ter_interaction <- c(
  "HOMER2",
  "COLA1A2",
  "RPL3L"
)

volcano_ter_interaction <- volcano_plot(results_ter_interaction, labels_ter_interaction) +
  labs(title = "Interaction effect (Type I vs. Type II)")

# Collect volcanoes into single layout
all_volcanoes <- (volcano_main | volcano_i | volcano_ii | volcano_interaction) /
  (volcano_res_main | volcano_res_i | volcano_res_ii | volcano_res_interaction) /
  (volcano_ter_main | volcano_ter_i | volcano_ter_ii | volcano_ter_interaction)

ggplot2::ggsave(snakemake@output[["fig_volcanoes"]], width = 180, height = 150, units = "mm", plot = all_volcanoes)

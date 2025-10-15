source("renv/activate.R")
source(snakemake@input[["functions"]])
source(snakemake@input[["settings"]])

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(PhosR)
  library(SummarizedExperiment)
  library(eulerr)
  library(pcaMethods)
  library(ggrepel)
})

# Load required files
results_myh_ii_vs_i <- readRDS(snakemake@input[["results_myh_ii_vs_i"]])
df_long <- readRDS(snakemake@input[["df_long"]])
se_raw <- readRDS(snakemake@input[["se_raw"]])
se <- readRDS(snakemake@input[["se"]])
metadata <- readRDS(snakemake@input[["metadata"]])


# MYH distribution -------------------------------------------------------
print_section_header("Creating MYH distribution figure (Fig. 2)")

# Undo log2transformation
df_long_myh <- df_long %>%
  dplyr::select(protein, sample, abundance, intervention, fiber_type) %>%
  dplyr::filter(grepl("MYH", protein)) %>%
  dplyr::mutate(abundance = 2^abundance) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(relative_abundance = abundance / sum(abundance, na.rm = TRUE) * 100) %>%
  dplyr::mutate(protein = factor(protein, levels = c(
    "MYH1",
    "MYH2",
    "MYH3",
    "MYH4",
    "MYH6",
    "MYH7",
    "MYH8",
    "MYH9",
    "MYH11",
    "MYH13",
    "MYH14"
  )))

# Export dataset
readr::write_csv(df_long_myh, snakemake@output[["df_long_myh"]])

# Create palette
palette <- c(
  "MYH1" = "#12466e",
  "MYH2" = myh2_color,
  "MYH3" = "#23829f",
  "MYH4" = "#e9db53",
  "MYH6" = "#e99022",
  "MYH7" = myh7_color,
  "MYH8" = "#7b7b7b",
  "MYH9" = "#b4b4b4",
  "MYH11" = "#f1f1f1",
  "MYH13" = "#e9e9e9",
  "MYH14" = "#e1e1e1"
)

# Plot
myh_plot <- df_long_myh %>%
  ggplot(aes(x = sample, y = relative_abundance, fill = protein)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  facet_grid(~fiber_type,
    scales = "free_x",
    labeller = labeller(fiber_type = c("I" = "Type I pools", "II" = "Type II pools"))
  ) +
  # scale_fill_viridis_d(option = "turbo")+
  scale_fill_manual(values = palette) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Sample") +
  ylab("MYH expression (%)")

ggplot2::ggsave(
  filename = snakemake@output[["fig_myh_distribution"]],
  plot = myh_plot,
  width = 100,
  height = 50,
  units = "mm"
)


# Venn diagram -----------------------------------------------------------
print_section_header("Creating venn diagram (Fig. 2)")

# Subset for MYH7
myh_7 <- se_raw[, se_raw$fiber_type == "I"]

# Filter for 70% valid values
myh_7_filtered <- PhosR::selectGrps(myh_7, colData(myh_7)$time, 0.7, n = 1)

# Same for MYH2
myh_2 <- se_raw[, se_raw$fiber_type == "II"]
myh_2_filtered <- PhosR::selectGrps(myh_2, colData(myh_2)$time, 0.7, n = 1)

# Plot euler
venn <- plot(
  eulerr::euler(
    list(
      MYH7 = rownames(SummarizedExperiment::assay(myh_7_filtered)),
      MYH2 = rownames(SummarizedExperiment::assay(myh_2_filtered))
    )
  ),
  fills = c("white", "firebrick3"), quantities = list(type = "counts"), legend = TRUE
)

# Save plot
ggplot2::ggsave(snakemake@output[["fig_venn"]], width = 200, height = 200, units = "mm", plot = venn)


# Ranked protein intensities ---------------------------------------------
print_section_header("Creating ranked protein intensities plot (Fig. 2)")

# Undo log2transformation
df_long_raw <- df_long %>%
  dplyr::mutate(abundance = 2^abundance)

# Select labels
labels <- c(
  "ACTA1",
  "MB",
  "CKM",
  "MYH7",
  "MYH2",
  "KLHL41",
  "MYBPH",
  "AKR1C1",
  "AKR1C2",
  "AKR1C3",
  "S100A13",
  "COL1A1",
  "COL1A2"
)

# Compute log10 median intensity
df_long_log10 <- df_long_raw %>%
  dplyr::group_by(protein) %>%
  dplyr::summarize(median_abundance = median(abundance, na.rm = T)) %>%
  dplyr::mutate(log10median = log10(median_abundance)) %>%
  dplyr::arrange(desc(log10median)) %>%
  dplyr::mutate(rank = row_number()) %>%
  dplyr::mutate(labels = case_when(protein %in% labels ~ protein)) %>%
  dplyr::ungroup()

# Plot
plot <- df_long_log10 %>%
  ggplot(aes(x = rank, y = log10median, label = labels)) +
  geom_point(size = 2, alpha = 0.1, stroke = 0) +
  geom_text_repel(
    point.size = 3,
    size = 2,
    min.segment.length = 0.1,
    segment.size = 0.1,
    family = "Source Sans 3"
  ) +
  theme() +
  labs(y = "Log 10 median intensity", x = "Rank")

# Save plot
ggplot2::ggsave(snakemake@output[["fig_ranked_intensities"]], width = 75, height = 50, units = "mm", plot = plot)


# MYH2 vs. MYH7 volcano --------------------------------------------------
print_section_header("Creating MYH2 vs. MYH7 volcano (Fig. 2)")

labels <- c(
  "ATP2A2",
  "MYH3",
  "MYH7",
  "MYH6",
  "MYL3",
  "MYL1",
  "TNNC2",
  "MYL1",
  "MYH2",
  "ACTN3",
  "MYH1"
)

myh_volcano <- volcano_plot_myh(results_myh_ii_vs_i, labels, ylim = c(0, 45), xlim = c(-6, 6)) +
  theme(aspect.ratio = NULL)

# Save plot
ggplot2::ggsave(snakemake@output[["fig_myh_volcano"]], height = 50, width = 66, units = "mm", plot = myh_volcano)

# PCA & loadings ---------------------------------------------------------
print_section_header("Creating PCA plot (Fig. 2)")

# Subset SummarizedExperiment for only "pre" samples and filter for 70% valid values in at least one fiber-type
se_pca <- PhosR::selectGrps(se[, se$time == "pre"], se[, se$time == "pre"]$fiber_type, 1, n = 2)

# Prepare dataset
df_pca <- SummarizedExperiment::assay(se_pca) %>%
  t() %>%
  scale() %>% 
  as.data.frame() %>%
  merge(metadata, by = 0, all.x = TRUE) %>%
  dplyr::select(-c("time", "intervention")) %>%
  dplyr::relocate(id, .after = Row.names) %>%
  dplyr::relocate(fiber_type, .after = Row.names) %>%
  tibble::column_to_rownames("Row.names")

# Set seed for reproducible imputation
set.seed(99)

# Run PCA analysis
pca <- pcaMethods::pca(df_pca, method = "ppca", nPcs = 2)

# Merge analysis wit PCA data frame. Also creates a new column which concatenates fiber_type and id column for PCA coloring.
df_pca_plot <- merge(df_pca, scores(pca), by = 0)

#Retrieve convex hull (for plotting)
hulls <- df_pca_plot %>%
  dplyr::group_by(fiber_type) %>%
  dplyr::slice(chull(PC1, PC2))

# Visualize
plot_pca <- ggplot(df_pca_plot, aes(PC1, PC2, color = fiber_type)) +
  geom_polygon(data = hulls, aes(fill = fiber_type), alpha = 0.2, color = NA) + #Convex hull
  geom_point(size = 4) +
  scale_color_manual(values = c(I = myh7_color, II = myh2_color)) +
  scale_fill_manual(values = c(I = myh7_color, II = myh2_color)) +
  theme(
    legend.position = "none"
  ) +
  xlab(paste("PC1 (", round(pca@R2[1] * 100, digits = 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca@R2[2] * 100, digits = 1), "%)", sep = ""))

# Save plot
ggplot2::ggsave(snakemake@output[["fig_pca"]], width = 66, height = 50, units = "mm", plot = plot_pca)

print_section_header("Creating PC loadings plot (Fig. 2)")

# Retrieve loadings
loadings <- as.data.frame(loadings(pca)) %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(protein = rowname)

# Loading labels
loading_labels <- c(
  "TNNC1",
  "MYH7",
  "MYL3",
  "MYL2",
  "MYL1",
  "MYH6",
  "MYH3",
  "TNNI1",
  "TNNT1",
  "MYL6B",
  "LDHB",
  "MYH2",
  "TNNC2",
  "MYBPC2",
  "MYH8",
  "MYL11",
  "TNNI2",
  "ACTN3",
  "MYH1"
)

# Plot
loadings_plot <- loadings %>%
  dplyr::mutate(
    color = case_when(
      protein %in% loading_labels & PC1 > 0 ~ "myh7",
      protein %in% loading_labels & PC1 < 0 ~ "myh2",
      T ~ "na"
    ),
    label_col = ifelse(protein %in% loading_labels, protein, NA)
  ) %>%
  ggplot(aes(x = PC1, y = PC2, text = protein)) +
  geom_point(aes(color = color, alpha = color), size = 1, shape = 16) +
  geom_text_repel(aes(label = label_col),
    point.size = 1,
    size = 1,
    min.segment.length = 0.1,
    force = 0.3,
    segment.size = 0.1,
    na.rm = T, 
    family = "Source Sans 3"
  ) +
  scale_color_manual(values = c(myh7 = myh7_color, myh2 = myh2_color, na = "grey")) +
  scale_alpha_manual(breaks = c("myh7", "myh2", "na"),
                       values=c(1, 1, 0.25))+
  theme(
    legend.position = "none"
  )

# Save plot
ggplot2::ggsave(snakemake@output[["fig_pca_loadings"]], width = 66, height = 50, units = "mm", plot = loadings_plot)

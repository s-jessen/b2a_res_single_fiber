# ---- Setup ----
source("renv/activate.R")
source(snakemake@input[["functions"]])
source(snakemake@input[["settings"]])

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggvenn)
})

# Load data --------------------------------------------------------------
df_long <- readRDS(snakemake@input[["df_long"]])
df_long_l2fc_mean <- readRDS(snakemake@input[["df_long_l2fc_mean"]])
gsea_all <- readRDS(snakemake@input[["gsea_all"]])

# Load results  ----------------------------------------------------------
# Terbutaline
results_ter_i <- readRDS(snakemake@input[["results_ter_i"]])
results_ter_ii <- readRDS(snakemake@input[["results_ter_ii"]])

# Resistance training
results_res_i <- readRDS(snakemake@input[["results_res_i"]])
results_res_ii <- readRDS(snakemake@input[["results_res_ii"]])

# Venn diagram -----------------------------------------------------------
print_section_header("Creating Venn diagram plot (Fig. 4)")

plot_venn <- ggvenn::ggvenn(
  list(
    "ter_i" = dplyr::filter(results_ter_i, xiao < 0.05)$protein,
    "ter_ii" = dplyr::filter(results_ter_ii, xiao < 0.05)$protein,
    "res_i" = dplyr::filter(results_res_i, xiao < 0.05)$protein,
    "res_ii" = dplyr::filter(results_res_ii, xiao < 0.05)$protein
  ),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_elements = F
)

ggplot2::ggsave(snakemake@output[["fig_venn"]], width = 100, height = 100, units = "mm", plot = plot_venn)


# Top 25 regulated proteins ----------------------------------------------
# Identify the top 25 most upregulated proteins
print_section_header("Creating Top25 upregulated proteins plot (Fig. 4)")

sum_most_up <- df_long_l2fc_mean %>%
  dplyr::filter(n > 4) %>%
  dplyr::group_by(protein) %>%
  dplyr::summarize(sum = sum(mean_l2fc, na.rm = T)) %>%
  dplyr::slice_max(order_by = sum, n = 25) %>%
  dplyr::pull(protein)

# Dataframe for plot
df_top25 <- df_long_l2fc_mean %>%
  dplyr::filter(protein %in% sum_most_up) %>%
  dplyr::mutate(color = paste(intervention, fiber_type, sep = "_")) %>% 
  dplyr::mutate(protein = factor(protein, levels = rev(sum_most_up)))

# Save csv
readr::write_csv(df_top25, snakemake@output[["top25_upregulated"]])

# Plot
plot_top25 <- ggplot2::ggplot(df_top25, aes(x = protein, y = mean_l2fc, fill = color)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", size = 7)
  ) +
  scale_fill_manual(
    values = c(
      resistance_I = res_i_color,
      resistance_II = res_ii_color,
      terbutaline_I = ter_i_color,
      terbutaline_II = ter_ii_color
    ),
    labels = c(
      "resistance_I" = "Resistance Type I",
      "resistance_II" = "Resistance Type II",
      "terbutaline_I" = "B2A, Type I",
      "terbutaline_II" = "B2a, Type II"
    )
  ) +
  labs(
    y = "Summed log2fold change",
    x = NULL)


ggplot2::ggsave(snakemake@output[["fig_top25_upregulated"]], width = 120, height = 50, units = "mm", plot = plot_top25)


# MYH6 plot --------------------------------------------------------------
print_section_header("Creating MYH6 plot (Fig. 4)")

myh6_abundances <- df_long %>%
  dplyr::filter(protein == "MYH6")

readr::write_csv(myh6_abundances, snakemake@output[["myh6_abundances"]])

plot_myh6 <- myh6_abundances %>%
  dplyr::mutate(fiber_intervention = interaction(intervention, fiber_type, sep = "_")) %>%
  ggplot(aes(x = factor(time, levels = c("pre", "post")), y = abundance, color = fiber_intervention)) +
  geom_point(size = 1) +
  geom_path(aes(group = interaction(id, fiber_type))) +
  facet_grid(~intervention) +
  scale_color_manual(
    values = c(
      resistance_I = res_i_color,
      resistance_II = res_ii_color,
      terbutaline_I = ter_i_color,
      terbutaline_II = ter_ii_color
    ),
    labels = c(
      "resistance_I" = "Resistance Type I",
      "resistance_II" = "Resistance Type II",
      "terbutaline_I" = "B2A, Type I",
      "terbutaline_II" = "B2A, Type II"
    )
  ) +
  scale_x_discrete(labels = c(pre = "Pre", post = "Post")) +
  theme(
    legend.position = "none",
    strip.text = element_blank()
  ) +
  labs(x = "", y = "log2 (abundance)")

ggplot2::ggsave(snakemake@output[["fig_myh6"]], width = 35, height = 40, units = "mm", plot = plot_myh6)

# Correlations - type II vs. type I --------------------------------------
print_section_header("Creating correlation plot (type II vs. type I) (Fig. 2)")
# Create plotting data frame
df <- df_long_l2fc_mean %>%
  dplyr::select(protein, intervention, fiber_type, mean_l2fc) %>%
  tidyr::pivot_wider(
    names_from = fiber_type,
    values_from = mean_l2fc
  )

readr::write_csv(df, snakemake@output[["cor_i_vs_ii"]])

# Prepare labels
labels <- c(
  "MYH6",
  "MYL2",
  "MYL3",
  "TNNC1",
  "TNNI1",
  "ANXA1",
  "CTSG",
  "MYH11",
  "MYH8"
)

# Plot
plot_cor_i_vs_ii <- df %>%
  dplyr::mutate(label_col = ifelse(protein %in% labels, protein, NA)) %>%
  ggplot(aes(x = I, y = II)) +
  geom_point(alpha = 0.1, shape = 16, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  labs(
    x = "Log2fold change, Type 1",
    y = "Log2fold change, Type 2",
    title = ""
  ) +
  geom_smooth(
    method = "lm",
    linewidth = 0.25,
    color = "black",
    se = F
  ) +
  geom_text_repel(aes(label = label_col),
    point.size = 1,
    size = 1,
    min.segment.length = 0.1,
    force = 0.3,
    segment.size = 0.1,
    family = "Source Sans 3"
  ) +
  facet_grid(~intervention,
    labeller = as_labeller(c(resistance = "RES", terbutaline = "B2A"))
  )

# Save plot
ggplot2::ggsave(snakemake@output[["fig_cor_i_vs_ii"]], width = 100, height = 50, units = "mm", plot = plot_cor_i_vs_ii)


# Correlations - RES vs. B2A ---------------------------------------------
print_section_header("Creating correlation plot (RES vs. B2a) (Fig. 4)")

# Create plotting data frame
df <- df_long_l2fc_mean %>%
  dplyr::select(protein, intervention, fiber_type, mean_l2fc) %>%
  tidyr::pivot_wider(
    names_from = intervention,
    values_from = mean_l2fc
  )

readr::write_csv(df, snakemake@output[["cor_res_vs_b2a"]])

# Prepare labels
labels <- c(
  "MYH6",
  "MYL2",
  "MYL3",
  "TNNC1",
  "TNNI1",
  "ANXA1",
  "CTSG",
  "MYH11",
  "MYH8"
)

# Plot
plot_cor_res_vs_b2a <- df %>%
  dplyr::mutate(label_col = ifelse(protein %in% labels, protein, NA)) %>%
  ggplot(aes(x = resistance, y = terbutaline)) +
  geom_point(alpha = 0.1, shape = 16, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  labs(
    x = "Log2fold change, RES",
    y = "Log2fold change, B2A",
    title = ""
  ) +
  geom_smooth(
    method = "lm",
    linewidth = 0.25,
    color = "black",
    se = F
  ) +
  geom_text_repel(aes(label = label_col),
    point.size = 1,
    size = 1,
    min.segment.length = 0.1,
    force = 0.3,
    segment.size = 0.1,
    family = "Source Sans 3"
  ) +
  facet_grid(~fiber_type,
    labeller = as_labeller(c(I = "Type I", II = "Type II"))
  )

# Save plot
ggplot2::ggsave(snakemake@output[["fig_cor_res_vs_b2a"]], width = 100, height = 50, units = "mm", plot = plot_cor_res_vs_b2a)


# GSEA results -----------------------------------------------------------
print_section_header("Creating GSEA plot (Fig. 4)")

# Convert to data frame
df_gsea <- gsea_all@compareClusterResult

# Mark terms of interest
terms <- c(
  "cellular respiration",
  "mitochondrion organization",
  "cytoplasmic translation",
  "cytoskeleton organization",
  "muscle contraction",
  "muscle cell development",
  "striated muscle contration",
  "fatty acid beta-oxidation",
  "tricarboxylic acid cycle",
  "gluconeogenesis",
  "cell differentiation",
  "cell development",
  "translation",
  "cell population proliferation",
  "ribosome assembly",
  "actin filament bundle assembly",
  "striated muscle cell differentation",
  "ribosome biogenesis",
  "actin cytoskeleton organization",
  "muscle structure development",
  "G protein-coupled receptor signaling pathway",
  "striated muscle adaptation",
  "myofibril assembly",
  "actin filament organization"
)

# Clean up data frame
df_gsea_plot <- df_gsea %>%
  dplyr::filter(Description %in% terms) %>%
  dplyr::mutate(direction = ifelse(NES > 0, "Enriched", "Depleted")) %>%
  dplyr::mutate(gene_ratio = (str_count(core_enrichment, "/") + 1) / setSize)

# Determine plotting order
order <- df_gsea_plot %>%
  dplyr::group_by(Description) %>%
  dplyr::summarize(mean = mean(NES)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(abs(mean))

# Set order
df_gsea_plot <- df_gsea_plot %>%
  dplyr::mutate(
    Description = factor(Description, levels = (order$Description)),
    Cluster = factor(Cluster, levels = c(
      "ter_i",
      "ter_ii",
      "res_i",
      "res_ii"
    ))
  )

# Plot
plot_gsea <- ggplot(df_gsea_plot, aes(x = Cluster, y = Description, size = gene_ratio)) +
  geom_point(aes(fill = p.adjust), shape = 21, color = "black", stroke = 0.25) +
  scale_size(
    range = c(2, 10),
    name = "Gene Ratio"
  ) +
  scale_fill_gradient(
    low = "#e26664",
    high = "#327ebd",
    name = "p.adjust"
  ) +
  scale_x_discrete(labels = c(
    ter_i = "B2A\nType I",
    ter_ii = "B2A\nType II",
    res_i = "RES\nType I",
    res_ii = "RES\nType II"
  )) +
  facet_grid(~ factor(direction, levels = c("Enriched", "Depleted"))) +
  theme(
    legend.title = element_text(),
    axis.text.x = element_text(size = 7)
  ) +
  labs(y = NULL, x = NULL)

# Save plot
ggplot2::ggsave(snakemake@output[["fig_gsea"]], height = 95, width = 145, units = "mm", plot = plot_gsea)

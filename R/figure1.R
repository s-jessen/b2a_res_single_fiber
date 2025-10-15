# ---- Setup ----
source("renv/activate.R")
source(snakemake@input[["functions"]])
source(snakemake@input[["settings"]])

# Load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
})

# Load data --------------------------------------------------------------
functional_data <- readxl::read_xlsx(snakemake@input[["functional_data"]])

# Functional measures -----------------------------------------------------
print_section_header("Creating functional measures plots (Fig. 1)")

df <- functional_data %>% 
    dplyr::group_by(id, treatment) %>% 
    dplyr::summarize(
        mpo_change = mpo[time == "post"] - mpo[time == "pre"],
        ppo_change = ppo[time == "post"] - ppo[time == "pre"],
        lean_mass_change = lean_mass[time == "post"] - lean_mass[time == "pre"],
        one_rm_change = one_rm_leg_extension[time == "post"] - one_rm_leg_extension[time == "pre"],
        .groups = "drop"
        ) %>% 
    tidyr::pivot_longer(
        cols = c("ppo_change", "mpo_change", "lean_mass_change", "one_rm_change"),
        names_to = "variable",
        values_to = "change"
    ) %>% 
    dplyr::mutate(variable = factor(variable, levels = c(
        "lean_mass_change",
        "one_rm_change",
        "ppo_change",
        "mpo_change"
    )))

# Plot
plot_functional_measures <- ggplot(df, aes(x = treatment, y = change, fill = treatment))+
    geom_bar(stat = "summary", fun = "mean")+
    #geom_boxplot()+
    geom_point(size = 1, alpha = 0.6, stroke = 0)+
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.25)+
    scale_fill_manual(
        values = c(
            "terbutaline" = terbutaline_color,
            "resistance" = resistance_color
            ),
        labels = c(
            "terbutaline" = "B2A",
            "resistance" = "RES"
            )
        )+
    scale_x_discrete(labels = c(
        "terbutaline" = "B2A",
        "resistance" = "RES"
            )
        )+
    labs(
        x = NULL,
        y = "Change",
        ) +
    scale_y_continuous(
        expand = expansion(mult = c(0.1, 0.5))
        )+
    theme(
        legend.position = "none"
    ) +
    facet_wrap(
        ~variable,
        ncol = 5,
        scales = "free",
        labeller = as_labeller(c(
            lean_mass_change = "Lean mass (kg)",
            one_rm_change = "1RM\nLeg Extension (kg)",
            mpo_change = "Mean\nPower Output (W)",
            ppo_change = "Peak\nPower Output (W)"
        ))
    )

ggsave(snakemake@output[["fig_functional_measures"]], width = 85, height = 40, units = "mm", plot = plot_functional_measures)

# Export to supplementary files ------------------------------------------
writexl::write_xlsx(functional_data, snakemake@output[["supplementary_functional"]])
suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

#Colors
upregulated_color = "#a51e22"
downregulated_color = "#1e5886"

myh7_color <- "#d44a4f"
myh2_color <- "#7b98cf"

terbutaline_color = "#C06C84"
resistance_color = "#355C7D"

res_i_color <- "#e65878"
res_ii_color <- "#f2a629"
ter_i_color <- "#113f5b"
ter_ii_color <- "#785396"

#Set ggplot theme
theme_set(
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black",
                                fill=NA,
                                linewidth = 0.25),
    panel.grid.minor=element_blank(),
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "white", linewidth = 0),
    plot.title = element_text(size = 8,
                              margin = margin(0, 0, 0, 0),
                              hjust = 0),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.25),
    axis.text = element_text(color = "black", family = "Source Sans 3"),
    axis.title = element_text(size = 8, family = "Source Sans 3"),
    text = element_text(size = 8,
                        color = "black",
                        family = "Source Sans 3"),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.size = (unit(3, "mm")),
    legend.position = "right",
    legend.justification = c(0, 1),
    legend.box.just = "left",
    legend.margin = margin(0, 0, -5, 0),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    title = element_text(size = 6, family = "Source Sans 3"),
    strip.background = element_blank()
  )
)

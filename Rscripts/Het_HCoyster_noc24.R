#!/usr/bin/env Rscript
rm(list = ls(all = TRUE))

## Load required packages
library(tidyverse)
library(cowplot)

## Define the paths and variable directly for running without command-line arguments
indir <- "~/Desktop/GitHub/SteinOedulis/figures_1strun/heterozygosity/sfs"
outdir <- "~/Desktop/GitHub/SteinOedulis/figures_1strun/heterozygosity"
sample_table_path <- "~/Desktop/GitHub/Gskimming/01_infofiles/SteinOedulis_loco.tsv"

## read in the data
sample_table <- read_tsv(sample_table_path)
het <- NULL
for (i in 1:nrow(sample_table)){
  sfs_i <- read_delim(str_c(indir, "/", sample_table$sample_name[i], ".sfs"), col_names = FALSE) %>%
    as.matrix() %>%
    as.vector()
  het_i <- sfs_i[2]/sum(sfs_i, na.rm = TRUE)
  het <- c(het, het_i)
}
het
## plot a boxplot
het_plot <- mutate(sample_table, het=het) %>%
  ggplot(aes(x=population, y=het, fill=population)) +
  geom_boxplot() +
  labs(
    x = "Population",
    y = "Heterozygosity") +
  geom_jitter(height = 0, alpha = 1) +
  coord_flip() +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "right",
        legend.title = element_text(color = "#000000", size = 13),
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 13, face = "bold"),
        axis.ticks = element_line(color = "#000000", size = 0.3),
        axis.line = element_line(colour = "#000000", size = 0.3)) +
  guides(fill = guide_legend(title = "Population", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))
het_plot

summary(het)
## Save the plot
ggsave(filename = str_c(outdir, "/heterozygosity_plot.png"), plot = het_plot)

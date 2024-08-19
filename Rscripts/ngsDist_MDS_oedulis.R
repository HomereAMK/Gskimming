### The BEGINNING ~~~~~
##
# ~ Plots BSG_AtlanticHerring--MDS | First written by Filipe G. Vieira with later modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())

# Loads required packages ~
pacman::p_load(optparse, tidyverse, plyr, RColorBrewer, extrafont, ggforce, ggstar, plotly)

# Reads data ~
data <- read.table("~/Desktop/GitHub/Gskimming/00_data/ngsDist/Oedulis/Raw_cov/MDS/combined.subsetted.filtered_chr4568.mds")
n <- ncol(data)
annot <- read_tsv("GitHub/Gskimming/01_infofiles/SteinOedulis_loco.tsv")
colnames(annot) <- c("Sample_ID","species", "Population", "bam")
data <- data %>% rownames_to_column(var = "Sample_ID")
data

#D1_9.24236936859321 D2_3.29265946212272
#D1_8.37412651756432 D2_2.67957562546938

# Merge MDS with AnnotUp
data <- merge(data, annot, by = "Sample_ID")
MDS_12 <-
  ggplot(data, aes_string(x = "D1_8.37412651756432", y = "D2_2.67957562546938", fill = "Population")) +
  geom_point(alpha = .9, size = 2.75, shape = 21, colour = "#000000") +
  #scale_fill_manual(values = c("#377eb8", "#4daf4a", "#d01c8b", "#e66101")) +
  scale_x_continuous("Dimension 1 (8.37%)",
                     #breaks = c(-0.075, -0.05, -0.025, 0, 0.025),
                     #labels = c("-0.075", "-0.05", "-0.025", "0", "0.025"),
                     #limits = c(-0.073, 0.03),
                     expand = c(.015, .015)) +
  scale_y_continuous("Dimension 2 (2.68%)",
                     #breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     #labels = c("-0.05", "-0.025", "0", "0.025", "0.05"), 
                     #limits = c(-0.0525, 0.0525),
                     expand = c(.015, .015)) +
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

ggsave(MDS_12, file = "~/Desktop/GitHub/Gskimming/02_figures/Oedulis/ngsdist/MDS/Oedulis_noinv_ngsdist_mds_percountry_pop_PC1PC2.png",scale = 1, dpi = 600)



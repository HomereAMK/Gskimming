### The BEGINNING ~~~~~
##
# ~ Plots BSG_AtlanticHerring--MDS | First written by Filipe G. Vieira with later modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())

# Loads required packages ~
pacman::p_load(optparse, tidyverse, plyr, RColorBrewer, extrafont, ggforce, ggstar, plotly)

# Reads data ~
data <- read.table("~/Desktop/GitHub/Gskimming/00_data/ngsDist/Oedulis/Oedulis_beagle.mds")
n <- ncol(data)
annot <- read_csv("~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/oedulis_sra.csv")
colnames(annot) <- c("Sample_ID","species", "Population", "bam")
data <- data %>% rownames_to_column(var = "Sample_ID")
str(data)
#Create $population
data <- data %>%
  mutate(population = substr(Sample_ID, 1, 4))

#D1_3.95925929224189 D2_2.49150038136162

# Merge MDS with AnnotUp
data <- merge(data, annot, by = "Sample_ID")
MDS_12 <-
  ggplot(data, aes_string(x = "D1_3.95925929224189", y = "D2_2.49150038136162", fill = "population")) +
  geom_point(alpha = .9, size = 2.75, shape = 21, colour = "#000000") +
  #scale_fill_manual(values = c("#377eb8", "#4daf4a", "#d01c8b", "#e66101")) +
  scale_x_continuous("Dimension 1 (3.96%)",
                     #breaks = c(-0.075, -0.05, -0.025, 0, 0.025),
                     #labels = c("-0.075", "-0.05", "-0.025", "0", "0.025"),
                     #limits = c(-0.073, 0.03),
                     expand = c(.015, .015)) +
  scale_y_continuous("Dimension 2 (2.49%)",
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

ggsave(MDS_12, file = "~/Desktop/GitHub/Gskimming/02_figures/Oedulis/ngsdist/Oedulis_ngsdist_mds__pop_PC1PC2.png",scale = 1, dpi = 600)



### The BEGINNING ~~~~~
##
# ~ Plots PopGenEstimates | By Marie-Christine Rufener & George Pacheco modified by Homère J. Alves Monteiro


# Cleans the environment ~ 
rm(list=ls())
rm(list=ls(all.names=TRUE))


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(scales, extrafont, dplyr, grid, lubridate, cowplot, egg, tidyverse, lemon, stringr, reshape, gghighlight)


# Load helper function ~
#source("utilities.R")



# Loads datasets ~
#PopGen <- read.table("BSG_Turbot_Ind66.PopGenEstimates.txt", sep = "\t", header = FALSE)
#PopGen <- read.table("~/Desktop/Scripts/Data/PopGenEstimates/Dataset_I/ThetaSummaries_Feb22--Unfolded_PopGen_ALLpop.PopGenEstimates.txt", sep = "\t", header = FALSE)
PopGen <- read.table("~/Desktop/Scripts/Data/PopGenEstimates/EUostrea/Dec22--Ind581.PopGenEstimates.txt", sep = "\t", header = FALSE)

#Hets <- read.table("FPG--Heterozygosity.txt", sep = "\t", header = FALSE); head(Hets)


# Adds column names to datasets ~
colnames(PopGen) <- c("Population", "NSites", "Nucleotide_Diversity", "Watterson_Theta", "Tajima_D")
#colnames(Hets) <- c("Sample_ID", "Population", "Het", "DataType"); head(Hets)

# Stats ~
summary(PopGen)
means <- aggregate(PopGen[,2:5], by=list(Population=PopGen$Population), FUN=mean)
summary(means)
PopGenMin_Max <- means %>%
  group_by(Population) %>%
  summarise(
    NSites_min = min(NSites),
    NSites_max = max(NSites),
    Nucleotide_Diversity_min = min(Nucleotide_Diversity),
    Nucleotide_Diversity_max = max(Nucleotide_Diversity),
    Watterson_Theta_min = min(Watterson_Theta),
    Watterson_Theta_max = max(Watterson_Theta),
    Tajima_D_min = min(Tajima_D),
    Tajima_D_max = max(Tajima_D)
  )
print(PopGenMin_Max)
# Corrects Population names in Hets:
#levels(Hets$Population <- sub("FeralUT", "SaltLakeCity", Hets$Population))
#levels(Hets$Population <- sub("FeralVA", "Virginia", Hets$Population))


# Tidies up data frames ~
levels(PopGen$Population)
PopGen$Population <- as.factor(gsub(" ", "", PopGen$Population))


# Sets BioStatus from character -> factor (better for data manipulation & necessary for plotting)
#PopGen$BioStatus <- as.factor(PopGen$BioStatus)
#Hets$BioStatus <- as.factor(Hets$BioStatus)


# Converts DF from wide into long ~
PopGenUp <- gather(PopGen, Estimate, Value, "Nucleotide_Diversity", "Watterson_Theta", "Tajima_D")


# Adds data ID column to each DF ~
PopGenUp$ID <- factor(paste("PopGen"))
#Hets$ID <- factor(paste("Hets"))


# Binds the 2 DFs based on common columns ~
#fulldf <- mybind(PopGenUp, Hets)


# Includes label for empty factor level (related to PHS) ~
#idx <- which(fulldf$ID == "Hets")
#fulldf[idx,"Estimate"] <- rep("PHS", length(idx))
#fulldf$Estimate <- factor(fulldf$Estimate)


# Reorders factor levels ~
PopGenUp$Estimate <-
  factor(PopGenUp$Estimate, ordered = T, levels = c("Nucleotide_Diversity",
                                                  "Watterson_Theta",
                                                  "Tajima_D"))

# Corrects the facet labels ~
ylabel <- c("Nucleotide_Diversity" = "Nucelotide Diversity",
            "Tajima_D"= "Tajima's D",
            "Watterson_Theta" = "Watterson's Theta")


# Corrects population names ~
#levels(PopGenUp$Population <- sub("NorthSea", "North Sea", PopGenUp$Population))
#levels(PopGenUp$Population <- sub("BalticSea", "Baltic Sea", PopGenUp$Population))



# Reorders populations ~
PopGenUp$Population <- factor(PopGenUp$Population, ordered = T,
                              levels = c("MOLU", "ZECE", "CRES",
                                         "ORIS","CORS", "PONT",  "RIAE",
                                         "MORL",
                                         "USAM",
                                         "TOLL", "COLN", "BARR",
                                         "TRAL", "CLEW",
                                         "RYAN",
                                         "GREV", "WADD",
                                         "NISS","LOGS","VENO", "HALS", "THIS",
                                         "KALV", "HYPP",
                                         "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                                         "INNE","VAGS", "AGAB", "OSTR"))


# Custom y-axis breaks ~
breaks_fun <- function(x){
  if (max(x) > 0.1){
    c(0.20, 0.40, 0.60)}
  else{
    c(0.0012, 0.0013, 0.0014, 0.0015)}}


# Custom y-axis labels ~
plot_index_labels <- 0
label_fun <- function(x) {
  plot_index_labels <<- plot_index_labels + 1L
  switch(plot_index_labels,
         scales::label_number(accuracy = 0.1, suffix = "X")(x),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(x),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(x),
         scales::label_number(scale = 1/1000000, accuracy = 1, big.mark = "", suffix = "M")(x))}


# Creates plot ~
PopGenEstimates_Plot <- 
 ggplot() +
  geom_point(data = PopGenUp,
             aes(x = Population, y = Value, fill = Population), colour = "#000000", shape = 21, size = 3.5, alpha = .9) +
  facet_rep_grid(Estimate ~. , scales = "free", labeller = labeller(Estimate = ylabel)) +
  scale_fill_manual(values =c(  "#A02353", "#A02353", "#A02353",
                                "#AD5B35",
                                "#ad7358",
                                "#CC480C",  "#CC480C",
                                "#969696",
                                "#000000",
                                "#D38C89", "#D38C89", "#D38C89",
                                "#C89AD1", "#C89AD1",
                                "#7210A0",
                                "#91BD96", "#91BD96",
                                "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                                "#45D1F7", "#45D1F7",
                                "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                                "#240377", "#240377", "#240377", "#240377" ))+
  #scale_colour_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 16),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(fill = guide_legend(title = "Populations:", title.theme = element_text(size = 12, face = "bold"),
                             subtitle ="Genome-wide",
                             label.theme = element_text(size = 14),
                             override.aes = list(size = 5, alpha = .9)), colour = "none")


# Saves plot ~
ggsave(PopGenEstimates_Plot, file = "~/Desktop/Scripts/EUostrea/Figures/PopGenEstimates/pi_theta_Taj_PopGenEstimates_14mar23.pdf",
  device = cairo_pdf, width = 12, height = 8, scale = 1.35, dpi = 600)















##### HET #####
Het <- read.table("~/Desktop/Scripts/Data/PopGenEstimates/EUostrea/GEO_EUostrea--AllSamples_setMinDepth600_setMaxDepth1200--HET.23jan23.txt", sep = "\t", header = FALSE)
Het <- Het[,-2]
colnames(Het) <- c("Sample_ID", "Heterozygosity")
# Drop the columns of the dataframe

#remove path bamfile
Het$Sample_ID=sub(".bam","",Het$Sample_ID)
Het$Sample_ID=sub(".+/","",Het$Sample_ID)
length(Het$Sample_ID)
Het$Sample_ID
Het$Sample_ID=sub(".depth.gz","",Het$Sample_ID)
length(Het$Sample_ID)
pop=substr(Het$Sample_ID,0,4) #name only 4characters -> pop
unique(pop)
#combine
Het_pop = as.data.frame(cbind(pop, Het)) 
unique(pop)
#order pop
Het_pop$pop <- factor(Het_pop$pop, ordered = T,
                      levels = c("MOLU", "ZECE", "CRES",
                                 "ORIS","CORS", "PONT",  "RIAE",
                                 "MORL",
                                 "USAM",
                                 "TOLL", "COLN", "BARR",
                                 "TRAL", "CLEW",
                                 "RYAN",
                                 "GREV", "WADD",
                                 "NISS","LOGS","VENO", "HALS", "THIS",
                                 "KALV", "HYPP",
                                 "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                                 "INNE","VAGS", "AGAB", "OSTR"))

#Plot
Het_plot <- 
  ggplot() +
  geom_point(data = Het_pop,
             aes(x = pop, y = Heterozygosity, fill = pop), colour = "#000000", shape = 21, size = 3.5, alpha = 0.6) +
  scale_fill_manual(values =c( "#A02353", "#A02353", "#A02353",
                                 "#AD5B35",
                                 "#ad7358",
                                 "#CC480C",  "#CC480C",
                                 "#969696",
                                 "#000000",
                                 "#D38C89", "#D38C89", "#D38C89",
                                 "#C89AD1", "#C89AD1",
                                 "#7210A0",
                                 "#91BD96", "#91BD96",
                                 "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                                 "#45D1F7", "#45D1F7",
                                 "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                                 "#240377", "#240377", "#240377", "#240377" ))+
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 16),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(fill = guide_legend(title = "Populations:", title.theme = element_text(size = 12, face = "bold"),
                             label = TRUE,
                             label.theme = element_text(size = 14),
                             override.aes = list(size = 5, alpha = .9)), colour = "none")

ggsave(Het_plot, file = "~/Desktop/Scripts/EUostrea/Figures/PopGenEstimates/Genome-wide_Heterozygosity_dotplot_Feb23.pdf",
       device = cairo_pdf, width = 12, height = 8, scale = 1.35, dpi = 600)



#geom_violin()
#Plot
Het_plot_violin <- 
  ggplot() +
  geom_violin(data = Het_pop,
             aes(x = pop, y = Heterozygosity, fill = pop), colour = "#000000", size = 0.3, alpha = .8) +
  scale_fill_manual(values =c( "#A02353", "#A02353", "#A02353",
                               "#AD5B35",
                               "#ad7358",
                               "#CC480C",  "#CC480C",
                               "#969696",
                               "#000000",
                               "#D38C89", "#D38C89", "#D38C89",
                               "#C89AD1", "#C89AD1",
                               "#7210A0",
                               "#91BD96", "#91BD96",
                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                               "#45D1F7", "#45D1F7",
                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                               "#240377", "#240377", "#240377", "#240377" ))+
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 16),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(fill = guide_legend(title = "Populations:", title.theme = element_text(size = 12, face = "bold"),
                             label = TRUE,
                             label.theme = element_text(size = 14),
                             override.aes = list(size = 0.3, alpha = .9)), colour = "none")

ggsave(Het_plot_violin, file = "~/Desktop/Scripts/EUostrea/Figures/PopGenEstimates/Genome-wide_Heterozygosity_violin_Feb23.pdf", 
       device = cairo_pdf, width = 12, height = 8, scale = 1.35, dpi = 600)












##### HET Inversion regions #####
Het <- read.table("~/Desktop/Scripts/Data/PopGenEstimates/EUostrea/InvReg/Reg04_24mar23--AllSamples_setMinDepth600_setMaxDepth1200--HET.txt", sep = "\t", header = FALSE)
Het <- Het[,-2]
colnames(Het) <- c("Sample_ID", "Heterozygosity")
# Drop the columns of the dataframe

#remove path bamfile
Het$Sample_ID=sub(".bam","",Het$Sample_ID)
Het$Sample_ID=sub(".+/","",Het$Sample_ID)
length(Het$Sample_ID)
Het$Sample_ID
Het$Sample_ID=sub(".depth.gz","",Het$Sample_ID)
length(Het$Sample_ID)
pop=substr(Het$Sample_ID,0,9) #name only 4characters -> pop
unique(pop)
pop <- substr(pop, nchar(pop)-3, nchar(pop))

#combine
Het_pop = as.data.frame(cbind(pop, Het)) 
unique(pop)
#order pop
Het_pop$pop <- factor(Het_pop$pop, ordered = T,
                      levels = c("MOLU", "ZECE", "CRES",
                                 "ORIS","CORS", "PONT",  "RIAE",
                                 "MORL",
                                 "USAM",
                                 "TOLL", "COLN", "BARR",
                                 "TRAL", "CLEW",
                                 "RYAN",
                                 "GREV", "WADD",
                                 "NISS","LOGS","VENO", "HALS", "THIS",
                                 "KALV", "HYPP",
                                 "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                                 "INNE","VAGS", "AGAB", "OSTR"))

head(Het_pop)
library(ggridges)
library(ggjoy)
library(ggplot2)
library(ggthemes)


#Plot
box <- ggplot(data = Het_pop, aes(x = pop, y = Heterozygosity, fill = pop)) +
  #geom_violin(alpha = 0.6, size = 0.3, color = "#000000") +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("#A02353", "#A02353", "#A02353",
                               "#AD5B35",
                               "#ad7358",
                               "#CC480C",  "#CC480C",
                               "#969696",
                               "#000000",
                               "#D38C89", "#D38C89", "#D38C89",
                               "#C89AD1", "#C89AD1",
                               "#7210A0",
                               "#91BD96", "#91BD96",
                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                               "#45D1F7", "#45D1F7",
                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                               "#240377", "#240377", "#240377", "#240377"))+
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.text.x = element_text(colour = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 16),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "none",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(colour = "none")+
  coord_flip()

last_plot()
ggsave(box, file = "~/Desktop/Scripts/EUostrea/Figures/PopGenEstimates/InvReg/Reg04_wide_boxplotHeterozygosity_mar23_flip.pdf", 
       device = cairo_pdf, width = 12, height = 8, scale = 1.35, dpi = 600)






#
H<-ggplot(data = Het_pop, aes(x = Heterozygosity , y = pop, fill = pop)) +
  geom_density_ridges2( aes( point_fill = pop), 
                                alpha = .2, point_alpha = 1, jittered_points = TRUE, scale = 3, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#A02353", "#A02353", "#A02353",
                               "#AD5B35",
                               "#ad7358",
                               "#CC480C",  "#CC480C",
                               "#969696",
                               "#000000",
                               "#D38C89", "#D38C89", "#D38C89",
                               "#C89AD1", "#C89AD1",
                               "#7210A0",
                               "#91BD96", "#91BD96",
                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                               "#45D1F7", "#45D1F7",
                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                               "#240377", "#240377", "#240377", "#240377"))+
  labs(x = "Population", y = "Heterozygosity",
       fill = "Population", color = "Population") +
  theme_joy() +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 16),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(fill = guide_legend(title = "Populations for Reg08:", title.theme = element_text(size = 12, face = "bold"),
                             label = TRUE,
                             label.theme = element_text(size = 14),
                             override.aes = list(size = 5, alpha = .9)), colour = "none")+

  coord_flip()
last_plot()
ggsave(H, file = "~/Desktop/Scripts/EUostrea/Figures/PopGenEstimates/InvReg/Reg05_wide_Heterozygosity_ggridges_24mar23.pdf",
       device = cairo_pdf, width = 12, height = 8, scale = 1.35, dpi = 600)



#Plot
box <- ggplot(data = Het_pop, aes(x = pop, y = Heterozygosity, fill = pop)) +
  #geom_violin(alpha = 0.6, size = 0.3, color = "#000000") +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("#A02353", "#A02353", "#A02353",
                               "#AD5B35",
                               "#ad7358",
                               "#CC480C",  "#CC480C",
                               "#969696",
                               "#000000",
                               "#D38C89", "#D38C89", "#D38C89",
                               "#C89AD1", "#C89AD1",
                               "#7210A0",
                               "#91BD96", "#91BD96",
                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                               "#45D1F7", "#45D1F7",
                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                               "#240377", "#240377", "#240377", "#240377"))+
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.text.x = element_text(colour = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 16),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "none",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(colour = "none")+
  coord_flip()

last_plot()
ggsave(box, file = "~/Desktop/Scripts/EUostrea/Figures/PopGenEstimates/InvReg/Reg05_wide_boxplotHeterozygosity_mar23_flip.pdf", 
       device = cairo_pdf, width = 12, height = 8, scale = 1.35, dpi = 600)



#geom_violin()
#Plot
Het_plot_violin <- 
  ggplot() +
  geom_violin(data = Het_pop,
              aes(x = pop, y = Heterozygosity, fill = pop), colour = "#000000", size = 0.3, alpha = .8) +
  scale_fill_manual(values =c( "#A02353", "#A02353", "#A02353",
                               "#AD5B35",
                               "#ad7358",
                               "#CC480C",  "#CC480C",
                               "#969696",
                               "#000000",
                               "#D38C89", "#D38C89", "#D38C89",
                               "#C89AD1", "#C89AD1",
                               "#7210A0",
                               "#91BD96", "#91BD96",
                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                               "#45D1F7", "#45D1F7",
                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                               "#240377", "#240377", "#240377", "#240377" ))+
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 16),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(fill = guide_legend(title = "Populations:", title.theme = element_text(size = 12, face = "bold"),
                             label = TRUE,
                             label.theme = element_text(size = 14),
                             override.aes = list(size = 0.3, alpha = .9)), colour = "none")

ggsave(Het_plot_violin, file = "~/Desktop/Scripts/EUostrea/Figures/PopGenEstimates/InvReg/Reg08_wide_boxplotHeterozygosity_mar23.pdf", 
       device = cairo_pdf, width = 12, height = 8, scale = 1.35, dpi = 600)













#
##
### The END ~~~~~
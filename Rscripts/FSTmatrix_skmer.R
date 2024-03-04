### The BEGINNING ~~~~
##
# ~ Plots FST heatmap from fst calculator module | By Eduardo Charvel  and  Homère J. Alves Monteiro
# Cleans the environment ~ 
rm(list=ls())


#install.package("pacman")
# Loads required packages ~
pacman::p_load(pheatmap, tidyverse, reshape2)

# Loads Fst table ~
# Don't forget to remove the "sample" string in the .txt file
data <- read.table("/Users/sjr729/Desktop/GitHub/Gskimming/01_infofiles/Fst/Skmer/jc-24.02-dist-mat_4x_Clupea_wcfst.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE,)
str(data)

#Remove the redundant first column, but the scalability of this need to be checked
long_df <- data %>%
  pivot_longer(
    cols = -sample,  # Select all columns except 'sample' to pivot
    names_to = "Pop2",  # Column that will contain the original column names
    values_to = "Weighted"  # Column that will contain the values from the pivoted columns
  ) %>%
  rename(Pop1 = sample)  # Rename 'sample' column to 'Pop1'
str(long_df)

# Loads Fst table from Angsd
Fst_angsd <- read.table("01_infofiles/Fst/Angsd/Mar24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist--Fst.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Fst_angsd) <- c("Pop1", "Pop2", "NumberOfSites", "Unweighted", "Weighted")
Fst_angsd <- Fst_angsd %>% dplyr::select(Pop1, Pop2, Weighted)
str(Fst_angsd)

# Ensure 'melted_Fst' and 'Fst_angsd' are loaded as data.frames
melted_Fst <- as.data.frame(long_df)
Fst_angsd <- as.data.frame(Fst_angsd)

# Create a new dataframe by joining 'Fst_angsd' with 'melted_Fst'
Fst <- left_join(Fst_angsd, melted_Fst, by = c("Pop1", "Pop2"))
Fst <- Fst %>% dplyr::select(Pop1, Pop2, Weighted.y)
colnames(Fst) <- c("Pop1", "Pop2","Weighted")

# Melts datasets ~
Fst_Pops = union(Fst$Pop1, Fst$Pop2)
n = length(Fst_Pops)


# Creates Fst-Fst matrix ~
Fst_Fst <- matrix(0, nrow = n, ncol = n, dimnames = list(Fst_Pops, Fst_Pops))
for (i in 1:nrow(Fst)) {
  Fst_Fst[Fst[i, "Pop1"], Fst[i, "Pop2"]] = Fst[i, "Weighted"]
  Fst_Fst[Fst[i, "Pop2"], Fst[i, "Pop1"]] = Fst[i, "Weighted"]}


# Gets Fst midpoint ~
Pops <- which(upper.tri(Fst_Fst), arr.ind = TRUE)
Fst_df <- data.frame(Site1 = dimnames(Fst_Fst)[[2]][Pops[, 2]],
                     Site2 = dimnames(Fst_Fst)[[1]][Pops[, 1]],
                     Value = Fst_Fst[Pops] %>% round(digits = 6))
Fstmiddle = max(Fst_df$Value) / 2


# Gets Fst label ~
Fst.label = expression(italic("F")[ST])


# Creates plot ~
Fst_Plot <-
  ggplot(data = Fst, aes(x = Pop1, y = Pop2, fill = Weighted)) +
  geom_tile(color = "#ffffff", lwd = 1.5, linetype = 1) +
  coord_fixed() +
  geom_text(aes(label = round(Weighted,digits = 3)), color = "#000000", size = 7) +
  scale_fill_gradient2(low = "#80b1d3", mid = "#fccde5", high = "#fb8072",
                       midpoint = Fstmiddle, name = Fst.label,
                       breaks = c(0.02, 0.04, 0.06, 0.08, 0.10),
                       labels = c("0.02", "0.04", "0.06", "0.08", "0.10"),
                       limits = c(0, .12)) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0.005, b = 0.005, r = .2, l = .2, unit = "cm"),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 18, angle = 60, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 18 ),
        axis.ticks = element_line(color = "#000000", size = .3),
        legend.position = "right",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.title = element_text(colour = "#000000", size = 18, face = "bold"),
        legend.text = element_text(colour = "#000000", size = 12))
Fst_Plot
 

#print(Fst_Plot)

ggsave(Fst_Plot, file = "~/Desktop/GitHub/Gskimming/02_figures/ClupeaAtmore/4x/Skmer/4x_Fst_n45_locality_JCcorr_Skmer_rawcoverage.png",
       scale = 1, width = 12, height = 12, dpi = 300)



### The BEGINNING ~~~~
##
# ~ Plots FST matrix from Angsd pipeline | Raw coverage | Homère J. Alves Monteiro

# Cleans the environment ~ 
rm(list=ls())

# Sets working directory ~
setwd(dir = "~/Desktop/GitHub/Gskimming/")

# Loads required packages ~
pacman::p_load(tidyverse, reshape2, Hmisc)

# Loads Fst table ~
Fst <- read.table("01_infofiles/ClupeaAtmore/Fst/Angsd/Mar24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist--Fst.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)


# Adds column names ~
colnames(Fst) <- c("Pop1", "Pop2", "NumberOfSites", "Unweighted", "Weighted")


# Melts datasets ~
Fst_Pops = union(Fst$Pop1, Fst$Pop2)
n = length(Fst_Pops)


# Creates Fst-Sites matrix ~
FstSites <- matrix(0, nrow = n, ncol = n, dimnames = list(Fst_Pops, Fst_Pops))
for (i in 1:nrow(Fst)) {
  FstSites[Fst[i, "Pop1"], Fst[i, "Pop2"]] = Fst[i, "NumberOfSites"]
  FstSites[Fst[i, "Pop2"], Fst[i, "Pop1"]] = Fst[i, "Weighted"]}


# Writes Fst-Sites matrix ~
write.table(FstSites, "01_infofiles/ClupeaAtmore/Fst/Angsd/Mar24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist--Fst-Sites.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Creates Fst-Fst matrix ~
Fst_Fst <- matrix(0, nrow = n, ncol = n, dimnames = list(Fst_Pops, Fst_Pops))
for (i in 1:nrow(Fst)) {
  Fst_Fst[Fst[i, "Pop1"], Fst[i, "Pop2"]] = Fst[i, "Weighted"]
  Fst_Fst[Fst[i, "Pop2"], Fst[i, "Pop1"]] = Fst[i, "Weighted"]}


# Writes Fst-Fst matrix ~
write.table(Fst_Fst, "01_infofiles/ClupeaAtmore/Fst/Angsd/Mar24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist--Fst-Fst.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


# Gets Fst midpoint ~
Pops <- which(upper.tri(Fst_Fst), arr.ind = TRUE)
Fst_df <- data.frame(Site1 = dimnames(Fst_Fst)[[2]][Pops[, 2]],
                     Site2 = dimnames(Fst_Fst)[[1]][Pops[, 1]],
                     Value = Fst_Fst[Pops] %>% round(digits = 6))
Fstmiddle = max(Fst_df$Value) / 2


# Gets Fst label ~
Fst.label = expression(italic("F")[ST])

str(Fst)

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

# Saves plot ~
ggsave(Fst_Plot, file = "02_figures/ClupeaAtmore/RawCov/Angsd/FST/Mar24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist--Fst.png",
        scale = 1, width = 12, height = 12, dpi = 600)

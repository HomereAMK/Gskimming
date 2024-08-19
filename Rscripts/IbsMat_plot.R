# Cleans the environment ~ 
rm(list=ls())

# Sets working directory ~
setwd("~/Desktop/")
#install.packages("RcppCNPy")
# Loads required packages ~
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)
source("GitHub/Gskimming/Rscripts/individual_pca_functions.R")

annot <- read_tsv("GitHub/Gskimming/01_infofiles/SteinOedulis_loco.tsv")


#Plot genome-wide PCoA with the ibsMat matrix
ibs_mat <- read_tsv("GitHub/Gskimming/00_data/Angsd/Oedulis/RawCov/MDS/NC_079172.1.ibsMat", col_names = F) %>% 
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
PCoA(ibs_mat, annot$sample_name, annot$population, 4, 1, 2, show.ellipse = F)

#For context this is Skmer distance matrix 
dist_skmer <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Skmer/Oedulis/RawCov/echarvel_skimming_pipeline/distance_matrix.txt", col_names = FALSE) %>%
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
summary(ibs_mat)



pheatmap(dist_skmer)
dimnames(ma)=list(indivi,indivi)



summary(dist_skmer)
dist_matrix <- as.matrix(dist_skmer[-1, -1])
dist_matrix <- apply(dist_matrix, 2, as.numeric)

# Extract all the distances excluding the diagonal
distances <- dist_matrix[lower.tri(dist_matrix)]

# Calculate the mean Jaccard distance
mean_jaccard_distance <- mean(distances)
mean_jaccard_distance

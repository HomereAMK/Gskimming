### The BEGINNING ~~~~~
##
# ~ Plot a PCoa/MDS with jc distance matrix from skmer -distance | from Homère J. Alves Monteiro and Eduardo Charvel
# ~ also Plot a PCoa/MDS with .ibs matrix from angsd -doIBS | from Homère J. Alves Monteiro and Eduardo Charvel

# R 4.3.2 GUI 1.80 Big Sur ARM build (8281)

### Clearing environment and setting working directory
rm(list = ls(all = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Loading necessary packages
library(tidyverse)
library(dplyr)
library(ape)
library(cowplot)
library(knitr)
library(scales)
library(miscTools)
if (!requireNamespace("dbscan", quietly = TRUE)) install.packages("dbscan")
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
library(dbscan)
library(mclust)

### Sourcing required functions
source("individual_mdsGskim_functions_hjam.R")


### Reading and processing ibs matrix
ibs_mat <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Angsd/Clupea/RawCov/MDS/screenfeb24_Atmore_modern_SNP_Inv_free_regions_minInd0.25.ibsMat", col_names = FALSE) %>%
  dplyr::select(1:nrow(.)) %>%
  as.matrix()

# Load complete modern Clupea atmore annotation file 
annot <- read_csv("../../Gskimming/01_infofiles/ClupeaAtmore/ClupeaModern_annot.csv")

### Performing PCoA with initial annotation
mds_plot <- PCoA(ibs_mat, annot$cleaned_id, annot$population,40, 1, 2, show.ellipse = FALSE, show.label = TRUE)
pcoa_table_genome_wide <- pcoa_table
write.csv(pcoa_table, "~/Desktop/GitHub/Gskimming/00_data/Angsd/Clupea/RawCov/MDS/Angsd_rawcov_AtmoreClupea_pcoa_table-mat_Mar24.csv", row.names = FALSE)


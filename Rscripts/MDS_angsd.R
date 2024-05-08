### The BEGINNING ~~~~~
##
# ~ Plot a PCoa/MDS with jc distance matrix from skmer -distance | from Homère J. Alves Monteiro and Eduardo Charvel
# ~ also Plot a PCoa/MDS with .ibs matrix from angsd -doIBS | from Homère J. Alves Monteiro and Eduardo Charvel

# R 4.3.2 GUI 1.80 Big Sur ARM build (8281)

### Clearing environment and setting working directory
rm(list = ls(all = TRUE))

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
source("~/Desktop/GitHub/Gskimming/Rscripts/individual_mdsGskim_functions_hjam.R")

#### Herring ####
### Reading and processing ibs matrix
ibs_mat <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Angsd/Clupea/RawCov/MDS/screenfeb24_Atmore_modern_SNP_Inv_free_regions_minInd0.25.ibsMat", col_names = FALSE) %>%
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
str(ibs_mat)
# Load complete modern Clupea atmore annotation file 
annot <- read_csv("~/Desktop/GitHub/Gskimming/01_infofiles/ClupeaAtmore/ClupeaModern_annot.csv")

### Performing PCoA with initial annotation
mds_plot <- PCoA(ibs_mat, annot$cleaned_id, annot$population,40, 1, 2, show.ellipse = FALSE, show.label = TRUE)
pcoa_table_genome_wide <- pcoa_table
write.csv(pcoa_table, "~/Desktop/GitHub/Gskimming/00_data/Angsd/Clupea/RawCov/MDS/Angsd_rawcov_AtmoreClupea_pcoa_table-mat_Mar24.csv", row.names = FALSE)




#### Magpie ####
ibs_mat <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Angsd/Magpie/RawCov/MDS/mar24_Magpie_SNP_minInd0.25_MinDepth450_MaxDepth1600.ibsMat", col_names = FALSE) %>%
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
# Load complete Magpie annot
annot <- read_csv("~/Desktop/GitHub/Gskimming/01_infofiles/Magpie/Magpie_annot.csv")
### Performing PCoA with initial annotation
mds_plot <- PCoA(ibs_mat, annot$sample_id, annot$population,40, 3, 4, show.ellipse = FALSE, show.label = TRUE)
ggsave(mds_plot, file = "~/Desktop/GitHub/Gskimming/02_figures/Magpie/RawCov/Angsd/MDS/Mds_1vs2_mar24_Magpie_SNP_minInd0.25_MinDepth450_MaxDepth1600.png",scale = 1, dpi = 600)

pcoa_table_genome_wide <- pcoa_table
write.csv(pcoa_table, "~/Desktop/GitHub/Gskimming/00_data/Angsd/Magpie/RawCov/MDS/Angsd_rawcov_Magpie_pcoa_table-mat_Apr24.csv", row.names = FALSE)

pcoa_table
pcoa_table_genome_wide_joined <- pcoa_table_genome_wide %>%
  dplyr::select(1:6) %>%
  left_join(annot_df_final, by=c("individual"="sample_id", "population"="population"))
mds_plot_sum <-pcoa_table_genome_wide_joined %>%
  ggplot(aes(x=dist_1, y=dist_2)) +
  geom_point(data=dplyr::select(pcoa_table_genome_wide_joined, -population), size = 0.1, color="grey") +
  geom_point(size=1, mapping = aes(color=population)) +
  facet_wrap(~population) +
  theme_cowplot() +
  theme(axis.text = element_blank())

## The BEGINNING ~~~~~
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
library(purrr)
library(scales)
library(miscTools)
if (!requireNamespace("dbscan", quietly = TRUE)) install.packages("dbscan")
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
library(dbscan)
library(mclust)
library(ggrepel)
### Sourcing required functions
source("~/Desktop/GitHub/Gskimming/Rscripts/individual_mdsGskim_functions_hjam.R")

### Reading and processing distance matrix
dist_jc <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Skmer2/Oedulis/RawCov/skmer2_skmer1lib_oedulis_bbmapreads_mappedtoGenome_p2_Herring_15.09.txt", col_names = FALSE) %>%
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
str(dist_jc)

## formatting
pairwise_distances <- dist_jc[-1, -1]
ibs_mat <- as.data.frame(as.matrix(pairwise_distances))
ibs_mat[is.na(ibs_mat) | ibs_mat == "nan"] <- 0
ibs_mat[is.na(ibs_mat) | ibs_mat == "-nan"] <- 0
ibs_mat[is.na(ibs_mat) | ibs_mat == "None"] <- 0
ibs_mat[is.na(ibs_mat) | ibs_mat == "none"] <- 0

#str(ibs_mat)

## Convert in numerical | necessary only for distance matrix from skmer2
# Assuming ibs_mat is already loaded and is as shown in the str() output
#ibs_mat <- ibs_mat %>%
#  mutate(across(everything(), ~as.numeric(.))) %>%  # Convert all columns to numeric
#  mutate(across(everything(), ~replace_na(., 0)))   # Replace NA with 0

# Check the structure and summary to confirm changes
str(ibs_mat)
summary(ibs_mat)

# Check for any remaining NA values after conversion
sum(is.na(ibs_mat))

### Creating annotation file
sample_ids <- dist_jc[-1, 1]
#cleaned_ids <- gsub("unclassified-kra_", "", sample_ids) #beware that the prefix might changed
cleaned_ids <- gsub("_bbmap", "", sample_ids) #beware that the suffix might changed
species_names <- rep("Oedulis", length(sample_ids))
#cleaned_ids <- sample_ids
### If you already have an tailored annotation file 
#annot_ed <- read_tsv(file = "../../Skmer-2/with_fam_nas/geno_fam_annotation.tsv" )

annotation_df <- data.frame(sample_id = sample_ids, cleaned_id = cleaned_ids, species = species_names)
#write.csv(annotation_df, "~/Desktop/GitHub/Gskimming/00_data/Skmer/Oedulis/RawCov/echarvel_oedulis_skimpip_results_17.08_nonrandlib_tresh0.2_0.2/17.08_nonrandlib_tresh0.2_0.2_oedulis_skmer1-dist-mat_annot.csv", row.names = FALSE)

### Link the master list to retrieve supplementary info
master_list<- read_csv("~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/oedulis_sra.csv")
str(annotation_df)
str(master_list)
# Perform the merge
annot_df_final <- annotation_df %>%
  left_join(master_list, by = c("cleaned_id" = "Run"))

### identify run_accession that are not present in annot_df_final
missing_ind <- master_list %>% 
  filter(!Run %in% annot_df_final$cleaned_id)
for (accession in missing_ind$Run) {
  cat(accession, "\n")
}

# Renaming the columns
annot_modern<- merge(annotation_df, annot_df_final, by = "cleaned_id")
#write.csv(annot_modern, "../../Gskimming/01_infofiles/ClupeaAtmore/ClupeaModern_annot.csv", row.names = FALSE)
fst_annot <- annot_df_final %>% dplyr::select(sample_id, `Sample Name`)
str(fst_annot)
# Assuming your data frame is named 'fst_annot'
fst_annot <- fst_annot %>%
  mutate(`Sample Name` = substr(`Sample Name`, 1, 4))

colnames(fst_annot) <- c("genome", "family")
# fst_annot$family <- gsub(",", "", fst_annot$family) # Replace commas with underscores
# fst_annot$family <- gsub(" ", "", fst_annot$family)  # Remove spaces

fst_annot$genome
# Writing the data frame to a TSV file
#write.table(fst_annot, "~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/17.08_nonrandlib_tresh0.2_0.2_oedulis_skmer1_Fstannot.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Performing PCoA with initial annotation
mds_plot <- PCoA(ibs_mat, annot_df_final$cleaned_id, annot_df_final$Locality,3, 1, 2, show.ellipse = TRUE, show.label = FALSE)
mds_plot
#ggsave(mds_plot, file = "~/Desktop/GitHub/Gskimming/02_figures/Oedulis/MDS/Kraken_n_consult_preprocessedOedulis_mds_skmer1_n214_pop_PC1PC2.png",scale = 1, dpi = 600)


### Extract outliers values
str(pcoa_table)
# Load necessary library
str(ibs_mat)

# Function to find outliers
find_outliers <- function(df, column) {
  Q1 <- quantile(df[[column]], 0.25)
  Q3 <- quantile(df[[column]], 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  outliers <- df[df[[column]] < lower_bound | df[[column]] > upper_bound, ]
  return(outliers)
}

# Applying the function to each dist_ column
outliers_dist_1 <- find_outliers(pcoa_table, 'dist_1')
outliers_dist_2 <- find_outliers(pcoa_table, 'dist_2')
outliers_dist_3 <- find_outliers(pcoa_table, 'dist_3')

# Combining the outliers from all dist_ columns
combined_outliers <- unique(rbind(outliers_dist_1, outliers_dist_2))
outlier_table <- combined_outliers[, c('individual', 'population')]
print(outlier_table)




# Extract the individual IDs from the outliers dataframe
outlier_ids <- outliers_dist_1$individual

# Assuming 'ibs_mat' is your original matrix with individual IDs as row and column names
# Remove the rows and columns in 'ibs_mat' corresponding to the outlier IDs

#mds_plot_ecot <- PCoA(ibs_mat, annot_df_final$cleaned_id, annot_df_final$sample_description,4, 1, 2, show.ellipse = FALSE, show.label = TRUE)
#ggsave(mds_plot_ecot, file = "~/Desktop/GitHub/Gskimming/02_figures/ClupeaAtmore/4x/Skmer/4x_ClupeaPCA_n45_ecotype_JCcorr_Skmer_rawcoverage_PC1PC2.png",scale = 1, dpi = 600)


# Remove outliers from pcoa_table
pcoa_table_clean <- pcoa_table %>%
  anti_join(combined_outliers, by = "individual")

# View the cleaned dataframe
pcoa_table_clean


#write.csv(pcoa_table, "~/Desktop/GitHub/Gskimming/01_infofiles/DBScan/4x_jc-24.02-pcoa_table-mat_4x_Clupea.csv", row.names = FALSE)
#write.csv(pcoa_table, "~/Desktop/GitHub/Gskimming/00_data/Skmer/Clupea/RawCov/MDS/pcoa_table-mat_40pc_4x_jc-24.02.24_Clupea.csv", row.names = FALSE)

pcoa_table
pcoa_table_genome_wide_joined_clean <- pcoa_table_clean %>%
  dplyr::select(1:5) %>%
  left_join(annot_df_final, by=c("individual"="cleaned_id", "population"="geo_loc_name"))
mds_plot_sum <-pcoa_table_genome_wide_joined_clean %>%
  ggplot(aes(x=dist_1, y=dist_2)) +
  geom_point(data=dplyr::select(pcoa_table_genome_wide_joined_clean, -population), size = 0.1, color="grey") +
  geom_point(size=1, mapping = aes(color=population)) +
  facet_wrap(~population) +
  theme_cowplot() +
  theme(axis.text = element_blank())
mds_plot_sum
#ggsave(mds_plot_sum, file = "~/Desktop/GitHub/Gskimming/02_figures/Oedulis/MDS/Oedulis_mds_percountry_skmer1_n214_pop_PC1PC2.png",scale = 1, dpi = 600)

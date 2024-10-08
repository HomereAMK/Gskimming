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
dist_jc <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Skmer2/Clupea/4x_respect_oct24/ref-dist-mat.txt", col_names = FALSE) %>%
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
cleaned_ids <- gsub("unclassified-kra_unclassified-seq_", "", sample_ids) #beware that the prefix might changed
cleaned_ids <- gsub("__merged", "", cleaned_ids) #beware that the suffix might changed
species_names <- rep("Clupea", length(sample_ids))
cleaned_ids
### If you already have an tailored annotation file
#annot_ed <- read_tsv(file = "../../Skmer-2/with_fam_nas/geno_fam_annotation.tsv" )

annotation_df <- data.frame(sample_id = sample_ids, cleaned_id = cleaned_ids, species = species_names)
#write.csv(annotation_df, "../../Skmer-2/data/ClupeaAtmore3feb/03.02-dist-mat_annot.csv", row.names = FALSE)

### Link the master list to retrieve supplementary info
master_list<- read_csv("~/Desktop/GitHub/Gskimming/Metadata/ClupeaAtmore/Master_list_ClupeaAtmore.csv")
annot_df_final <- merge(annotation_df, master_list, by.x = "cleaned_id", by.y = "run_accession", all.x = TRUE)


### identify run_accession that are not present in annot_df_final
missing_ind <- master_list %>% 
  filter(!run_accession %in% annot_df_final$cleaned_id)
for (accession in missing_ind$run_accession) {
  cat(accession, "\n")
}

# Renaming the columns
annot_modern<- merge(annotation_df, annot_df_final, by = "cleaned_id")
#write.csv(annot_modern, "../../Gskimming/01_infofiles/ClupeaAtmore/ClupeaModern_annot.csv", row.names = FALSE)
fst_annot <- annot_modern %>% dplyr::select(sample_id.x, population)


colnames(fst_annot) <- c("genome", "family")
# fst_annot$family <- gsub(",", "", fst_annot$family) # Replace commas with underscores
# fst_annot$family <- gsub(" ", "", fst_annot$family)  # Remove spaces

fst_annot$genome
# Writing the data frame to a TSV file
write.table(fst_annot, "~/Desktop/GitHub/Gskimming/01_infofiles/ClupeaAtmore/Fst/Skmer2/4x_respect_skmer2_02.10-dist-mat_Clupea_Fstannot.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Performing PCoA with initial annotation
mds_plot <- PCoA(ibs_mat, annot_df_final$cleaned_id, annot_df_final$population,3, 1, 2, show.ellipse = FALSE, show.label = TRUE)
mds_plot
ggsave(mds_plot, file = "~/Desktop/GitHub/Gskimming/02_figures/ClupeaAtmore/Mappedreads/Skmer/MDS/Mapped_reads_skmer_rawcov_ClupeaPCA_n42_pop_PC1PC2.png",scale = 1, dpi = 600)

str(ibs_mat)



#### Outliers removal ####
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


pcoa_table_genome_wide_joined_clean <- pcoa_table_clean %>%
  dplyr::select(1:5) %>%
  left_join(annot_df_final, by=c("individual"="cleaned_id", "population"="population"))
mds_plot_sum <-pcoa_table_genome_wide_joined_clean %>%
  ggplot(aes(x=dist_1, y=dist_2)) +
  geom_point(data=dplyr::select(pcoa_table_genome_wide_joined_clean, -population), size = 0.1, color="grey") +
  geom_point(size=1, mapping = aes(color=population)) +
  facet_wrap(~population) +
  theme_cowplot() +
  theme(axis.text = element_blank())
mds_plot_sum



# Create a filtered data frame excluding specified individuals
filtered_pcoa_table <- pcoa_table %>%
  filter(!individual %in% c("ERR9709251", "ERR9709230", "ERR9709203", "ERR9709248"))
# Perform PCA
pca_result <- prcomp(filtered_pcoa_table[, c("dist_1", "dist_2", "dist_3")], center = TRUE, scale. = TRUE)
# Extract PCA scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$individual <- filtered_pcoa_table$individual
pca_scores$population <- filtered_pcoa_table$population
# Plot PCA results
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = population)) +
  geom_point(size = 2) +
  #geom_text(aes(label = individual), size = 2, vjust = 1.5, hjust = 0.5) +
  facet_wrap(~ population) +
  theme_cowplot() +
  theme(axis.text = element_blank())

# Print the plot
print(pca_plot)




#PCoA_label(ibs_mat, annotation_df$cleaned_id, annotation_df$species,3, 1, 3, show.ellipse = F, show.label = F)

### Extract outliers values
str(pcoa_table)
# Load necessary library

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
combined_outliers <- unique(rbind(outliers_dist_1, outliers_dist_2, outliers_dist_3))
outlier_table <- combined_outliers[, c('individual', 'population')]
print(outlier_table)








### Reading metadata and updating annotation_df
metadata_clupea <- read_tsv("~/Desktop/GitHub/Gskimming/Metadata/filereport_read_run_PRJEB52723_tsv (5).txt", col_names = TRUE)
annotation_df <- merge(annotation_df, metadata_clupea, by.x = "cleaned_id", by.y = "run_accession")

### Adding 'type' column based on sample_alias
# annotation_df <- annotation_df %>%
#   mutate(
#     type = case_when(
#       substr(sample_alias, 1, 8) == "BAL_22_M" ~ "modern",
#       substr(sample_alias, 1, 8) == "BAL_22_H" ~ "ancient",
#       TRUE ~ "unknown"
#     )
#   )

### Performing PCoA with updated annotation and type information
mds_type <- PCoA(ibs_mat, annot_df_final$cleaned_id, annot_df_final$tissue_type, 3, 1, 2, show.ellipse = TRUE, show.label = FALSE)
#ggsave(mds_type, file = "~/Desktop/GitHub/Gskimming/02_figures/ClupeaAtmore/4x/Skmer/4x_ClupeaPCA_n45_tissue_JCcorr_Skmer_rawcoverage_PC1PC2.png",scale = 1, dpi = 600)

### Reading skims stats file
# stat_clupea <- read_tsv("../01_infofiles/ClupeaAtmore/Stats/stats-postprocess_03.02.csv", col_names = FALSE)
# str(stat_clupea)
# stat_clupea_df <- stat_clupea %>%
#   pivot_wider(names_from = X2, values_from = X3) %>%
#   rename(
#     coverage = `coverage`,
#     genome_length = `genome_length`,
#     read_length = `read_length`,
#     error_rate = `error_rate`
#   )%>%
#   mutate(X1 = str_extract(X1, "[^/]+$")) %>%
#   mutate(X1 = gsub("unclassified-kra_", "",X1))%>%
#   mutate(X1 = gsub("_", "",X1))
# annotation_df_tot <- merge(annot_df_final, stat_clupea_df, by.x = "cleaned_id", by.y = "X1")
#write.csv(annotation_df_tot, "../../Gskimming/01_infofiles/ClupeaAtmore/Total_ClupeaModern_annot.csv", row.names = FALSE)


annotation_df_tot <- read_csv("~/Desktop/GitHub/Gskimming/01_infofiles/ClupeaAtmore/Total_ClupeaModern_annot.csv")
# Range of coverage by adding CoverageCategory ~
annotation_df_tot$CoverageCategory <- ifelse(annotation_df_tot$coverage <= 0.7, "< 1x",
                                             ifelse(annotation_df_tot$coverage <= 2, "< 2x",
                                                    ifelse(annotation_df_tot$coverage <= 5, "< 5x",
                                                           ifelse(annotation_df_tot$coverage <= 10, "< 10x",
                                                                  ifelse(annotation_df_tot$coverage <= 20, "< 20x",
                                                                         ifelse(annotation_df_tot$coverage <= 30, "< 30x",
                                                                                ifelse(annotation_df_tot$coverage <= 50, "< 50x", "Above 50x")))))))




annotation_df_tot$error_rateCategory <- ifelse(annotation_df_tot$error_rate <= 0.005, "< 0.5%",
                                               ifelse(annotation_df_tot$error_rate <= 0.01, "< 1%",
                                                      ifelse(annotation_df_tot$error_rate <= 0.05, "< 5%", "error rate dangerously high")))


annotation_df_tot$genome_lengthCategory <- ifelse(annotation_df_tot$genome_length <= 7e+08, "< 7e+08 bp",
                                               ifelse(annotation_df_tot$genome_length <= 8e+08, "< 8e+08 bp",
                                                      ifelse(annotation_df_tot$genome_length <= 9e+08, "< 9e+08 bp",
                                                             ifelse(annotation_df_tot$genome_length <= 10e+08, "< 10e+08 bp", "genome length dangerously high"))))




### 
PCoA(ibs_mat, annotation_df_tot$cleaned_id, annotation_df_tot$CoverageCategory, 3, 1, 2, show.ellipse = FALSE, show.label = FALSE)
PCoA(ibs_mat, annotation_df_tot$cleaned_id, annotation_df_tot$error_rateCategory, 3, 1, 2, show.ellipse = FALSE, show.label = FALSE)
#ggsave(mds_errate, file = "~/Desktop/GitHub/Gskimming/02_figures/ClupeaAtmore/4x/Skmer/4x_ClupeaPCA_n45_errorrate_JCcorr_Skmer_rawcoverage_PC1PC2.png",scale = 1, dpi = 600)

PCoA(ibs_mat, annotation_df_tot$cleaned_id, annotation_df_tot$genome_lengthCategory, 3, 1, 2, show.ellipse = FALSE, show.label = FALSE)
#ggsave(mds_geno, file = "~/Desktop/GitHub/Gskimming/02_figures/ClupeaAtmore/4x/Skmer/4x_ClupeaPCA_n45_genomelength_JCcorr_Skmer_rawcoverage_PC1PC2.png",scale = 1, dpi = 600)








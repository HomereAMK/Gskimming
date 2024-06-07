### Clearing environment and setting working directory
rm(list = ls(all = TRUE))

# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)


D1 <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Skmer/Clupea/4x/MDS/subsampled_4.txt", col_names = FALSE) %>%
  dplyr::select(1:nrow(.)) %>% 
  as.matrix()

sample_ids_D1 <- D1[-1, 1]
D1 <- D1[-1, -1]
D1 <- as.data.frame(as.matrix(D1))
summary(D1)

D2 <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Skmer2/Clupea/RawCov/MDS/skmer2_Ref-09.05-dist-mat.txt", col_names = FALSE) %>%
  dplyr::select(1:nrow(.)) %>%
  as.matrix() 
sample_ids_D2 <- D2[-1, 1]

D2 <- D2[-1, -1]
D2 <- as.data.frame(as.matrix(D2))
summary(D2)

# Create a data frame with lower triangular elements of both matrices
df <- data.frame(Skmer = as.numeric(D1[lower.tri(D1)]), Skmer2 = as.numeric(D2[lower.tri(D2)]))
lower_tri_indices <- which(lower.tri(D1), arr.ind = TRUE)

# Display summary of the data frame
summary(df)
ggplot(df,aes(x=Skmer,y=Skmer2)) + geom_point()
# Add the row names corresponding to these indices
df$skmerx1 <- sample_ids_D1[lower_tri_indices[, 1]]
df$skmerx2 <- sample_ids_D2[lower_tri_indices[, 2]] 
# Remove the specified substring
df$skmerx1 <- sub("unclassified-kra_unclassified-seq_", "", df$skmerx1)
df$skmerx2 <- sub("unclassified-kra_unclassified-seq_", "", df$skmerx2)
df$skmerx1 <- sub("__merged", "", df$skmerx1)
df$skmerx2 <- sub("__merged", "", df$skmerx2)

str(df)


### Creating annotation file
cleaned_ids <- gsub("unclassified-kra_unclassified-seq_", "", sample_ids_D1) #beware that the prefix might changed
cleaned_ids <- gsub("__merged", "", cleaned_ids) #beware that the suffix might changed
species_names <- rep("Clupea", length(sample_ids_D1))
cleaned_ids
annotation_df <- data.frame(sample_id = sample_ids_D1, cleaned_id = cleaned_ids, species = species_names)
### Link the master list to retrieve supplementary info
master_list<- read_csv("~/Desktop/GitHub/Gskimming/Metadata/ClupeaAtmore/Master_list_ClupeaAtmore.csv")
annot_df_final <- merge(annotation_df, master_list, by.x = "cleaned_id", by.y = "run_accession", all.x = TRUE)





# Load required library
library(dplyr)

# Assuming your data frames are named df and master_list
# Step 1: Create a lookup table for populations based on run_accession
lookup_table <- annot_df_final %>% dplyr::select(cleaned_id, population)
# Step 2: Join df with the lookup_table to get population for skmerx2
df_with_population <- df %>%
  left_join(lookup_table, by = c("skmerx2" = "cleaned_id")) %>%
  rename(population_skmerx2 = population)
# Step 3: Join df_with_population with the lookup_table again to get population for skmerx1
df_with_population <- df_with_population %>%
  left_join(lookup_table, by = c("skmerx1" = "cleaned_id")) %>%
  rename(population_skmerx1 = population)
# Step 4: Create the comp variable by concatenating the population values
df_with_population <- df_with_population %>%
  mutate(comp = paste(population_skmerx2, population_skmerx1, sep = "_"))
# View the resulting data frame
print(df_with_population)
df_with_population$comp

ggplot(df_with_population,aes(x=Skmer,y=Skmer2, color=comp)) + geom_point()


# Step 1: Modify the comp variable to replace identical population values with "same"
df_with_population <- df_with_population %>%
  mutate(comp = ifelse(substr(comp, 1, 4) == substr(comp, 6, 9), 
                       "same", comp))

# Step 2: Ensure that combinations like "IDEO_RISO" and "RISO_IDEO" are treated the same
df_with_population <- df_with_population %>%
  mutate(comp = ifelse(substr(comp, 1, 4) > substr(comp, 6, 9), 
                       paste(substr(comp, 6, 9), substr(comp, 1, 4), sep = "_"), 
                       comp))

ggplot(df_with_population,aes(x=Skmer,y=Skmer2, color=comp)) + geom_point()


df_with_population$comp[ df_with_population$comp != "_same" ] <- "Mixed"
ggplot(df_with_population,aes(x=Skmer,y=Skmer2, color=comp)) + geom_point() + facet_grid(comp~.)



##density plot
# A single scatter plot
scatter <- ggplot(df_with_population,aes(x=Skmer,y=Skmer2, color=comp)) + geom_point() + 
  theme(legend.position=c(1,1),legend.justification=c(1,1)) 
#top
plot_top <- ggplot( df_with_population, aes(x=Skmer,fill=comp)) + geom_density(alpha=0.5) + theme(legend.position="none")
plot_side <- ggplot( df_with_population, aes(x=Skmer2, fill=comp)) + geom_density(alpha=0.5) + coord_flip() + theme(legend.position="none")

empty <- ggplot() + theme_void()

library(gridExtra)
grid.arrange(plot_top, empty, scatter, plot_side, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))


write_tsv(annot_df_final, "~/Desktop/GitHub/Gskimming/01_infofiles/Clupea_n42_13var.tsv")



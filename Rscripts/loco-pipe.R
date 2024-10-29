### Clearing environment and setting working directory
rm(list = ls(all = TRUE))

library(readr)
library(dplyr)

#### Herring ####
# Read the TSV file
data <- read_csv("~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/oedulis_sra.csv", col_names = TRUE)

# # Select the required columns and add an empty column between the first and second column
# selected_data <- data %>%
#   dplyr::select(sample_name = sample_name, empty_var = sample_id, species = sample_description, population) %>%
#   mutate(empty_var = "")
# 
# #empty_var should be called "bam" and values should be for each row the combination of /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_realigned/" then $sample_name value for each row then "_minq20.nocig.realigned.bam"

# Read the TSV file
data <- read_tsv("~/Desktop/GitHub/Gskimming/01_infofiles/Clupea_n42_13var.tsv", col_names = TRUE)

# Select the required columns and create the new bam column with the specified values
selected_data <- data %>%
  dplyr::select(sample_name = sample_name, species=sample_description, population) %>%
  mutate(bam = paste0("/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned/", sample_name, "_minq20.nocig.realigned.bam"))

# View the resulting dataframe
print(selected_data)
write_tsv(selected_data, "~/Desktop/GitHub/Gskimming/01_infofiles/sample_table_Herring_0.5x_loco.tsv")

#raw cov
selected_data_raw <- data %>%
  dplyr::select(sample_name = sample_name, species=sample_description, population) %>%
  mutate(bam = paste0("/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned/unclassified-seq_", sample_name, ".nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned.bam"))
write_tsv(selected_data_raw, "~/Desktop/GitHub/Gskimming/01_infofiles/sample_table_Herring_raw_loco.tsv")









#### oedulis ####
### Clearing environment and setting working directory
rm(list = ls(all = TRUE))

data <- read_csv("~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/oedulis_sra.csv", col_names = TRUE)
Stein_bam <- read_table("~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/oedulis_SteinMortensen.tsv", col_names = FALSE) %>% 
              rename(bam = X1)
str(Stein_bam)
str(data)
# Select the required columns and add an empty column between the first and second column
selected_data <- data %>%
  dplyr::select( sample_name = Run, species = geo_loc_name , population =Locality) %>%
  mutate(bam = paste0("/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/Stein_Mortensen_fqs/realigned/", sample_name, "_minq20.nocig.realigned.bam")) %>%
  semi_join(Stein_bam, by = "bam")
str(selected_data)

selected_data <- selected_data %>%
  mutate(population = if_else(population == "Aga_Bømlo", "AGAB", toupper(substr(population, 1, 4))))
write_tsv(selected_data, "~/Desktop/GitHub/Gskimming/01_infofiles/SteinOedulis_loco.tsv")
chr_table <- read_tsv("~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/oedulis_chr_table.tsv", col_names = TRUE)
str(chr_table)


# Read the original TSV file
chr_table <- read_tsv("~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/oedulis_chr_table.tsv", col_names = TRUE)

# Select the required columns and remove column names
chr_table_selected <- chr_table %>%
  dplyr::select(`RefSeq seq accession`, `Chromosome name`)


# Write the modified data back to the original file without column names
write_tsv(chr_table_selected, "~/Desktop/GitHub/Gskimming/01_infofiles/oedulis/oedulis_chr_table.tsv", col_names = FALSE)

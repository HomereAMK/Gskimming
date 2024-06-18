### Clearing environment and setting working directory
rm(list = ls(all = TRUE))

library(readr)
library(dplyr)

# Read the TSV file
data <- read_tsv("~/Desktop/GitHub/Gskimming/01_infofiles/Clupea_n42_13var.tsv", col_names = TRUE)

# Select the required columns and add an empty column between the first and second column
selected_data <- data %>%
  dplyr::select(sample_name = sample_name, empty_var = sample_id, species = sample_description, population) %>%
  mutate(empty_var = "")

#empty_var should be called "bam" and values should be for each row the combination of /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_realigned/" then $sample_name value for each row then "_minq20.nocig.realigned.bam"


library(readr)
library(dplyr)

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

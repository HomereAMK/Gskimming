### Clearing environment and setting working directory
rm(list = ls(all = TRUE))

library(readr)
library(dplyr)

# Read the TSV file
data <- read_tsv("~/Desktop/GitHub/Gskimming/01_infofiles/Clupea_n42_13var.tsv", col_names = TRUE)

# Select the required columns and add an empty column between the first and second column
selected_data <- data %>%
  dplyr::select(sample_name = cleaned_id, empty_var = sample_id, species = sample_description, population) %>%
  mutate(empty_var = "")

#empty_var should be called "bam" and values should be for each row the combination of /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_realigned/" then $sample_name value for each row then "_minq20.nocig.realigned.bam"


library(readr)
library(dplyr)

# Read the TSV file
data <- read_tsv("~/Desktop/GitHub/Gskimming/01_infofiles/Clupea_n42_13var.tsv", col_names = TRUE)

# Select the required columns and create the new bam column with the specified values
selected_data <- data %>%
  dplyr::select(sample_name = cleaned_id, species, population) %>%
  mutate(bam = paste0("/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_realigned/", sample_name, "_minq20.nocig.realigned.bam"))

# View the resulting dataframe
print(selected_data)
write_tsv(selected_data, "~/Desktop/GitHub/Gskimming/01_infofiles/sample_table_Herring_4x_loco.tsv")

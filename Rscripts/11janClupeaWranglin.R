### The BEGINNING ~~~~
##
# ~
# Cleans the environment ~ 
rm(list=ls())

#install.package("pacman")
# Loads required packages ~
pacman::p_load(pheatmap, tidyverse, reshape2)

help("read_csv")
# Reading the TSV file into a data frame
data1 <- read_tsv("~/Desktop/GitHub/Gskimming/Metadata/Clupea/11janENA_clupeaAtmore.tsv")

# Reading the CSV file into a data frame
data2 <- read_delim("~/Desktop/GitHub/Gskimming/Metadata/Clupea/11janModern_Inds_Clupea.csv", delim = ";")

# Merging the data frames on the 'sample_title' variable
merged_data <- merge(data1, data2, by = "sample_title")

# Viewing the first few rows of the merged data frame
head(merged_data)
write.csv(merged_data, "~/Desktop/GitHub/Gskimming/Metadata/Clupea/Master_list_ClupeaAtmore.csv", row.names = FALSE)

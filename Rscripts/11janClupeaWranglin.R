### The BEGINNING ~~~~
##
# ~
# Cleans the environment ~ 
rm(list=ls())

#install.package("pacman")
# Loads required packages ~
pacman::p_load(pheatmap, tidyverse, reshape2)

#### ClupeaAtmore ####
# Reading the TSV file into a data frame
data1 <- read_tsv("~/Desktop/GitHub/Gskimming/Metadata/ClupeaAtmore/11janENA_clupeaAtmore.tsv")

# Reading the CSV file into a data frame
data2 <- read_delim("~/Desktop/GitHub/Gskimming/Metadata/ClupeaAtmore/11janModern_Inds_Clupea.csv", delim = ";")

# Merging the data frames on the 'sample_title' variable
merged_data <- merge(data1, data2, by = "sample_title")


# Extracting the first four letters and converting to uppercase
merged_data$population <- toupper(substr(merged_data$`geographic location (region and locality)`, 1, 4))

# Now, your merged_data will have an additional column 'new_variable' with the desired values
# Viewing the first few rows of the merged data frame
str(merged_data)
write.csv(merged_data, "~/Desktop/GitHub/Gskimming/Metadata/ClupeaAtmore/Master_list_ClupeaAtmore.csv", row.names = FALSE)






#### ClupeaHan ####
# Reading the TSV file into a data frame
Han_raw <- read_tsv("~/Desktop/GitHub/Gskimming/Metadata/ClupeaHan/filereport_read_run_PRJNA642736_ClupeaHan.tsv")

# Reading the CSV file into a data frame
data2 <- read_delim("~/Desktop/GitHub/Gskimming/Metadata/ClupeaAtmore/11janModern_Inds_Clupea.csv", delim = ";")


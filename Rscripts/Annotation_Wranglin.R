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



#### Magpie ####
magpie_annot <- read_tsv("~/Desktop/GitHub/Gskimming/01_infofiles/Magpie/Magpie_bam_list.txt", col_names = FALSE)
magpie_annot$sample_id <- magpie_annot$X1
magpie_annot$sample_id <- gsub("/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/realigned/", "", magpie_annot$sample_id)
magpie_annot$sample_id <- gsub(".nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned.bam", "", magpie_annot$sample_id)
# create the pop var
magpie_annot <- magpie_annot %>%
  mutate(population = sub("_.*", "", sample_id))
write.csv(magpie_annot,"~/Desktop/GitHub/Gskimming/01_infofiles/Magpie/Magpie_annot.csv")



#### Ragweed #####

# Reading the TSV file into a data frame
data1 <- read_delim("~/Desktop/GitHub/Gskimming/01_infofiles/Ragweed/sciadv.abo5115_table_s1.csv", col_names = TRUE)

# Reading the CSV file into a data frame
data2 <- read_table("~/Desktop/GitHub/Gskimming/01_infofiles/Ragweed/plant-modern-filenames.txt", col_names = FALSE)
colnames(data2) <- c("SRA accession code")


# Merging the data frames on the 'sample_title' variable
merged_data <- merge(data1, data2, by = "SRA accession code")
write.csv(merged_data, "~/Desktop/GitHub/Gskimming/01_infofiles/Ragweed/Master_list_Ragweed_modern.csv", row.names = FALSE)

### The BEGINNING ~~~~~
##
# ~ Plot  Heterozygosity compute by angsd | from Homère J. Alves Monteiro and Eduardo Charvel modify from Nicolas R. Lou

# R 4.3.2 GUI 1.80 Big Sur ARM build (8281)

### Clearing environment and setting working directory
rm(list = ls(all = TRUE))
library(tidyverse)
library(cowplot)
library(RcppCNPy)
library(knitr)
master_annot <- read_csv("~/Desktop/Github/Gskimming/01_infofiles/ClupeaAtmore/Total_ClupeaModern_annot.csv")
# wrangling
master_annot$sample_id <- gsub("unclassified-kra", "unclassified-seq", master_annot$sample_id)
master_annot$sample_id <- sub("_$", "", master_annot$sample_id)


# Excluding transitions -> Initialize an empty data frame
all_theta_notrans <- data.frame()

for (i in 1:nrow(master_annot)){
  sample_seq_id <- master_annot$sample_id[i]
  sample_id <- master_annot$cleaned_id[i]
  population <- master_annot$population[i]
  data_type <- master_annot$tissue_type[i]
  group <- master_annot$sample_description[i]
  path <- str_c("~/Desktop/Github/Gskimming/00_data/Angsd/Clupea/RawCov/Theta/noTrans/", sample_seq_id,  ".nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned_mindp2_maxdp800_minq33_minmapq20")
  
  # Check if the file exists before trying to read it
  if (file.exists(str_c(path, ".average_thetas.pestPG"))) {
    theta_notrans <- read_tsv(str_c(path, ".average_thetas.pestPG")) %>%
      janitor::clean_names() %>%
      dplyr::select(chr, t_w, n_sites) %>%
      mutate(sample_id=sample_id, population=population, group=group, type="Excluding transitions")
    
    # Append the data from this iteration to the all_theta_notrans data frame
    all_theta_notrans <- rbind(all_theta_notrans, theta_notrans)
  } else {
    message("File not found for sample_seq_id: ", sample_seq_id)
  }
}

# With transitions -> Initialize an empty data frame
all_theta_trans <- data.frame()

for (i in 1:nrow(master_annot)){
  sample_seq_id <- master_annot$sample_id[i]
  sample_id <- master_annot$cleaned_id[i]
  population <- master_annot$population[i]
  data_type <- master_annot$tissue_type[i]
  group <- master_annot$sample_description[i]
  path <- str_c("~/Desktop/Github/Gskimming/00_data/Angsd/Clupea/RawCov/Theta/withTrans/", sample_seq_id,  ".nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned_mindp2_maxdp800_minq33_minmapq20_withTrans")
  
  # Check if the file exists before trying to read it
  if (file.exists(str_c(path, ".average_thetas.pestPG"))) {
    theta_trans <- read_tsv(str_c(path, ".average_thetas.pestPG")) %>%
      janitor::clean_names() %>%
      dplyr::select(chr, t_w, n_sites) %>%
      mutate(sample_id=sample_id, population=population, group=group, type="Including transitions")
    
    # Append the data from this iteration to the all_theta_notrans data frame
    all_theta_trans <- rbind(all_theta_trans, theta_trans)
  } else {
    message("File not found for sample_seq_id: ", sample_seq_id)
  }
}



## Calculate average heterozygosity (with filtering of inversions if needed)
het <- all_theta_notrans %>%
  filter(! is.na(chr)) %>% 
  #filter(! chr %in% c("LR535869.1")) %>% filter out some chr if needed
  group_by(sample_id, population, group, type) %>%
  summarise(sum_t_w=sum(t_w), sum_n_sites=sum(n_sites), heterozygosity=sum_t_w/sum_n_sites) 
set.seed(42)
het %>%ggplot(aes(x=population, y=heterozygosity, color=population)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0) +
  facet_wrap(~type) +
  coord_flip() +
  theme_cowplot()+
  theme(legend.position = "none")
ggsave(het, file = "~/Desktop/GitHub/Gskimming/02_figures/ClupeaAtmore/RawCov/Angsd/Theta/Het_Clupea_n45_Invfree_minq20_minq20.nocig.realigned_mindp2_maxdp800_minq33_minmapq20.average_thetas.png",scale = 1, dpi = 600)


#mergin with and excluding trans and replotting
theta_final <- bind_rows(all_theta_notrans, all_theta_trans)
het_tot <- theta_final %>%
  filter(! is.na(chr)) %>% 
  #filter(! chr %in% c("LG01", "LG02", "LG07", "LG12")) %>%
  group_by(sample_id, population, group, type) %>%
  summarise(sum_t_w=sum(t_w), sum_n_sites=sum(n_sites), heterozygosity=sum_t_w/sum_n_sites) %>%
  mutate(type=factor(type,levels=c("Including transitions", "Excluding transitions")))
set.seed(42)
het_tot %>%
  ggplot(aes(x=population, y=heterozygosity, color=population)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0) +
  facet_wrap(~type) +
  coord_flip() +
  theme_cowplot() +
  theme(legend.position = "none")
het
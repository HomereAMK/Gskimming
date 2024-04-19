### The BEGINNING ~~~~~
##
# ~ Plot read and raw bases counts before and after mapping with angsd| from Homère J. Alves Monteiro and Eduardo Charvel modify from Nicolas R. Lou

# R 4.3.2 GUI 1.80 Big Sur ARM build (8281)
pacman::p_load(tidyverse, cowplot, knitr, RcppCNPy, ggrepel, MetBrewer, ade4)

### Clearing environment and setting working directory
master_annot <- read_csv("~/Desktop/Github/Gskimming/01_infofiles/ClupeaAtmore/Total_ClupeaModern_annot.csv")
stats_clupea <- read_tsv("~/Desktop/GitHub/Gskimming/01_infofiles/ClupeaAtmore/Stats/Summary_modernClupea_stats_9apr23.txt", col_names = FALSE)
colnames(stats_clupea) <- c("cleaned_id", "rawreads", "rawbases", "adapterclipbases", "mappedbases", "dedupmappedbases", "realignedmappedbases")
stats_clupea <- subset(stats_clupea, select = -ncol(stats_clupea))

ann_stats_clupea <- merge(stats_clupea, master_annot, by.x = "cleaned_id")
write.csv(ann_stats_clupea, ""~/Desktop/GitHub/Gskimming/01_infofiles/ClupeaAtmore/Stats/ann_stats_clupea.csv")
select_relevant_rows <- function(x){
  dplyr::select(x, cleaned_id, population, rawbases, adapterclipbases, mappedbases, dedupmappedbases, realignedmappedbases)
}


sum_count <- function(x){
  summarize(x, rawbases=sum(as.numeric(rawbases)), 
            adapterclipbases=sum(as.numeric(adapterclipbases)), 
            mappedbases=sum(as.numeric(mappedbases)), 
            dedupmappedbases=sum(as.numeric(dedupmappedbases)),
            realignedmappedbases=sum(as.numeric(realignedmappedbases)))%>%
    ungroup()
}

# In the data
count_by_pop <- group_by(ann_stats_clupea, population) %>%
  sum_count() %>%
  gather(key="steps", value = base_count, 2:6) %>%
  mutate(steps=fct_reorder(steps,base_count, mean, .desc = T))
sums <- group_by(ann_stats_clupea, population) %>%
  summarise(sum=sum(as.numeric(dedupmappedbases)))
means <- group_by(ann_stats_clupea, population) %>%
  summarise(mean=mean(dedupmappedbases))
proportion_retained <-group_by(ann_stats_clupea, population) %>%
  summarise(proportion=mean(realignedmappedbases/rawbases))

# Note , adapterclipbases correspond to fastq files after consult decontamination.

#Total number of bases per individual, grouped by populations (in GB)
means <- group_by(ann_stats_clupea, population) %>%
  summarise(mean=mean(dedupmappedbases))
set.seed(1)
ggplot(ann_stats_clupea, aes(x=population, y=dedupmappedbases/10^9)) +
  geom_jitter(aes(color=population), height=0) +
  geom_boxplot(alpha=0, outlier.colour=NA) +
  geom_label(data=means, aes(label = round(mean/10^9, 3), y = mean/10^9), alpha=1) +
  geom_point(data=subset(ann_stats_clupea, genome_length > 1000000000), aes(x=population, y=dedupmappedbases/10^9), color="red", size=4) +
  geom_text(data=subset(ann_stats_clupea, genome_length > 1000000000), aes(x=population, y=dedupmappedbases/10^9, label=cleaned_id), size=2, vjust=-1) +
  coord_flip() +
  theme_cowplot()

# Sum of bases retained in each step per population (in GB)
ggplot(count_by_pop) +
  geom_bar(mapping=aes(x=population, y=base_count/10^9, fill=steps), stat="identity", pos="identity", color="black", width = 0.8) +
  geom_label(data=sums, aes(label=round(sum/10^9,2), x=population, y=sum/10^9-0.8)) +
  coord_flip() +
  theme_cowplot()


# Percentage of bases retained after consult grouped by populations (in percent)
set.seed(1)
ggplot(ann_stats_clupea, aes(x=population, y=adapterclipbases/rawbases*100)) +
  geom_jitter(aes(color=population), height=0) +
  geom_boxplot(alpha=0, outlier.colour=NA) +
  geom_point(data=subset(ann_stats_clupea, genome_length > 1000000000), aes(x=population, y=adapterclipbases/10^9), color="red", size=4) +
  geom_text(data=subset(ann_stats_clupea, genome_length > 1000000000), aes(x=population, y=adapterclipbases/10^9, label=cleaned_id), size=2, vjust=-1) +
  coord_flip() +
  theme_cowplot()



# Percentage of bases retained after consult + (not kraken) +mapping + rm deduplicated reads + realign reads around indels, grouped by populations (in percent)
set.seed(1)
ggplot(ann_stats_clupea, aes(x=population, y=realignedmappedbases/rawbases*100)) +
  geom_jitter(aes(color=population), height=0) +
  geom_boxplot(alpha=0, outlier.colour=NA) +
  geom_label(data=proportion_retained, aes(label = round(proportion*100, 1), y = proportion*100), alpha=1) +
  geom_point(data=subset(ann_stats_clupea, genome_length > 1000000000), aes(x=population, y=realignedmappedbases/10^9), color="red", size=4) +
  geom_text(data=subset(ann_stats_clupea, genome_length > 1000000000), aes(x=population, y=realignedmappedbases/10^9, label=cleaned_id), size=2, vjust=-1) +
  coord_flip() +
  theme_cowplot()


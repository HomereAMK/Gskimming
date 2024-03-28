### The BEGINNING ~~~~~
##
# ~ Plots --per-position depth distribution in order to set up reasonable depth filters. | First written by Nicolas Lou with later modifications by Homère J. Alves Monteiro

rm(list = ls(all = TRUE))

pacman::p_load(tidyverse, cowplot, knitr, RcppCNPy, scales, MetBrewer, ade4)

setwd("~/Desktop/Scripts/EUostrea")
#### depth distribution obtained from ANGSD to establish depth filters for SNP calling. ####
depth_hist_angsd <- read_lines("~/Desktop/Scripts/EUostrea/02_data/SetAngsdFilters/Nlou_filtered_minq20_minmapq20_angsd0.929_htslib1.9.depthGlobal") %>%
  str_split(pattern = "\t") %>%
  .[[1]] %>%
  .[1:6001] %>%
  as.integer() %>%
  tibble(depth=0:6000, count=., method="angsd")
histo_distrib_2<- depth_hist_angsd %>%
  filter(depth>0, depth<6100) %>%
  ggplot(aes(x=depth, y=count, color=method)) +
  geom_freqpoly(stat = "identity") +
  theme_cowplot()
ggsave(histo_distrib_2, file = "Figures/SetAngsdFilters/per-position_depth_distrib.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()
# 923 is the mode of the peak using ANGSD
filter(depth_hist_angsd, depth>5) %>%
  group_by(method) %>%
  slice_max(order_by = count)

#### Find the standard deviation of the peak by fitting the right part of it to a normal distribution. #### 
sd_table <- NULL
# Calculate residual sum of squares by iterating through a range of possible standard devation values
for (standard_deviation in 5:95){
  x <- filter(depth_hist_angsd, depth>=923, depth<=8000) %>%
    mutate(temp=dnorm(depth, mean=923, sd=standard_deviation)) %>%
    mutate(count_fitted=temp*1225196/max(temp)) %>%
    summarise(rss=sum((count-count_fitted)^2)) %>%
    mutate(sd=standard_deviation)
  sd_table=bind_rows(sd_table, x)
}
# Find the standard deviation value with the lowest RSS
optimal_sd <-slice_min(sd_table, rss, n = 1)$sd
optimal_sd
#> optimal_sd
#[1] 80
sd_table %>%
  ggplot(aes(x=sd, y=rss)) +
  geom_line()
# fit right part or the distrib to normal distribution
fit_norm <- filter(depth_hist_angsd, depth>1, depth<=2000) %>%
  mutate(temp=dnorm(depth, mean=923, sd=optimal_sd)) %>%
  mutate(count_fitted=temp*1225196/max(temp)) %>%
  ggplot(aes(x=depth, y=count)) +
  geom_freqpoly(stat = "identity") +
  geom_freqpoly(aes(y=count_fitted), stat = "identity", color="blue") +
  theme_cowplot()
ggsave(fit_norm, file = "Figures/SetAngsdFilters/Fit_normal_per-position_depth_distrib.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()
#The fitted normal distribution doesn't match the empricial depth distribution quite well.

#### Use two standard deviations from the mode as depth cutoffs #### 
depth_hist_angsd <- depth_hist_angsd %>% drop_na()
max_depth <- 923 + 2*optimal_sd # =1121
sum(filter(depth_hist_angsd, depth>max_depth)$count)/sum(depth_hist_angsd$count)
# #sites
# x <- filter(depth_hist_angsd, depth>max_depth)
# count <- x$count
# a <- sum(count, na.rm = TRUE)
# 
# y <- depth_hist_angsd$count
# b <- sum(y, na.rm = TRUE)
# a/b
# [1] 0.1586364
#final mapped data
depth_hist_angsd <- read_lines("~/Desktop/Scripts/EUostrea/02_data/SetAngsdFilters/Nlou_filtered_minq20_minmapq20_angsd0.929_htslib1.9.depthGlobal") %>%
  str_split(pattern = "\t") %>%
  .[[1]] %>%
  .[1:2001] %>%
  as.integer() %>%
  tibble(depth=0:2000, count=., method="angsd")
as.double(sum(filter(depth_hist_angsd, depth>max_depth)$count*filter(depth_hist_angsd, depth>max_depth)$depth)/sum(depth_hist_angsd$count*depth_hist_angsd$depth))
## If a MaxDepth=1121 filter is used (shown below), ~11.7% of all sites and ~19,1% of the final mapped data will be lost

min_depth <- 923 - 2*optimal_sd
sum(filter(depth_hist_angsd, depth<min_depth)$count)/sum(depth_hist_angsd$count)
#[1] 0.4087959
sum(filter(depth_hist_angsd, depth<min_depth)$count*filter(depth_hist_angsd, depth<min_depth)$depth)/sum(depth_hist_angsd$count*depth_hist_angsd$depth)
#[1] 0.2023289
## If a MinDepth=725 filter is used (shown below), ~41% of all sites and ~20% of the final mapped data will be lost


## If these filters are used
cut_off_norm <-filter(depth_hist_angsd, depth>1, depth<=2000) %>%
  ggplot(aes(x=depth, y=count)) +
  geom_freqpoly(stat = "identity") +
  geom_vline(xintercept = c(max_depth,min_depth), color="red") +
  theme_cowplot()
ggsave(cut_off_norm, file = "Figures/SetAngsdFilters/cutoffMinMaxDepthFit_normal_per-position_depth_distrib.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()

# Way too much read loss. Difficulties to fit a normal distribution to the read depth distribution.


#### Settling for a MaxDepth=1200 (~6.2% of all sites and ~12% of the mapped data loss) and MinDepth=600 (~32% of all sites and ~12% of the final data mapped lost) #### 
max_depth <- 1200# =1121
sum(filter(depth_hist_angsd, depth>max_depth)$count)/sum(depth_hist_angsd$count)
as.double(sum(filter(depth_hist_angsd, depth>max_depth)$count*filter(depth_hist_angsd, depth>max_depth)$depth)/sum(depth_hist_angsd$count*depth_hist_angsd$depth))
min_depth <- 600
sum(filter(depth_hist_angsd, depth<min_depth)$count)/sum(depth_hist_angsd$count)
#[1] 0.4087959
sum(filter(depth_hist_angsd, depth<min_depth)$count*filter(depth_hist_angsd, depth<min_depth)$depth)/sum(depth_hist_angsd$count*depth_hist_angsd$depth)
#[1] 0.2094492
viz_cutoff <- filter(depth_hist_angsd, depth>1, depth<=2000) %>%
  ggplot(aes(x=depth, y=count)) +
  geom_freqpoly(stat = "identity") +
  geom_vline(xintercept = c(1200,600), color="purple") +
  labs(
    colour = "Number of cylinders",
    title = "MaxDepth=1200 and MinDepth=600",
    subtitle = "(~6.2% of all sites and ~12% of the mapped data loss and ~32% of all sites and ~12% of the final data mapped lost)")
ggsave(viz_cutoff, file = "Figures/SetAngsdFilters/cutoffMinMaxDepthFit_normal_per-position_depth_distrib.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()



#### Set up minimum individual filters #### 
#The distribution of numbers of individuals being present at each site is available from the mafs.gz file generated by ANGSD. Note that the mafs.gz file is very large and this step can be very memory consuming (~55GB for this dataset).
#On c2
nind_angsd <- data.table::fread("/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters/Nlou_filtered_minq20_minmapq20_angsd0.929_htslib1.9.mafs.gz") %>%
  count(nInd)
write_tsv(nind_angsd, "/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters/Nlou_filtered_minq20_minmapq20_angsd0.929_htslib1.9.presenceGlobal")

nind_angsd <- read_tsv("../Data/SetAngsdFilters/Nlou_filtered_minq20_minmapq20_angsd0.929_htslib1.9.presenceGlobal")
nind_angsd %>%
  slice_max(n)
# A tibble: 1 × 2
#nInd       n
#<dbl>   <dbl>
#  1   450 6116239

#The mode of the second peak is 450 and the sample size is 582 Let’s use 1/4 of the sample size as our minInd filter.
## minInd=145
nind_angsd %>%
  ggplot(aes(x=nInd, y=n)) +
  geom_freqpoly(stat = "identity") +
  geom_vline(xintercept = 582/4, color="purple") +
  theme_cowplot()
ggsave(nind_angsd, file = "Figures/SetAngsdFilters/minInd_cutoffdistrib_numbers_individualperSite.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()


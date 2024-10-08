### The BEGINNING ~~~~~
## 
# ~ Plots -- mapping stats Gskimming
rm(list = ls(all = TRUE))

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Load the data
cstats <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Mapping_stats/combined_3dataset_stats_mapping.tsv")
str(cstats)

# Add species variable based on sample names
cstats <- cstats %>%
  mutate(species = case_when(
    substr(Sample, 1, 3) == "ERR" ~ "Herring",
    substr(Sample, 1, 3) == "SRR" ~ "Oyster",
    substr(Sample, 1, 3) == "L68" ~ "Sandeel",
    TRUE ~ "Unknown"
  ))

# Pivot the data into a longer format
cstats_long <- cstats %>%
  pivot_longer(cols = c(Percentage_Mapped_Reads_BWA, Mapped_Percentage_Post_BWA_BBMAP, 
                        Mosdepth_Avg_Depth, Prop_g_cov_1,  Prop_g_cov_4),
               names_to = "Metric", values_to = "Value")

# Reset the label index
plot_index_labels <- 0

# Custom y-axis labels function
labels_fun <- function(x) {
  plot_index_labels <<- plot_index_labels + 1L
  switch(plot_index_labels,
         # For percentage values
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(x),  # Percentage_Mapped_Reads_BWA
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(x),  # Mapped_Percentage_Post_BWA_BBMAP
         
         # For coverage depth values (e.g., Mosdepth_Avg_Depth) add "X"
         scales::label_number(accuracy = 0.1, suffix = "X")(x),  # Mosdepth_Avg_Depth
         
         # For genomic coverage proportions (Prop_g_cov_1 and Prop_g_cov_4) add "%"
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(x),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(x)
  )
}

# Define ylabels for each Metric
ylabels <- c(
  "Percentage_Mapped_Reads_BWA" = "Percentage Mapped Reads BWA",
  "Mapped_Percentage_Post_BWA_BBMAP" = "Mapped Percentage Post BWA BBMAP",
  "Mosdepth_Avg_Depth" = "Mosdepth Average Depth (X)",
  "Prop_g_cov_1" = "Proportion Genomic Coverage 1",
  "Prop_g_cov_4" = "Proportion Genomic Coverage 4"
)

# Create the combined plot with violin plot and summary statistics
BSG_Combined <- ggplot() +
  geom_violin(data = cstats_long %>% 
                filter(!(Metric == "Mosdepth_Avg_Depth" & Value > 10)), # Exclude outliers for Mosdepth_Avg_Depth
              aes(x = species, y = Value), 
              fill = "#ffffff", colour = "#000000", show.legend = FALSE, alpha = .9, size = .3, width = .5) +
  stat_summary(data = cstats_long, aes(x = species, y = Value),  
               fun = mean, geom = "point", shape = 21, size = 2.5, alpha = 1, colour = "#000000", fill = "#df65b0") +
  facet_grid(Metric ~ ., scales = "free", labeller = labeller(Metric = ylabels)) + # Apply ylabels
  scale_y_continuous(labels = labels_fun) +
  theme(
    panel.background = element_rect(fill = "#ffffff"),
    panel.grid.major = element_line(color = "#ededed", linetype = "dashed", size = .00005),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.spacing.y = unit(1, "cm"),
    axis.line = element_line(colour = "#000000", size = .3),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "#000000", size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "#000000", size = 8, face = "bold"),
    axis.ticks.x = element_line(color = "#000000", size = .3),
    axis.ticks.y = element_line(color = "#000000", size = .3),
    strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
    strip.text = element_text(colour = "#000000", size = 6, face = "bold"),
    legend.position = "top",
    legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
    legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
    legend.key = element_rect(fill = NA),
    legend.background = element_blank()
  )

# Display the plot
print(BSG_Combined)

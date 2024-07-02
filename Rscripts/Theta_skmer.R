# Load necessary libraries
library(tidyverse)

# Theta Skmer1
data <- read_tsv("~/Desktop/GitHub/Gskimming/00_data/Skmer/theta_summary_jun24.tsv")

# Convert data from wide to long format
data_long <- data %>%
  pivot_longer(cols = -population, names_to = "coverage", values_to = "heterozygosity")

# Convert 'coverage' to a factor with a specific order
data_long$coverage <- factor(data_long$coverage, levels = c("raw_cov", "8x", "4x", "2x", "1x", "0.5x", "0.25x"))

# Plot the data
ggplot(data_long, aes(x = coverage, y = heterozygosity, group = population, color = population)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.0042, linetype = "dashed", color = "red") +
  labs( x = "Coverage", 
       y = "Heterozygosity") +
  theme_minimal()

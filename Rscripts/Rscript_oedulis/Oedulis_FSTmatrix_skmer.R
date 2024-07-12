### The BEGINNING ~~~~
##
# ~ Plots FST heatmap from fst calculator module | By Eduardo Charvel  and  Homère J. Alves Monteiro
# Cleans the environment ~ 
rm(list=ls())


#install.package("pacman")
# Loads required packages ~
pacman::p_load(pheatmap, tidyverse, reshape2)

# Loads Fst table ~
# Don't forget to remove the "sample" string in the .txt file
data <- read.table("/Users/sjr729/Desktop/GitHub/Gskimming/00_data/Skmer/Oedulis/RawCov/Fst/oedulis_skmer1_n214_wcfst.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE,)
str(data)

#Remove the redundant first column, but the scalability of this need to be checked
long_df <- data %>%
  pivot_longer(
    cols = -sample,  # Select all columns except 'sample' to pivot
    names_to = "Pop2",  # Column that will contain the original column names
    values_to = "Weighted"  # Column that will contain the values from the pivoted columns
  ) %>%
  rename(Pop1 = sample)  # Rename 'sample' column to 'Pop1'
str(long_df)


# Ensure 'melted_Fst' is loaded as data.frame
melted_Fst <- as.data.frame(long_df)


# Gets Fst label ~
Fst.label = expression(italic("F")[ST])


# Convert population names to factors to preserve order
melted_Fst$Pop1 <- factor(melted_Fst$Pop1, levels = unique(melted_Fst$Pop1))
melted_Fst$Pop2 <- factor(melted_Fst$Pop2, levels = unique(melted_Fst$Pop2))

# Filter to get only the upper triangular part
upper_tri_melted_Fst <- melted_Fst %>%
  filter(as.numeric(Pop1) <= as.numeric(Pop2))

# Create a new column for coloring purposes, keeping the negative values as zero for coloring
upper_tri_melted_Fst$WeightedColor <- pmax(upper_tri_melted_Fst$Weighted, 0)

# Create a new column for identifying categories
upper_tri_melted_Fst$Category <- ifelse(upper_tri_melted_Fst$Weighted > 1, "greater_than_1", "within_0_to_1")


# Plotting the upper triangular heatmap with custom color breaks, including NA values in grey
Fst_Plot <- ggplot() +
  geom_tile(data = upper_tri_melted_Fst, 
            aes(x = Pop1, y = Pop2, fill = Weighted), color = "#ffffff", lwd = 1.5, linetype = 1, na.rm = FALSE) +
  geom_text(data = upper_tri_melted_Fst, 
            aes(x = Pop1, y = Pop2, label = ifelse(is.na(Weighted), "NA", round(Weighted, digits = 2))), 
            color = "white", size = 4) +
  coord_fixed() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), position = "right") +
  scale_fill_gradientn(colors = c("black", "blue", "red"), 
                       values = scales::rescale(c(0, 0.1, 1)),
                       limits = c(0, 1),
                       breaks = c(0, 0.1, 1),
                       labels = c("0", "0.1", "1"),
                       na.value = "black",  # Set NA values to grey
                       name = "Weighted") +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 4)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0.005, b = 0.005, r = .2, l = .2, unit = "cm"),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 18, angle = 60, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 18),
        axis.ticks = element_line(color = "#000000", size = .3),
        legend.position = "right",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.title = element_text(colour = "#000000", size = 18, face = "bold"),
        legend.text = element_text(colour = "#000000", size = 10))


# Save the plot
ggsave(Fst_Plot, file = "~/Desktop/GitHub/Gskimming/02_figures/Oedulis/FST/oedulis_skmer1_n214_locality.png",
       scale = 1, width = 20, height = 15, dpi = 300)



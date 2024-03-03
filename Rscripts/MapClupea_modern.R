### The BEGINNING ~~~~~
##
# ~ Plots --MAP by Homère J. Alves Monteiro
rm(list = ls(all = TRUE))
# Loads required packages ~
pacman::p_load(devtools, tidyverse, ggrepel,geosphere, knitr, sf, rgdal, marmap, sdmpredictors, ggcorrplot, raster, gdistance, ade4, cowplot)
pacman::p_load(devtools, rnaturalearth, rnaturalearthdata, rgeos, sf, ggspatial,
               tidyverse, ggforce, ggstar, ggrepel, extrafont, cowplot)

setwd(dir = "~/Desktop/GitHub/Gskimming/")
# First download shape data in your working directory
#download.file(“http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip”, “countries.zip”)
# Then unzip
#unzip(“countries.zip”)
pop_info <- read_csv("01_infofiles/ClupeaAtmore/ClupeaModern_annot.csv")
set.seed(1)
str(pop_info)

# Imports .shp files ~
Global <- ne_countries(scale = "large", returnclass = "sf")
MR <- st_read("01_infofiles/ClupeaAtmore/Map_shapefile/ne_10m_admin_0_countries.shp")
MR$name <- as.factor(MR$name)

# Loads coordinates ~
Coords_Global<-pop_info
Coords_Global$Longitude <- as.numeric(pop_info$longitude)
Coords_Global$Latitude <- as.numeric(pop_info$latitude)

# Transforms coordinates ~
Coords_Global_sf <- st_as_sf(Coords_Global, coords = c("Longitude", "Latitude"), crs = 4326)



Clupea_Map <-
  ggplot() +
  geom_sf(data = MR, fill = "#f7fcfd", colour = NA) +
  geom_sf(data = Global, fill = "#fff7f3", colour = "#000000") +
  geom_sf(data = Coords_Global_sf, aes(fill = population, shape = sample_description), colour = "#000000", size = 6) +
  coord_sf(xlim = c(-11, 31), ylim = c(67, 49.75), expand = FALSE)+
  scale_shape_manual(values = c(21, 22, 24, 25, 23)) +
  scale_x_continuous(breaks = seq(-10, 32, by = 10)) +
  scale_y_continuous(breaks = seq(-20, 70, by = 10)) +
  annotation_north_arrow(location = "br", which_north = "false", style = north_arrow_fancy_orienteering,
                         height = unit(2.5, "cm"), width = unit(3, "cm"),
                         pad_x = unit(13.5, "in"), pad_y = unit(.2, "in")) +
  annotation_scale(location = "br", line_width = 2, text_cex = 1.35, style = "ticks", pad_y = unit(.3, "in")) +
  theme(panel.background = element_rect(fill = "#bfd3e6"),
        panel.border = element_rect(colour = "black", size = .5, fill = NA),
        panel.grid.major = element_line(color = "#000000", linetype = "dashed", size = .00005),
        plot.margin = margin(t = .005, b = .005, r = .2, l = .2, unit = "cm"),
        legend.background = element_rect(fill = "#fff7f3", size = 0.15, color = "#5e5e5e", linetype = "dotted"),
        legend.key = element_blank(),
        legend.position = c(.143, .825),
        axis.title = element_blank(),
        axis.text = element_text(colour = "#000000", size = 15, face = "bold"),
        axis.ticks = element_line(color ="black", size = .5)) +
  guides(colour = guide_legend(title = "Location", title.theme = element_text(size = 14, face = "bold"),
                               label.theme = element_text(size = 12),
                               override.aes = list(shape = 21, size = 4)),
         fill = guide_legend(title = "Location", title.theme = element_text(size = 14, face = "bold"),
                             label.theme = element_text(size = 12),
                             override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(title = "Species", title.theme = element_text(size = 14, face = "bold"),
                              label.theme = element_text(size = 12), override.aes = list(size = 4)))
Clupea_Map


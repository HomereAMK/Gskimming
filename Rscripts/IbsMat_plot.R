# Cleans the environment ~ 
rm(list=ls())

# Sets working directory ~
setwd("~/Desktop/Scripts/")
#install.packages("RcppCNPy")
# Loads required packages ~
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)

### The BEGINNING ~~~~~
##
# ~ Plot a PCoa/MDS with jc distance matrix from skmer -distance | from Homère J. Alves Monteiro and Eduardo Charvel
# ~ also Plot a PCoa/MDS with .ibs matrix from angsd -doIBS | from Homère J. Alves Monteiro and Eduardo Charvel

# R 4.3.2 GUI 1.80 Big Sur ARM build (8281)

### Clearing environment and setting working directory
rm(list = ls(all = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Loading necessary packages
library(tidyverse)
library(dplyr)
library(ape)
library(cowplot)
library(knitr)
library(scales)
library(miscTools)
if (!requireNamespace("dbscan", quietly = TRUE)) install.packages("dbscan")
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
library(dbscan)
library(mclust)

### Sourcing required functions
source("individual_mdsGskim_functions_hjam.R")

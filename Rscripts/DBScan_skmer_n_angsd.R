### The BEGINNING ~~~~~
##
# ~ Eval of Clustering performances for MDS with distance matrix from ANGSD (-doIBS) and Skmer distance | from Homère J. Alves Monteiro and Eduardo Charvel
### Clearing environment and setting working directory
rm(list = ls(all = TRUE))

# R 4.3.2 GUI 1.80 Big Sur ARM build (8281)
library(mclust)
library(dbscan)

#### Clustering evaluation ####
##Skmer
pcoa_table_skmer <- read.csv("~/Desktop/GitHub/Gskimming/00_data/Angsd/Clupea/RawCov/MDS/Angsd_rawcov_AtmoreClupea_pcoa_table-mat_Mar24.csv")
#  We'll focus on dist_1 and dist_2 for clustering
data <- pcoa_table_skmer[, c("dist_1", "dist_2")]
ggplot(data, aes(x = dist_1, y = dist_2)) + geom_point()
dbscan::kNNdistplot(data, k =  2)

# Step 3: Apply the DBSCAN clustering algorithm
# adjust 'eps' and 'minPts'
dbscan_result <- dbscan(data, eps = 0.5, minPts = 2)

# Step 4: Determine the ground truth based on unique populations
ground_truth <- as.factor(pcoa_table_skmer$population)

# Step 5: Calculate the Adjusted Rand Index (ARI)
# Note: You need to convert clustering labels to a factor for ARI calculation
cluster_labels <- factor(dbscan_result$cluster)
ari_result <- adjustedRandIndex(cluster_labels, ground_truth)

# Print the ARI result
print(ari_result)



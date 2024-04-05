### The BEGINNING ~~~~~
##
# ~ Eval of Clustering performances for MDS with distance matrix from ANGSD (-doIBS) and Skmer distance | from Homère J. Alves Monteiro 
### Clearing environment and setting working directory
rm(list = ls(all = TRUE))

# R 4.3.2 GUI 1.80 Big Sur ARM build (8281)
library(mclust)
library(dbscan)
install.packages("chngpt")
library("chngpt")
#### Clustering evaluation ####
##Angsd
pcoa_table_angsd <- read.csv("~/Desktop/GitHub/Gskimming/00_data/Angsd/Clupea/RawCov/MDS/Angsd_rawcov_AtmoreClupea_pcoa_table-mat_Mar24.csv")
#  We'll focus on dist_1 and dist_2 for clustering
data <- pcoa_table_angsd[, c("dist_1", "dist_2")]
ggplot(data, aes(x = dist_1, y = dist_2)) + geom_point()
y<-sort(dbscan::kNNdist(data, k =  2))
dbscan::kNNdistplot(data, k =  2)
chngptm.xy(x=1:(length(y)),y=y,type="segmented")
print(y[chngptm.xy(x=1:(length(y)),y=y,type="segmented")$chngpt])

# Step 3: Apply the DBSCAN clustering algorithm
# adjust 'eps' and 'minPts'
dbscan_result <- dbscan(data, eps = y[chngptm.xy(x=1:(length(y)),y=y,type="segmented")$chngpt], 
                        minPts = 2*2) # rule of thumbs is minPts = 2 * dim, eps = y-axis value where curve make an elbow

# Step 4: Determine the ground truth based on unique populations
ground_truth <- as.factor(pcoa_table_angsd$population)

# Step 5: Calculate the Adjusted Rand Index (ARI)
# Note: You need to convert clustering labels to a factor for ARI calculation
cluster_labels <- factor(dbscan_result$cluster)
ari_result <- adjustedRandIndex(cluster_labels, ground_truth)
print(cluster_labels)
# Print the ARI result
print(ari_result)


##Skmer
pcoa_table_skmer <- read.csv("~/Desktop/GitHub/Gskimming/00_data/Skmer/Clupea/RawCov/MDS/pcoa_table-mat_40pc_4x_jc-24.02.24_Clupea.csv")
#  We'll focus on dist_1 and dist_2 for clustering
data <- pcoa_table_skmer[, c("dist_1", "dist_2")]
ggplot(data, aes(x = dist_1, y = dist_2)) + geom_point()
dbscan::kNNdistplot(data, k =  2)
y<-sort(dbscan::kNNdist(data, k =  2))
chngptm.xy(x=1:(length(y)),y=y,type="segmented")
print(y[chngptm.xy(x=1:(length(y)),y=y,type="segmented")$chngpt])

# Step 3: Apply the DBSCAN clustering algorithm
# adjust 'eps' and 'minPts'
dbscan_result <- dbscan(data, eps = y[chngptm.xy(x=1:(length(y)),y=y,type="segmented")$chngpt], minPts = 2*2) # rule of thumbs is minPts = 2 * dim, eps = y-axis value where curve make an elbow

# Step 4: Determine the ground truth based on unique populations
ground_truth <- as.factor(pcoa_table_angsd$population)

# Step 5: Calculate the Adjusted Rand Index (ARI)
# Note: You need to convert clustering labels to a factor for ARI calculation
cluster_labels <- factor(dbscan_result$cluster)
ari_result <- adjustedRandIndex(cluster_labels, ground_truth)
print(cluster_labels)
# Print the ARI result
print(ari_result)

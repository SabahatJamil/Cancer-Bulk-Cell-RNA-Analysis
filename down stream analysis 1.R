install.packages("openxlsx")
library(openxlsx)
#counts <- read.csv(.csv"
                   #, row.names = 1)


output_file_path <- "D:/Docs/Ayaas Lab Work/Assignment 18 bulk cell Analysis/Batch Analysis/"

input_file_path <- "D:/Docs/Ayaas Lab Work/Assignment 19/comb.xlsx"

#######################################
##HERE THE FILE IS LOADED AND THE PREPROCESSING IS DONE
#####################################

counts <- read.xlsx(input_file_path
                    ,sheet="Sheet1")
counts <- na.omit(counts)
counts <- counts[!duplicated(counts$Gene), ]
row.names(counts) <- counts$Gene
df <- subset(counts, select = -Gene)
counts <- subset(counts, select = -Gene)

##################################################################################
###HERE ITS FUTHER PROCESSED, NORMALIZED AND PCA ETC ARE MADE
#################################################################################

col_data <- data.frame(row.names = colnames(counts), condition = rep("condition", ncol(counts)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~1)
vsd <- vst(dds, blind = TRUE)
normalized_counts <- assay(vsd)

pca_result <- PCA(t(normalized_counts), graph = FALSE)

fviz_pca_ind(pca_result,
             geom.ind = "point",
             pointsize = 3,
             title = "PCA of Bulk RNA-seq Counts")

######################
k <- 5
km_res <- kmeans(pca_result$x, centers = k, nstart = 25)

# Visualize the PCA individuals with k-means clusters
fviz_cluster(km_res, 
             data = pca_result$x,
             geom = "point", 
             pointsize = 3, 
             ellipse = TRUE, 
             main = "PCA with K-means Clusters")

###############################################
## HERE THE K-MEANS CLUSTERING IS DONE TO MAKE CLUSTERS AND PLOTT THEM
#################################################
set.seed(123)  
kmeans_result <- kmeans(pca_result$x$coord[, 1:2], centers = 5)
custom_palette <- c("#00AFBB", "#E7B800","#68af2c","#efc39d","#f82927")
#, "#FC4E07", "#00BFC4", 
pca_plot <- fviz_pca_ind(pca_result,
             geom.ind = "point",
             pointsize = 3,
             title = "PCA with Clusters",
             col.ind = as.factor(kmeans_result$cluster),
             palette = custom_palette)  

pca_plot + geom_text(aes(label = ifelse(colnames(counts) %in% first_19_colnames, 
                                        colnames(counts), 
                                        "")), 
                     vjust = -1, 
                     hjust = 1.5, 
                     size = 3, 
                     color = "black")

#############################################################
#sequecing depth plotting (IF NEEDED)
############################################################

sequencing_depth <- colSums(counts)
print(sequencing_depth)

depth_df <- data.frame(Sample = names(sequencing_depth), Depth = sequencing_depth)

ggplot(depth_df, aes(x = Sample, y = Depth)) +
  geom_bar(stat = "identity", fill = "#00AFBB") +
  labs(title = "Sequencing Depth per Sample", x = "Sample", y = "Sequencing Depth") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##############################################################
#OUTLIER DETECTION CODE
#############################################################
pca_coords <- pca_result$ind$coord[, 1:2]  # Using first two principal components


center <- colMeans(pca_coords)


distances <- apply(pca_coords, 1, function(row) sqrt(sum((row - center)^2)))


# Here, we use 1.5 times the interquartile range (IQR) above the third quartile as the threshold
threshold <- quantile(distances, 0.75) + 1.5 * IQR(distances)


outliers <- which(distances > threshold)
outlier_samples <- rownames(pca_coords)[outliers]

cat("Outliers detected:", outlier_samples, "\n")  
 
df <- counts
###################################################################
##INTER-CLUSTER DEGS ARE CALCULTED HERE (MAKES ALL POSSIBLE COMBINATIONS) (TAKES SOME TIME)
####################################################################
col_data$cluster <- as.factor(kmeans_result$cluster)

clusters <- levels(col_data$cluster)
deg_results <- list()  

for (i in 1:(length(clusters) - 1)) {
  for (j in (i + 1):length(clusters)) {
    cluster1 <- clusters[i]
    cluster2 <- clusters[j]
    comparison_name <- paste("Cluster", cluster1, "vs", "Cluster", cluster2)
    subset_col_data <- col_data[col_data$cluster %in% c(cluster1, cluster2), ]
    subset_df <- df[, rownames(subset_col_data)]
    dds_subset <- DESeqDataSetFromMatrix(countData = subset_df,
                                         colData = subset_col_data,
                                         design = ~ cluster)
    dds_subset$cluster <- relevel(dds_subset$cluster, ref = cluster1) 
    
    dds_subset <- DESeq(dds_subset)
    res <- results(dds_subset)
    
    deg_results[[comparison_name]] <- res
  }
}

for (comparison_name in names(deg_results)) {
  cat("\n", comparison_name, ":\n")
  print(summary(deg_results[[comparison_name]]))
}

sig_deg_results <- lapply(deg_results, function(res) {
  res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
})

View(res)

################################################################################
##SAVES THE CLUSTER SAMPLE NAMES (MEANING WHICH SAMPLES BELONG TO WHICH CLUSTER)
################################################################################


# Create a new Excel workbook
wb <- createWorkbook()

# Loop through each cluster and add its sample names to the workbook
for (cluster in unique(col_data$cluster)) {
  # Get sample names in the current cluster
  samples_in_cluster <- rownames(col_data[col_data$cluster == cluster, ])
  
  # Create a data frame with sample names for easy writing to Excel
  samples_df <- data.frame(Sample = samples_in_cluster)
  
  # Add a new sheet with the cluster name
  addWorksheet(wb, sheetName = paste("Cluster", cluster))
  
  # Write the data frame to the sheet
  writeData(wb, sheet = paste("Cluster", cluster), x = samples_df)
}

# Save the workbook to an Excel file
saveWorkbook(wb, paste0(output_file_path,"Cluster_Samples.xlsx"), overwrite = TRUE)


###############################################################################
##SAVES THE ORIGINAL INTER-CLUSTER DEGS FILE
###############################################################################
wb <- createWorkbook()

# Loop through each DEG result and add it as a sheet to the workbook
for (comparison_name in names(deg_results)) {
  # Get the DEG result for the current comparison
  res <- deg_results[[comparison_name]]
  
  # Filter for significant DEGs (if needed)
  sig_res <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
  
  # Convert to data frame and add gene names as a column
  sig_res_df <- as.data.frame(sig_res)
  sig_res_df$Gene <- rownames(sig_res_df)  # Add gene names as a new column
  
  # Move 'Gene' column to the first position
  sig_res_df <- sig_res_df[, c("Gene", setdiff(names(sig_res_df), "Gene"))]
  
  # Add a new sheet with the comparison name
  addWorksheet(wb, sheetName = comparison_name)
  
  # Write the data frame to the sheet
  writeData(wb, sheet = comparison_name, x = sig_res_df)
}

# Save the workbook to an Excel file
saveWorkbook(wb, paste0(output_file_path,"DEGs.xlsx"), overwrite = TRUE)


#####################################################################################
#HIERARCHICAL CLUSTERING (OPTIONAL NEEDED, IF ASKED) ITS A SORT OF CONFIRMATION
######################################################################################
pca_result <- PCA(t(normalized_counts), graph = FALSE)

# Get PCA coordinates for the first few principal components
pca_coords <- pca_result$ind$coord[, 1:3]  # Adjust number of PCs as needed

# Calculate the distance matrix
dist_matrix <- dist(pca_coords)

# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plot the dendrogram
plot(hclust_result, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "")


cluster_assignments <- cutree(hclust_result, k = 5)  # Adjust k as needed

# Add cluster assignments to your data
pca_coords <- as.data.frame(pca_coords)
pca_coords$Cluster <- as.factor(cluster_assignments)

custom_palette <- c("#00AFBB", "#E7B800","#68af2c","#efc39d","#f82927","#c9ffc7")
#, 
# Visualize the clusters on PCA
fviz_pca_ind(pca_result,
             geom.ind = "point",
             pointsize = 3,
             col.ind = pca_coords$Cluster,
             palette = custom_palette,
             title = "PCA with Hierarchical Clusters")

######################################################################
#PLOTS THE OUTLIERS ON K-MEANS CLUSTERING PLOT
#####################################################################
pca_result <- PCA(t(normalized_counts), graph = FALSE)

# Get PCA coordinates for the first two components
pca_coords <- as.data.frame(pca_result$ind$coord[, 1:2])
colnames(pca_coords) <- c("PC1", "PC2")
pca_coords$Sample <- rownames(pca_coords)

# K-means clustering
set.seed(123)
kmeans_result <- kmeans(pca_coords[, 1:2], centers = 5)
pca_coords$Cluster <- as.factor(kmeans_result$cluster)

# Calculate center and distances for outlier detection
center <- colMeans(pca_coords[, 1:2])
distances <- apply(pca_coords[, 1:2], 1, function(row) sqrt(sum((row - center)^2)))
threshold <- quantile(distances, 0.75) + 1.5 * IQR(distances)
pca_coords$Outlier <- distances > threshold

# Define custom colors
custom_palette <- c("#00AFBB", "#E7B800","#68af2c","#efc39d","#f82927")

# Plot PCA with ggplot2, adding sample names and highlighting outliers
ggplot(pca_coords, aes(x = PC1, y = PC2, label = Sample, color = Cluster)) +
  geom_point(size = 3, aes(shape = Outlier)) +
  scale_color_manual(values = custom_palette) +
  geom_text(aes(label = ifelse(Outlier, Sample, "")), hjust = 1.1, vjust = 1.1, color = "red") +
  ggtitle("PCA with Clusters and Outliers") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(shape = "Outlier") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) 

##############################
#for adding the sample titles on the samples in the plot

p <- fviz_pca_ind(pca_result,
                  geom.ind = "point",
                  pointsize = 3,
                  title = "PCA with Clusters",
                  col.ind = as.factor(kmeans_result$cluster),
                  palette = custom_palette)

p + geom_text(aes(label = rownames(pca_result$ind$coord)), vjust = -1, size = 3)


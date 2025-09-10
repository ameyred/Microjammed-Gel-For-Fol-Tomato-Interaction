# This script was written on 21st May, 2025
# For making PCA plot after log2 transformation on counts matrix
# Comparing Axenic Fol vs Gamborg Broth vs Gamborg Microgel

library(ggplot2)

#Load count matrix
count_matrix <- read.csv("~/Amey_Lab/Microgel_Project/All_Combined_Counts.csv", row.names = 1)

# Convert to numeric matrix 
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "numeric"

# Remove genes with all zero expression across all samples
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

# Log2-transform with pseudocount to avoid log(0)
log_mat <- log2(count_matrix + 1)

#Perform PCA (transpose so samples are rows)
pca <- prcomp(t(log_mat), center = TRUE, scale. = TRUE)

# Extract treatment group from sample names
samples <- colnames(count_matrix)
treatment <- factor(sub("_.*", "", samples))

# Step 7: Create PCA plotting dataframe
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$treatment <- treatment
pca_df$sample <- samples

# Step 8: Plot PCA 
ggplot(pca_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 10) +  # Increased dot size
  scale_color_manual(values = c(
    "Axenic" = rgb(46, 37, 133, maxColorValue = 255),
    "Broth" = rgb(93, 168, 153, maxColorValue = 255),
    "Gel" = rgb(194, 106, 119, maxColorValue = 255)
  )) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(linewidth = 1.5, color = "black"),  # Use linewidth instead of size
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) + 
  xlab(paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  ggtitle("PCA Plot (log2-transformed expression values)\nAxenic vs Broth vs Microgel")










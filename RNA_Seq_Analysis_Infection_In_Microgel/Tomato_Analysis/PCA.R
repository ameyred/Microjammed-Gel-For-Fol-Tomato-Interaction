# This script was updated on 13th June, 2025

library(ggplot2)

# Load and prepare count matrix
count_matrix <- read.csv("~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/All_Combined_Counts.csv", row.names = 1)
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "numeric"
count_matrix <- count_matrix[rowSums(count_matrix > 0) >= 3, ]

# Log transform and run PCA
log_mat <- log2(count_matrix + 1)
pca <- prcomp(t(log_mat), center = TRUE, scale. = TRUE)

# Extract sample names
samples <- colnames(count_matrix)

# Parse condition from sample names
condition <- factor(sub("_[0-9]+$", "", samples), levels = c(
  "Mock_Tip", "Infected_Tip", "Mock_Root", "Infected_Root"
))

# Create PCA dataframe
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$sample <- samples
pca_df$condition <- condition

# Define fill color mapping using your custom RGB values
fill_colors <- c(
  "Mock_Tip" = rgb(46, 37, 133, maxColorValue = 255),
  "Infected_Tip" = rgb(51, 117, 56, maxColorValue = 255),
  "Mock_Root" = rgb(220, 205, 125, maxColorValue = 255),
  "Infected_Root" = rgb(126, 41, 84, maxColorValue = 255)
)

# Final PCA plot (all circles, filled)
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = condition), shape = 21, size = 7, stroke = 1.2, color = "black") +
  scale_fill_manual(values = fill_colors, name = "Condition") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1.5, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  xlab(paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  ggtitle("PCA Plot (log2-transformed expression values)\nMock and Fol-infected root tips and bulk")


# Load necessary libraries
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# Read the data from CSV
data <- read.csv("~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Analysis/Tip_VS_Root/Up_Genes_Heat_Map.csv")

# Handle duplicate row names (if any)
data <- data %>% distinct(data[, 1], .keep_all = TRUE)

# Set row names using the first column and remove it from data
rownames(data) <- data[, 1]
data <- data[, -1]

# Extract relevant columns for the heatmap
heatmap_data <- data[, c("Log2FC_Tip", "Log2FC_Root")]

# Remove rows with NA values
heatmap_data <- na.omit(heatmap_data)

# Define breaks for better color contrast
breaks <- breaks <- seq(-0.5, 10.5, by = 1.5)

pheatmap(heatmap_data, 
         cluster_cols = TRUE, 
         cluster_rows = TRUE, 
         treeheight_row = 0, 
         treeheight_col = 0, 
         show_rownames = TRUE, 
         show_colnames = TRUE,  
         border_color = "black", 
         cellwidth = 25, 
         cellheight = 10,
         fontsize = 8,               # Overall font size
         fontsize_row = 6,           # Row label font size
         fontsize_col = 8,           # Column label font size
         breaks = breaks, 
         color = brewer.pal(n = 5, name = "Greens"))

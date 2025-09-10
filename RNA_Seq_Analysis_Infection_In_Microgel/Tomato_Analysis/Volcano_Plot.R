library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)

# Set base paths
base_dir <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis"
image_dir <- file.path(base_dir, "Figures")

# List of comparisons
comparisons <- c("Mock_Root_VS_Infected_Root", "Mock_Tip_VS_Infected_Tip")

# Loop through each comparison
for (comp in comparisons) {
  file_path <- file.path(base_dir, "Analysis", comp, "DESeq_Output.csv")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  
  # Read data
  data <- read.csv(file_path)
  
  # Count upregulated and downregulated genes
  up_count <- sum(data$log2FoldChange > 2 & data$padj < 0.05, na.rm = TRUE)
  down_count <- sum(data$log2FoldChange < -2 & data$padj < 0.05, na.rm = TRUE)
  
  # Volcano plot with annotations
  p <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(
      size = 1,
      shape = 1,
      aes(color = factor(ifelse(padj < 0.05,
                                ifelse(log2FoldChange < 2,
                                       ifelse(log2FoldChange > -2, "black", "red"),
                                       "red"),
                                "black")))
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", colour = "black") +
    scale_y_continuous(limits = c(0, 300)) +
    scale_x_continuous(limits = c(-15, 15), breaks = seq(-15, 15, 3)) +
    scale_color_manual(values = c("black", "red"),
                       labels = c("Not significant", "Significant")) +
    annotate("text", x = -15, y = 290, label = paste("Downregulated:", down_count),
             hjust = 0, size = 8, color = "red", fontface = "bold") +
    annotate("text", x = 5, y = 290, label = paste("Upregulated:", up_count),
             hjust = 0, size = 8, color = "red", fontface = "bold") +
    theme(
      panel.background = element_rect(fill = 'transparent'),
      axis.line = element_line(size = 1.5, colour = "black"),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14)
    ) +
    ggtitle(paste("Volcano Plot -", comp))
  
  # Save the plot
  output_path <- file.path(image_dir, paste0(comp, ".png"))
  ggsave(output_path, plot = p, width = 10, height = 10, dpi = 300)
  
  message("Saved: ", output_path)
}


library(ggplot2)
library(tidyr)
library(dplyr)

# Read CSV
data <- read.csv("~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Analysis/Tip_VS_Root/Up_Genes_Heat_Map.csv")

# Rename first column if needed
colnames(data)[1] <- "Genes"

# Reshape
data_long <- data %>%
  pivot_longer(cols = c(Log2FC_Tip, Log2FC_Root),
               names_to = "Tissue",
               values_to = "Log2FC")

# Rank genes by Tip log2FC
ranked_genes <- data %>%
  arrange(desc(Log2FC_Tip)) %>%
  pull(Genes)

data_long$Genes <- factor(data_long$Genes, levels = ranked_genes)

#Plot
pdf("~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Figures/Ranked_Bar_Plot_Up_Genes.pdf", width = 20, height = 20) 

ggplot(data_long, aes(x = Genes, y = Log2FC, fill = Tissue)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Log2FC_Tip" = "#5DA899", "Log2FC_Root" = "#C26A77")) +
  coord_flip() +
  scale_y_continuous(limits = c(-2, 12), breaks = seq(-2, 12, by = 2)) +
  labs(x = NULL, y = "Log₂ Fold Change", fill = "Tissue") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    panel.grid.major.y = element_blank(),                              # No horizontal gridlines
    panel.grid.major.x = element_line(color = "gray80", size = 0.5),   # Vertical gridlines
    panel.grid.minor = element_blank(),                                # No minor gridlines
    
    axis.line.x = element_line(color = "black", size = 0.8),           # X-axis line
    axis.line.y = element_line(color = "black", size = 0.8),           # Y-axis line
    
    axis.text.x = element_text(size = 16),   # Bigger x-axis numbers
    axis.text.y = element_blank(),           # Remove gene names
    axis.ticks.y = element_blank(),
    
    axis.title.y = element_text(size = 16, face = "bold"),
    
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "top"
  )

dev.off()
# This script was written on 21st May, 2025
# For making Volcano plot on DESEq Output on microgel project

#Volcano plot Fol Axenic vs Broth 
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel)
library(ggplot2)

data<-read.csv("~/Amey_Lab/Microgel_Project/Axenic_VS_Broth/DESEeq_Output.csv")

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
  geom_vline(xintercept = -2, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 2, linetype = "dashed", colour = "black") +
  scale_y_continuous(limits = c(0, 300)) +
  scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
  scale_color_manual(values = c("black", "red"), 
                     labels = c("Not significant", "Significant")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    axis.line = element_line(size = 1.5, colour = "black"),  # Thicker and black axes
    axis.text.x = element_text(size = 18),  # Increase x-axis text size
    axis.text.y = element_text(size = 18),  # Increase y-axis text size
    axis.title.x = element_text(size = 20), # Optional: x-axis title size
    axis.title.y = element_text(size = 20), # Optional: y-axis title size
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)   # Optional: bigger legend text
  )

p
ggsave("~/Amey_Lab/Microgel_Project/Axenic_VS_Broth/Volcano_Plot.png", plot = p, width = 10, height = 10, dpi = 300)


#---------------------------------------------------------------------------------------------------------------



#Volcano plot Fol Axenic vs Gel 
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel)
library(ggplot2)

data<-read.csv("~/Amey_Lab/Microgel_Project/Axenic_VS_Gel/DESEeq_Output.csv")
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
  geom_vline(xintercept = -2, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 2, linetype = "dashed", colour = "black") +
  scale_y_continuous(limits = c(0, 300)) +
  scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
  scale_color_manual(values = c("black", "red"), 
                     labels = c("Not significant", "Significant")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    axis.line = element_line(size = 1.5, colour = "black"),  # Thicker and black axes
    axis.text.x = element_text(size = 18),  # Increase x-axis text size
    axis.text.y = element_text(size = 18),  # Increase y-axis text size
    axis.title.x = element_text(size = 20), # Optional: x-axis title size
    axis.title.y = element_text(size = 20), # Optional: y-axis title size
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)   # Optional: bigger legend text
  )

p
ggsave("~/Amey_Lab/Microgel_Project/Axenic_VS_Gel/Volcano_Plot.png", plot = p, width = 10, height = 10, dpi = 300)


#---------------------------------------------------------------------------------------------------------------



#Volcano plot Fol Broth vs Gel 
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel)
library(ggplot2)

data<-read.csv("~/Amey_Lab/Microgel_Project/Broth_VS_Gel/DESEeq_Output.csv")
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
  geom_vline(xintercept = -2, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 2, linetype = "dashed", colour = "black") +
  scale_y_continuous(limits = c(0, 300)) +
  scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
  scale_color_manual(values = c("black", "red"), 
                     labels = c("Not significant", "Significant")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    axis.line = element_line(size = 1.5, colour = "black"),  # Thicker and black axes
    axis.text.x = element_text(size = 18),  # Increase x-axis text size
    axis.text.y = element_text(size = 18),  # Increase y-axis text size
    axis.title.x = element_text(size = 20), # Optional: x-axis title size
    axis.title.y = element_text(size = 20), # Optional: y-axis title size
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)   # Optional: bigger legend text
  )

p
ggsave("~/Amey_Lab/Microgel_Project/Broth_VS_Gel/Volcano_Plot.png", plot = p, width = 10, height = 10, dpi = 300)







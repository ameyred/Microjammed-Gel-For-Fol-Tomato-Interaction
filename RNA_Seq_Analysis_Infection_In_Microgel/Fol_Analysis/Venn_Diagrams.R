#Written on 16th June for big experiment (Fol side) venn diagram between root tip and bulk root 


library(grid)

# Start a new plot
png("~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Analysis/Figures/Up_Fol_Tip_Fol_Root.png", width = 6, height = 4, units = "in", res = 300)

grid.newpage()

# Define circle style
circle_style <- gpar(fill = "#2e258580", col = "#2e2585", lwd = 2)
circle_style2 <- gpar(fill = "#5da89980", col = "#5da899", lwd = 2)

# Draw circles manually
grid.circle(x = 0.44, y = 0.5, r = 0.2, gp = circle_style)   # Circle A
grid.circle(x = 0.56, y = 0.5, r = 0.2, gp = circle_style2)  # Circle B

# Add gene counts (update x/y to position as needed)
text_style <- gpar(fontsize = 14)
grid.text("223", x = 0.37, y = 0.50, gp = text_style)  # Unique to A
grid.text("237",  x = 0.63, y = 0.50, gp = text_style)  # Unique to B
grid.text("351", x = 0.50, y = 0.50, gp = text_style)  # A & B

dev.off()

#--------------------------------------------------------------------


# Downregulated Genes

# Start a new plot
png("~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Analysis/Figures/Down_Fol_Tip_Fol_Root.png", width = 6, height = 4, units = "in", res = 300)

grid.newpage()

# Define circle style
circle_style <- gpar(fill = "#2e258580", col = "#2e2585", lwd = 2)
circle_style2 <- gpar(fill = "#5da89980", col = "#5da899", lwd = 2)

# Draw circles with more overlap
grid.circle(x = 0.44, y = 0.5, r = 0.2, gp = circle_style)   # Circle A
grid.circle(x = 0.56, y = 0.5, r = 0.2, gp = circle_style2)  # Circle B

# Add gene counts
text_style <- gpar(fontsize = 14)
grid.text("259", x = 0.37, y = 0.50, gp = text_style)  # Unique to A
grid.text("367", x = 0.63, y = 0.50, gp = text_style)  # Unique to B
grid.text("726", x = 0.50, y = 0.50, gp = text_style)  # Overlap

# Save the image
dev.off()
#Written on 21st June to make a Venn overlap to show upregulated genes between root tip and root bulk. 


library(grid)

# Start a new plot
png("~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Figures/Venn_Root_Tip_Up.png", width = 6, height = 4, units = "in", res = 300)

grid.newpage()

# Define circle style
circle_style <- gpar(fill = "#2e258580", col = "#2e2585", lwd = 2)
circle_style2 <- gpar(fill = "#5da89980", col = "#5da899", lwd = 2)

# Draw circles manually
grid.circle(x = 0.44, y = 0.5, r = 0.2, gp = circle_style)   # Circle A
grid.circle(x = 0.56, y = 0.5, r = 0.2, gp = circle_style2)  # Circle B

# Add gene counts (update x/y to position as needed)
text_style <- gpar(fontsize = 14)
grid.text("135", x = 0.37, y = 0.50, gp = text_style)  # Unique to A
grid.text("54",  x = 0.63, y = 0.50, gp = text_style)  # Unique to B
grid.text("43", x = 0.50, y = 0.50, gp = text_style)  # A & B

dev.off()


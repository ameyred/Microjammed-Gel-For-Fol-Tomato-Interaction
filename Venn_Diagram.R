library(VennDiagram)
library(grid)


# Upregulated Genes
grid.newpage()

# Venn Layout
draw.triple.venn(
  area1 = 100, area2 = 100, area3 = 100,  # equal size circles
  n12 = 20, n13 = 20, n23 = 20, n123 = 10,  # just enough to form visible overlaps
  category = rep("", 3),
  fill = c("#2e2585", "#5da899", "#c26a77"),
  alpha = 0.4,
  lwd = 2,
  cat.cex = 0,
  cex = 0
)

# Add gene counts 
text_style <- gpar(fontsize = 24)

grid.text("1737", x = 0.23, y = 0.70, gp = text_style)  # Unique to A
grid.text("317",  x = 0.77, y = 0.70, gp = text_style)  # Unique to B
grid.text("200",  x = 0.50, y = 0.20, gp = text_style)  # Unique to C

grid.text("1568", x = 0.50, y = 0.80, gp = text_style)  # A & B
grid.text("0",    x = 0.35, y = 0.45, gp = text_style)  # A & C
grid.text("75",   x = 0.65, y = 0.45, gp = text_style)  # B & C

grid.text("19",   x = 0.50, y = 0.60, gp = text_style)  # A & B & C

#--------------------------------------------------------------------


# Downregulated Genes
grid.newpage()

# Venn Layout
draw.triple.venn(
  area1 = 100, area2 = 100, area3 = 100,  # equal size circles
  n12 = 20, n13 = 20, n23 = 20, n123 = 10,  # just enough to form visible overlaps
  category = rep("", 3),
  fill = c("#2e2585", "#5da899", "#c26a77"),
  alpha = 0.4,
  lwd = 2,
  cat.cex = 0,
  cex = 0
)

# Add gene counts 
text_style <- gpar(fontsize = 24)

grid.text("422", x = 0.23, y = 0.70, gp = text_style)  # Unique to A
grid.text("365",  x = 0.77, y = 0.70, gp = text_style)  # Unique to B
grid.text("1044",  x = 0.50, y = 0.20, gp = text_style)  # Unique to C

grid.text("1054", x = 0.50, y = 0.80, gp = text_style)  # A & B
grid.text("0",    x = 0.35, y = 0.45, gp = text_style)  # A & C
grid.text("206",   x = 0.65, y = 0.45, gp = text_style)  # B & C

grid.text("61",   x = 0.50, y = 0.60, gp = text_style)  # A & B & C




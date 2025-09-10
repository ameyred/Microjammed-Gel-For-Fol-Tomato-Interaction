library(eulerr)

# For upregulated genes
fit <- euler(c(
  "A" = 1737,
  "B" = 317,
  "C" = 200,
  "A&B" = 1568,
  "A&C" = 0,
  "B&C" = 75,
  "A&B&C" = 19
))

plot(fit,
     fills = list(fill = c("#2e2585", "#5da899", "#c26a77"), alpha = 0.6),
     labels = FALSE,
     quantities = list(cex = 2),  # ← Increase text size here
     main = "Euler Diagram of Gene Set Overlaps")


# For downregulated genes
fit <- euler(c(
  "A" = 422,
  "B" = 365,
  "C" = 1044,
  "A&B" = 1054,
  "A&C" = 0,
  "B&C" = 206,
  "A&B&C" = 61
))

plot(fit,
     fills = list(fill = c("#2e2585", "#5da899", "#c26a77"), alpha = 0.6),
     labels = FALSE,
     quantities = list(cex = 2),  # ← Increase text size here too
     main = "Euler Diagram of Gene Set Overlaps")

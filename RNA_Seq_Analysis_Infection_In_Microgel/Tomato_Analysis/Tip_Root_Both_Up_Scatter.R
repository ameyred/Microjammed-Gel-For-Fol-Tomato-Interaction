# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tibble)

# Load the data
data <- read.csv("~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Analysis/Tip_VS_Root/Up_Both.csv", row.names = 1)


# Compute Kendall's correlation coefficient
kendall_correlation <- cor(data$Log2FC.Tip, data$Log2FC.Root, method = "kendall")

# Perform Kendall's correlation test for statistical significance
kendall_test <- cor.test(data$Log2FC.Tip, data$Log2FC.Root, method = "kendall")

# Print the results
cat("Kendall Correlation Coefficient:", kendall_correlation, "\n")
cat("\nKendall Test Results :\n")
print(kendall_test)

# Create a scatter plot
ggplot(data, aes(x = Log2FC.Root, y = Log2FC.Tip)) +
  geom_point(alpha = 0.6) +  # semi-transparent points
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +  # optional trend line
  labs(
    title = "Scatter Plot of Gene Expression: Tip vs Root",
    x = "Log2 Fold Change in Root",
    y = "Log2 Fold Change in Tip"
  ) +
  theme_minimal()



#Kendall's rank correlation tau

#data:  data$Log2FC.Tip and data$Log2FC.Root
#T = 573, p-value = 0.01072
#alternative hypothesis: true tau is not equal to 0
#sample estimates:
     #tau 
#0.269103 



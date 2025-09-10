# Load libraries
library(DESeq2)
library(tidyverse)


# Base directory
base_dir <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Analysis"

# Find all Combined_Counts.csv files
count_files <- list.files(path = base_dir, pattern = "Combined_Counts.csv", recursive = TRUE, full.names = TRUE)

# Loop over each dataset
for (count_file in count_files) {
  cat("Processing:", count_file, "\n")
  folder <- dirname(count_file)
  
  # Load counts file
  counts_data <- read.csv(count_file, header = FALSE)
  colnames(counts_data) <- as.character(counts_data[1, ])
  counts_data <- counts_data[-1, ]
  rownames(counts_data) <- as.character(counts_data[[1]])
  counts_data <- counts_data[-1]
  counts_data_numeric <- apply(counts_data, 2, function(x) as.numeric(as.character(x)))
  counts_data_table <- as.data.frame(counts_data_numeric)
  rownames(counts_data_table) <- rownames(counts_data)
  counts_data_table <- round(counts_data_table)
  
  # Display preview
  print(head(counts_data_table))
  print(colnames(counts_data_table))
  
  # Load sample info
  sample_info_file <- file.path(folder, "Sample_Info.csv")
  if (!file.exists(sample_info_file)) {
    cat("Sample_Info.csv not found in", folder, "- skipping...\n")
    next
  }
  
  colData <- read.csv(sample_info_file, header = FALSE)
  colnames(colData) <- as.character(colData[1, ])
  colData <- colData[-1, ]
  rownames(colData) <- as.character(colData[[1]])
  colData <- colData[-1]
  col_data_table <- as.data.frame(colData)
  
  # Display sample info preview
  print(col_data_table)
  print(rownames(col_data_table))
  
  # Check name alignment
  if (!all(colnames(counts_data_table) %in% rownames(col_data_table))) {
    cat("Mismatch between count column names and sample info row names in", folder, "\n")
    next
  }
  if (!all(colnames(counts_data_table) == rownames(col_data_table))) {
    cat("Reordering sample info to match count data column order\n")
    col_data_table <- col_data_table[colnames(counts_data_table), ]
  }
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts_data_table,
                                colData = col_data_table,
                                design = ~Treatment)
  
  # Filter low counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # Set reference level
  dds$Treatment <- relevel(dds$Treatment, ref = "Mock")
  
  # Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds)
  results_df <- as.data.frame(res)
  
  # Save output
  output_file <- file.path(folder, "DESeq_Output.csv")
  write.csv(results_df, file = output_file, row.names = TRUE)
  
  cat("Saved DESeq2 results to:", output_file, "\n\n")
}

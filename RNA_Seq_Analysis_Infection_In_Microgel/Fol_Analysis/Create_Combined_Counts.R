# Load necessary packages
library(tximport)
library(readr)

# Path to the transcript-to-gene mapping file
tx2gene_path <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Genome/Transcript_To_Gene.csv"
tx2gene <- read.csv(tx2gene_path, stringsAsFactors = FALSE)
colnames(tx2gene) <- c("transcript_id", "gene_id")

# Base directory containing all sample folders
base_dir <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Analysis"

# Find all Quant_SF_Files folders recursively
quant_dirs <- list.dirs(path = base_dir, recursive = TRUE, full.names = TRUE)
quant_dirs <- quant_dirs[grepl("Quant_Files$", quant_dirs)]

# Loop over each folder
for (quant_dir in quant_dirs) {
  cat("Processing folder:", quant_dir, "\n")
  
  # List quant.sf files
  quant_files <- list.files(path = quant_dir, pattern = ".sf$", full.names = TRUE)
  if (length(quant_files) == 0) {
    cat("No .sf files found in", quant_dir, "- skipping...\n")
    next
  }
  
  # Extract sample names from file paths without ".sf"
  sample_names <- sub("\\.sf$", "", basename(quant_files))
  
  # Perform tximport
  txi <- tximport(files = quant_files, type = "salmon", tx2gene = tx2gene)
  counts <- txi$counts
  combined_counts <- as.data.frame(counts)
  combined_counts$Gene_ID <- rownames(combined_counts)
  combined_counts <- combined_counts[, c("Gene_ID", colnames(combined_counts)[-ncol(combined_counts)])]
  colnames(combined_counts)[-1] <- sample_names
  
  # Remove the header name for Gene_ID column
  colnames(combined_counts)[1] <- ""
  
  # Save to CSV in the parent directory of Quant_SF_Files
  output_dir <- dirname(quant_dir)
  output_file <- file.path(output_dir, "Combined_Counts.csv")
  write.csv(combined_counts, file = output_file, row.names = FALSE)
  cat("Saved combined counts to:", output_file, "\n")
  
  # Create Sample Info
  # Create Sample Info
  sample_info <- data.frame(
    Sample = sample_names,
    Treatment = sub("(_[0-9]+dpi)?_[0-9]+$", "", sample_names),
    stringsAsFactors = FALSE
  )
  
  sample_info_file <- file.path(output_dir, "Sample_Info.csv")
  write.csv(sample_info, file = sample_info_file, row.names = FALSE)
  cat("Saved sample info to:", sample_info_file, "\n\n")
}


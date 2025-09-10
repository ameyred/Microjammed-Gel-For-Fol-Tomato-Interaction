#Written on 12th June to plot PCA for Microgel paper

# Load necessary packages
library(tximport)
library(readr)

# Path to the transcript-to-gene mapping file
tx2gene_path <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Analysis/Fol_Genome/Transcript_To_Gene.csv"
tx2gene <- read.csv(tx2gene_path, stringsAsFactors = FALSE)

# Ensure the column names are correct
colnames(tx2gene) <- c("transcript_id", "gene_id")

# Print the first few rows of the transcript-to-gene mapping to check
print("Transcript-to-Gene Mapping:")
print(head(tx2gene))

# Set the directory where the quant.sf files are located
quant_dir <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Analysis/All_Quant_Files"

# List all quant.sf files
quant_files <- list.files(path = quant_dir, pattern = ".sf$", full.names = TRUE)

# Extract sample names from quant file paths
sample_names <- sub(".*/|.sf$", "", quant_files)

# Perform tximport
txi <- tximport(files = quant_files, type = "salmon", tx2gene = tx2gene)

# Extract counts matrix
counts <- txi$counts

# Convert counts to a data frame
combined_counts <- as.data.frame(counts)

# Assign gene IDs as a separate column
combined_counts$Gene_ID <- rownames(combined_counts)

# Reorder columns to have Gene_ID first
combined_counts <- combined_counts[, c("Gene_ID", colnames(combined_counts)[-ncol(combined_counts)])]

# Rename count columns with the corresponding sample names
colnames(combined_counts)[-1] <- sample_names

# Write the combined counts matrix to a CSV file
output_file <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Analysis/All_Combined_Counts.csv"
write.csv(combined_counts, file = output_file, row.names = FALSE)

# Message to confirm the output
message("All quant files processed successfully and combined counts saved to: ", output_file)

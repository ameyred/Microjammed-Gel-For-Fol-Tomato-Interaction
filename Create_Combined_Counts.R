#Written on 21st May 2025 for Micorgel project 
#Convert Quant.sf files to count matrix for Axenic vs Broth, Axenic vs Gel and Gel vs Broth. 



#Fol Axenic VS Broth

# Load necessary packages
library(tximport)
library(readr)

# Path to the transcript-to-gene mapping file
tx2gene_path <- "~/Amey_Lab/Microgel_Project/Fol_Genome/Transcript_To_Gene.csv"
tx2gene <- read.csv(tx2gene_path, stringsAsFactors = FALSE)

# Ensure the column names are correct
colnames(tx2gene) <- c("transcript_id", "gene_id")

# Print the first few rows of the transcript-to-gene mapping to check
print("Transcript-to-Gene Mapping:")
print(head(tx2gene))

# Set the directory where the quant.sf files are located
quant_dir <- "~/Amey_Lab/Microgel_Project/Axenic_VS_Broth/Quant_Files"

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
output_file <- "~/Amey_Lab/Microgel_Project/Axenic_VS_Broth/Combined_Counts.csv"
write.csv(combined_counts, file = output_file, row.names = FALSE)

# Message to confirm the output
cat("All quant files processed successfully and combined counts saved to: ", output_file)

#----------------------------------------------------------------------------------------------


#Fol Axenic VS Gel

# Load necessary packages
library(tximport)
library(readr)

# Path to the transcript-to-gene mapping file
tx2gene_path <- "~/Amey_Lab/Microgel_Project/Fol_Genome/Transcript_To_Gene.csv"
tx2gene <- read.csv(tx2gene_path, stringsAsFactors = FALSE)

# Ensure the column names are correct
colnames(tx2gene) <- c("transcript_id", "gene_id")

# Print the first few rows of the transcript-to-gene mapping to check
print("Transcript-to-Gene Mapping:")
print(head(tx2gene))

# Set the directory where the quant.sf files are located
quant_dir <- "~/Amey_Lab/Microgel_Project/Axenic_VS_Gel/Quant_Files"

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
output_file <- "~/Amey_Lab/Microgel_Project/Axenic_VS_Gel/Combined_Counts.csv"
write.csv(combined_counts, file = output_file, row.names = FALSE)

# Message to confirm the output
cat("All quant files processed successfully and combined counts saved to: ", output_file)

#----------------------------------------------------------------------------------------------



#Fol Broth VS Gel

# Load necessary packages
library(tximport)
library(readr)

# Path to the transcript-to-gene mapping file
tx2gene_path <- "~/Amey_Lab/Microgel_Project/Fol_Genome/Transcript_To_Gene.csv"
tx2gene <- read.csv(tx2gene_path, stringsAsFactors = FALSE)

# Ensure the column names are correct
colnames(tx2gene) <- c("transcript_id", "gene_id")

# Print the first few rows of the transcript-to-gene mapping to check
print("Transcript-to-Gene Mapping:")
print(head(tx2gene))

# Set the directory where the quant.sf files are located
quant_dir <- "~/Amey_Lab/Microgel_Project/Broth_VS_Gel/Quant_Files"

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
output_file <- "~/Amey_Lab/Microgel_Project/Broth_VS_Gel/Combined_Counts.csv"
write.csv(combined_counts, file = output_file, row.names = FALSE)

# Message to confirm the output
cat("All quant files processed successfully and combined counts saved to: ", output_file)

#----------------------------------------------------------------------------------------------


#Fol Broth VS Gel VS Axenic (for PCA and other purposed)

# Load necessary packages
library(tximport)
library(readr)

# Path to the transcript-to-gene mapping file
tx2gene_path <- "~/Amey_Lab/Microgel_Project/Fol_Genome/Transcript_To_Gene.csv"
tx2gene <- read.csv(tx2gene_path, stringsAsFactors = FALSE)

# Ensure the column names are correct
colnames(tx2gene) <- c("transcript_id", "gene_id")

# Print the first few rows of the transcript-to-gene mapping to check
print("Transcript-to-Gene Mapping:")
print(head(tx2gene))

# Set the directory where the quant.sf files are located
quant_dir <- "~/Amey_Lab/Microgel_Project/Quant_SF_Files"

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
output_file <- "~/Amey_Lab/Microgel_Project/All_Combined_Counts.csv"
write.csv(combined_counts, file = output_file, row.names = FALSE)

# Message to confirm the output
cat("All quant files processed successfully and combined counts saved to: ", output_file)

#----------------------------------------------------------------------------------------------








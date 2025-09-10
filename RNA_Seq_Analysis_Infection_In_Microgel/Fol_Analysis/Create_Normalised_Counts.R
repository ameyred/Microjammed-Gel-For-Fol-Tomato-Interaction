library(DESeq2)
library(readr)

# Load counts
counts_path <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Analysis/All_Combined_Counts.csv"
count_df <- read.csv(counts_path, row.names = 1, check.names = FALSE)

# Load sample info
sample_info_path <- "~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Analysis/All_Sample_Info.csv"
sample_info <- read.csv(sample_info_path, row.names = 1, check.names = FALSE)

# Convert counts to matrix and round
count_matrix <- round(as.matrix(count_df))

# Reorder sample info to match columns in count matrix
sample_info <- sample_info[colnames(count_matrix), , drop = FALSE]

# Sanity check
stopifnot(all(colnames(count_matrix) == rownames(sample_info)))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Treatment)

# Estimate size factors and extract normalized counts
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Save to CSV
norm_df <- as.data.frame(norm_counts)
norm_df$Gene <- rownames(norm_df)
norm_df <- norm_df[, c("Gene", setdiff(names(norm_df), "Gene"))]

write.csv(norm_df, file = "~/Amey_Lab/Microgel_Project/Big_Experiment/Fol_Analysis/Normalized_Counts.csv", row.names = FALSE)

message("✅ Normalized counts saved.")


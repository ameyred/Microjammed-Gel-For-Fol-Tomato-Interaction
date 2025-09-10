#This script was written on 21st May 2025 for DESeq on Syona's microgel samples comparing axenic Fol (from my bulk RNA) and fungus grown on Gamborg gel vs broth.



# Fol Axenic vs Broth
# Load libraries
library(DESeq2)
library(tidyverse)

# Load the counts file 
counts_data <- read.csv("~/Amey_Lab/Microgel_Project/Axenic_VS_Broth/Combined_Counts.csv" , header = FALSE)

# Set the first row as column names
colnames(counts_data) <- as.character(counts_data[1, ])
counts_data <- counts_data[-1, ]  # Remove the first row (which contains column headers)

# Set the first column as row names
rownames(counts_data) <- as.character(counts_data[[1]])
counts_data <- counts_data[-1]  # Remove the first column (gene names) from the data frame

# Convert the remaining data to numeric
counts_data_numeric <- apply(counts_data, 2, function(x) as.numeric(as.character(x)))

# Convert to a data frame
counts_data_table <- as.data.frame(counts_data_numeric)

# Ensure row names are preserved
rownames(counts_data_table) <- rownames(counts_data)


# Round counts to integers
counts_data_table <- round(counts_data_table)

# Display the final table and column names
print(head(counts_data_table))
colnames(counts_data_table)


#Sample info
colData <- read.csv("~/Amey_Lab/Microgel_Project/Axenic_VS_Broth/Sample_Info.csv", header = FALSE)
#Set the first row as column names
colnames(colData) <- as.character(colData[1, ])
colData <- colData[-1, ]  # Remove the first row
#Set the first column as row names
rownames(colData) <- as.character(colData[[1]])
colData <- colData[-1]  # Remove the first column
#Convert to a table
col_data_table <- as.data.frame(colData)
# Display the final table
print(col_data_table)
rownames(col_data_table)

#Check if all colnames in counts data are the same as rownames in col_data and in the same order. 
all(colnames(counts_data_table) %in% rownames(col_data_table))
all(colnames(counts_data_table) == rownames(col_data_table))

#DESeq Object Define
dds <- DESeqDataSetFromMatrix(countData = counts_data_table, 
                              colData = col_data_table, 
                              design = ~Treatment)
dds 

#Filtering rows with less than 10 number of reads (Recommended) 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#Set level, as in who to normalise against whom. 
dds$Treatment <- relevel(dds$Treatment, ref = "Axenic")
dds$Treatment

#Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
results_df <- as.data.frame(res)
output_file <- "~/Amey_Lab/Microgel_Project/Axenic_VS_Broth/DESEeq_Output.csv"
write.csv(results_df, file = output_file, row.names = TRUE)

#----------------------------------------------------------------------


# Fol Axenic vs Gel
# Load libraries
library(DESeq2)
library(tidyverse)

# Load the counts file 
counts_data <- read.csv("~/Amey_Lab/Microgel_Project/Axenic_VS_Gel/Combined_Counts.csv" , header = FALSE)

# Set the first row as column names
colnames(counts_data) <- as.character(counts_data[1, ])
counts_data <- counts_data[-1, ]  # Remove the first row (which contains column headers)

# Set the first column as row names
rownames(counts_data) <- as.character(counts_data[[1]])
counts_data <- counts_data[-1]  # Remove the first column (gene names) from the data frame

# Convert the remaining data to numeric
counts_data_numeric <- apply(counts_data, 2, function(x) as.numeric(as.character(x)))

# Convert to a data frame
counts_data_table <- as.data.frame(counts_data_numeric)

# Ensure row names are preserved
rownames(counts_data_table) <- rownames(counts_data)


# Round counts to integers
counts_data_table <- round(counts_data_table)

# Display the final table and column names
print(head(counts_data_table))
colnames(counts_data_table)


#Sample info
colData <- read.csv("~/Amey_Lab/Microgel_Project/Axenic_VS_Gel/Sample_Info.csv", header = FALSE)
#Set the first row as column names
colnames(colData) <- as.character(colData[1, ])
colData <- colData[-1, ]  # Remove the first row
#Set the first column as row names
rownames(colData) <- as.character(colData[[1]])
colData <- colData[-1]  # Remove the first column
#Convert to a table
col_data_table <- as.data.frame(colData)
# Display the final table
print(col_data_table)
rownames(col_data_table)

#Check if all colnames in counts data are the same as rownames in col_data and in the same order. 
all(colnames(counts_data_table) %in% rownames(col_data_table))
all(colnames(counts_data_table) == rownames(col_data_table))

#DESeq Object Define
dds <- DESeqDataSetFromMatrix(countData = counts_data_table, 
                              colData = col_data_table, 
                              design = ~Treatment)
dds 

#Filtering rows with less than 10 number of reads (Recommended) 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#Set level, as in who to normalise against whom. 
dds$Treatment <- relevel(dds$Treatment, ref = "Axenic")
dds$Treatment

#Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
results_df <- as.data.frame(res)
output_file <- "~/Amey_Lab/Microgel_Project/Axenic_VS_Gel/DESEeq_Output.csv"
write.csv(results_df, file = output_file, row.names = TRUE)



#----------------------------------------------------------------------

# Fol Broth vs Gel
# Load libraries
library(DESeq2)
library(tidyverse)

# Load the counts file 
counts_data <- read.csv("~/Amey_Lab/Microgel_Project/Broth_VS_Gel/Combined_Counts.csv" , header = FALSE)

# Set the first row as column names
colnames(counts_data) <- as.character(counts_data[1, ])
counts_data <- counts_data[-1, ]  # Remove the first row (which contains column headers)

# Set the first column as row names
rownames(counts_data) <- as.character(counts_data[[1]])
counts_data <- counts_data[-1]  # Remove the first column (gene names) from the data frame

# Convert the remaining data to numeric
counts_data_numeric <- apply(counts_data, 2, function(x) as.numeric(as.character(x)))

# Convert to a data frame
counts_data_table <- as.data.frame(counts_data_numeric)

# Ensure row names are preserved
rownames(counts_data_table) <- rownames(counts_data)


# Round counts to integers
counts_data_table <- round(counts_data_table)

# Display the final table and column names
print(head(counts_data_table))
colnames(counts_data_table)


#Sample info
colData <- read.csv("~/Amey_Lab/Microgel_Project/Broth_VS_Gel/Sample_Info.csv", header = FALSE)
#Set the first row as column names
colnames(colData) <- as.character(colData[1, ])
colData <- colData[-1, ]  # Remove the first row
#Set the first column as row names
rownames(colData) <- as.character(colData[[1]])
colData <- colData[-1]  # Remove the first column
#Convert to a table
col_data_table <- as.data.frame(colData)
# Display the final table
print(col_data_table)
rownames(col_data_table)

#Check if all colnames in counts data are the same as rownames in col_data and in the same order. 
all(colnames(counts_data_table) %in% rownames(col_data_table))
all(colnames(counts_data_table) == rownames(col_data_table))

#DESeq Object Define
dds <- DESeqDataSetFromMatrix(countData = counts_data_table, 
                              colData = col_data_table, 
                              design = ~Treatment)
dds 

#Filtering rows with less than 10 number of reads (Recommended) 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#Set level, as in who to normalise against whom. 
dds$Treatment <- relevel(dds$Treatment, ref = "Broth")
dds$Treatment

#Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
results_df <- as.data.frame(res)
output_file <- "~/Amey_Lab/Microgel_Project/Broth_VS_Gel/DESEeq_Output.csv"
write.csv(results_df, file = output_file, row.names = TRUE)



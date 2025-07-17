##### Load Libraries for Conversion of Entrez Gene Id to Gene symbol ######
source("/data/kumarr9/Priya_Project/nextflow/working_data/DGE/Load_package.R")
library(org.Hs.eg.db)
library(annotate)

### Converting the ENTREZ ID to Gene Symbol ######
library(org.Hs.eg.db)
library(annotate)

# 1. Read data and convert to character vector
a <- read.csv("/data/kumarr9/KIF18A/TCGA/LUAD_TCGA_entrez_ID.csv", header = TRUE)
entrez_ids <- as.character(a$Entrez_Gene_Id)  # Convert to character vector

# 2. Get gene symbols
gene_symbols <- getSYMBOL(entrez_ids, data='org.Hs.eg')

# 3. Combine results with original data
result <- data.frame(Entrez_Gene_Id = a$Entrez_Gene_Id,
                    Gene_Symbol = gene_symbols)

# 4. Save/output (optional)
write.csv(result, "/data/kumarr9/KIF18A/TCGA/LUAD_TCGA_Entrez_to_Symbol.csv", row.names = FALSE)

# Entrez Gene IDs to Ensembl IDs using org.Hs.eg.db in R:
library(org.Hs.eg.db)

# 1. Read your data (assuming same file structure)
a <- read.csv("/data/kumarr9/KIF18A/TCGA/LUAD_TCGA_entrez_ID.csv", header = TRUE)

# 2. Convert Entrez to Ensembl ID
result <- data.frame(
  Entrez_Gene_Id = a$Entrez_Gene_Id,
  Ensembl_ID = mapIds(org.Hs.eg.db,
                     keys = as.character(a$Entrez_Gene_Id),
                     column = "ENSEMBL",        # Key difference from SYMBOL
                     keytype = "ENTREZID",
                     multiVals = "first")       # Takes first match if multiple
)

# 3. Remove rows where Ensembl ID is NA (no mapping exists)
result <- result[!is.na(result$Ensembl_ID), ]

# 4. Save results
write.csv(result, "/data/kumarr9/KIF18A/TCGA/LUAD_TCGA_Entrez_to_Ensembl.csv", row.names = FALSE)

###### Working on the entire dataset and adding the Gene symbol and Ensemble_ID to the entire dataset #####
library(org.Hs.eg.db)
library(dplyr)
library(readr)

# 1. Read the expression data file
input_file <- "/data/kumarr9/KIF18A/TCGA/data_mrna_seq_fpkm_LUAD_TCGA.txt"
output_file <- "/data/kumarr9/KIF18A/TCGA/data_mrna_seq_fpkm_LUAD_TCGA_with_ids.tsv"

# Read the original data (tab-delimited, with header)
exp_data <- read.delim(input_file, check.names = FALSE)

# 2. Create mappings for both symbol and ensembl
gene_mappings <- data.frame(
  Entrez_Gene_Id = as.character(exp_data$Entrez_Gene_Id),
  Gene_symbol = mapIds(org.Hs.eg.db,
                      keys = as.character(exp_data$Entrez_Gene_Id),
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first"),
  Ensembl_ID = mapIds(org.Hs.eg.db,
                     keys = as.character(exp_data$Entrez_Gene_Id),
                     column = "ENSEMBL",
                     keytype = "ENTREZID",
                     multiVals = "first"),
  stringsAsFactors = FALSE
)

# 3. Merge with original data
# First extract the sample columns (all columns except Entrez_Gene_Id)
sample_cols <- setdiff(colnames(exp_data), "Entrez_Gene_Id")

# Create new dataframe with desired column order
final_data <- cbind(
  Entrez_Gene_Id = exp_data$Entrez_Gene_Id,
  gene_mappings[, c("Gene_symbol", "Ensembl_ID")],
  exp_data[, sample_cols]
)

# 4. Save as TSV
write_tsv(final_data, output_file)

# Print confirmation
message("File saved successfully at:\n", output_file)
message("\nNew columns added between Entrez_Gene_Id and sample columns.")

###### LUSC Dataset #########
library(org.Hs.eg.db)
library(dplyr)
library(readr)

# 1. Read the expression data file
input_file <- "/data/kumarr9/KIF18A/TCGA/data_mrna_seq_fpkm_LUSC_TCGA.txt"
output_file <- "/data/kumarr9/KIF18A/TCGA/data_mrna_seq_fpkm_LUSC_TCGA_with_ids.tsv"

# Read the original data (tab-delimited, with header)
exp_data <- read.delim(input_file, check.names = FALSE)

# 2. Create mappings for both symbol and ensembl
gene_mappings <- data.frame(
  Entrez_Gene_Id = as.character(exp_data$Entrez_Gene_Id),
  Gene_symbol = mapIds(org.Hs.eg.db,
                      keys = as.character(exp_data$Entrez_Gene_Id),
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first"),
  Ensembl_ID = mapIds(org.Hs.eg.db,
                     keys = as.character(exp_data$Entrez_Gene_Id),
                     column = "ENSEMBL",
                     keytype = "ENTREZID",
                     multiVals = "first"),
  stringsAsFactors = FALSE
)

# 3. Merge with original data
# First extract the sample columns (all columns except Entrez_Gene_Id)
sample_cols <- setdiff(colnames(exp_data), "Entrez_Gene_Id")

# Create new dataframe with desired column order
final_data <- cbind(
  Entrez_Gene_Id = exp_data$Entrez_Gene_Id,
  gene_mappings[, c("Gene_symbol", "Ensembl_ID")],
  exp_data[, sample_cols]
)

# 4. Save as TSV
write_tsv(final_data, output_file)

# Print confirmation
message("File saved successfully at:\n", output_file)
message("\nNew columns added between Entrez_Gene_Id and sample columns.")


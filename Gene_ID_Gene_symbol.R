# Load libraries
library(dplyr)
library(biomaRt)

# Read your expression data
expr_df <- read.table("/data/kumarr9/SCLC_Jiang/GSE60052_raw_counts_GRCh38.p13_NCBI.tsv",
                      sep = "\t", header = TRUE, check.names = FALSE)

# Extract Entrez Gene IDs
entrez_ids <- as.character(expr_df$GeneID)

# Use biomaRt to map Entrez IDs to gene symbols
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve mapping: Entrez → HGNC symbol
id_mapping <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = entrez_ids,
  mart = mart
)

# Merge with expression matrix
expr_df_annotated <- expr_df %>%
  left_join(id_mapping, by = c("GeneID" = "entrezgene_id")) %>%
  relocate(hgnc_symbol, .after = GeneID)

# Optional: filter out rows with no gene symbol
expr_df_annotated <- expr_df_annotated %>% filter(hgnc_symbol != "")

# Save output
write.table(expr_df_annotated,
            file = "/data/kumarr9/SCLC_Jiang/GSE60052_with_symbols.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

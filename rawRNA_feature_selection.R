# This script performs differential expression analysis on RNA-seq raw counts
# using DESeq2, focusing on cases with RCB assessment and more than one cycle
# of therapy, as per the study design detailed in Methods (Statistical testing).
# It assumes the raw counts data is in a gzipped tab-separated file and metadata
# is in an Excel file.

# Load necessary libraries
# Install required packages if not already installed
packages <- c("data.table", "readxl")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}
# Bioconductor packages like DESeq2 may need BiocManager
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("DESeq2")
}
# glmnet for LASSO regression
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
# Ensure biomaRt is installed for gene annotation
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("biomaRt")
}
# EnhancedVolcano for visualization
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("EnhancedVolcano")
}

library(data.table)
library(readxl)
library(DESeq2)
library(glmnet)
library(biomaRt)
library(EnhancedVolcano)

# Load the data
# Assuming the file is in the current working directory or specify the full path
file_path <- "data/transneo-diagnosis-RNAseq-rawcounts.tsv.gz"
# Use fread from data.table for fast reading of gzipped tab-separated files
df <- fread(file_path, header = TRUE, sep = "\t", data.table = FALSE)
row.names(df) <- df[[1]]  # Set the first column as row names
df <- df[, -1, drop = FALSE]  # Remove the first column after setting row names

# Check if the data is loaded correctly
if (is.null(df) || nrow(df) == 0) {
  stop("Data loading failed or the file is empty.")
}

# Load metadata
# Assuming the metadata file is in the same directory as the raw counts file
metadata <- data.frame(read_excel("data/Supplementary_table.xlsx", sheet = 1))

# Perform analyses with cases that had an RCB assessment and received more than
# one cycle of therapy
# Filter metadata to include only cases with RCB assessment and more than one
# cycle of therapy
# and HER2-targeted therapy, as per the study design
# as detailed in Methods (Statistical testing)
metadata <- metadata[metadata$RCB.category != "NA", ]
metadata <- metadata[metadata$Chemo.cycles > 1 & metadata$aHER2.cycles > 1, ]

# Create a binary “response” column
metadata$Response <- ifelse(metadata$RCB.category == "pCR", "pCR", "RD")
# When build colData, it will carry that column through:
common   <- intersect(colnames(df), metadata$Donor.ID)
count_data <- df[, common, drop = FALSE]
col_data   <- metadata[match(common, metadata$Donor.ID), , drop = FALSE]

# Differential expression with DESeq2
# Input: un-normalized counts
dds <- DESeqDataSetFromMatrix(
  countData = count_data,   # raw counts
  colData   = col_data,     # metadata
  design    = ~Response     # design formula: Response as a factor
)

# Run DESeq
dds <- DESeq(dds)

# Check what you can test
resultsNames(dds)

# res <- results(dds, contrast = c("Response", "pCR", "RD"))
res <- results(dds)
summary(res)

# Filter results for significant genes
# Adjusted p-value < 0.05 and absolute log2 fold change > 0
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0)
# View top significant genes
head(sig_genes[order(sig_genes$padj), ])

# Prepare data: expression matrix of significant genes
gene_mat <- t(count_data[rownames(sig_genes), ])
response <- as.factor(col_data$Response)

# LASSO logistic regression
set.seed(123)
cvfit <- cv.glmnet(gene_mat, response, family = "binomial", alpha = 1)

# Get non-zero coefficients (selected predictors)
best_lambda <- cvfit$lambda.min
coef_lasso <- coef(cvfit, s = best_lambda)
selected_genes <- rownames(coef_lasso)[which(coef_lasso != 0)]
selected_genes <- selected_genes[selected_genes != "(Intercept)"]  # Remove intercept

# View top predictors
print(selected_genes)

# Assume your gene IDs are Ensembl IDs in rownames(sig_genes)
gene_ids <- rownames(sig_genes)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


# Connect to Ensembl BioMart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve annotations: gene symbol, description, etc.
annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart
)

# Merge annotation with your results
annotated_sig_genes <- merge(
  data.frame(ensembl_gene_id = gene_ids, sig_genes),
  annotations,
  by = "ensembl_gene_id",
  all.x = TRUE
)

# View annotated results
head(annotated_sig_genes)
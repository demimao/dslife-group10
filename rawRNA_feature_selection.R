# This script performs differential expression analysis on RNA-seq raw counts
# using DESeq2, focusing on cases with RCB assessment and more than one cycle
# of therapy, as per the study design detailed in Methods (Statistical testing).
# It assumes the raw counts data is in a gzipped tab-separated file and metadata
# is in an Excel file.

# Load necessary libraries
# Install required CRAN packages if not already installed
cran_packages <- c("data.table", "readxl", "glmnet")
installed <- cran_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(cran_packages[!installed])
}

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages if not already installed
bioc_packages <- c("DESeq2", "biomaRt")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

library(data.table)
library(readxl)
library(DESeq2)
library(glmnet)
library(biomaRt)
library(dplyr)

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
cvfit <- cv.glmnet(gene_mat, response, family = "binomial", alpha = 1, standardize = TRUE)

# Get non-zero coefficients (selected predictors)
best_lambda <- cvfit$lambda.min
coef_lasso <- coef(cvfit, s = best_lambda)
# convert to a named numeric vector and drop the intercept
coef_vec <- as.numeric(coef_lasso)
names(coef_vec) <- rownames(coef_lasso)
coef_vec <- coef_vec[names(coef_vec) != "(Intercept)"]

# select only those genes with non-zero coefficient
selected_genes <- names(coef_vec)[coef_vec != 0]

# Filter the matrix
important_genes <- sig_genes[selected_genes , , drop=FALSE]
# Add lasso coefficients
important_genes$lasso_coef <- coef_vec[ rownames(important_genes) ]
head(important_genes)

# Assume your gene IDs are Ensembl IDs in rownames(sig_genes)
gene_ids <- rownames(important_genes)

# Connect to Ensembl BioMart
# Or use a recent archived release if the current site is down:
mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host    = "https://oct2024.archive.ensembl.org/"
)

# Retrieve annotations: gene symbol, description, etc.
annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart
)

# Merge annotation with your results
annotated_sig_genes <- merge(
  data.frame(ensembl_gene_id = gene_ids, important_genes),
  annotations,
  by = "ensembl_gene_id"
)

# View annotated results
head(annotated_sig_genes)

# 1. Define “pCR” groups
res <- annotated_sig_genes
res$group_pcr <- "Not significant"
res$group_pcr[res$padj < 0.05 & res$log2FoldChange >  0] <- "Overexpressed in pCR"
res$group_pcr[res$padj < 0.05 & res$log2FoldChange <  0] <- "Underexpressed in pCR"
# Extract over-underexpressed genes
# 5 most under-expressed (smallest log2FC)
top5_under <- res %>% 
  slice_min(order_by = log2FoldChange, n = 5)

# 5 most over-expressed (largest log2FC)
top5_over  <- res %>% 
  slice_max(order_by = log2FoldChange, n = 5)

# Combine them into one table
top10 <- bind_rows(
  top5_under  %>% mutate(direction = "under"),
  top5_over   %>% mutate(direction = "over")
)
# This script performs differential expression analysis on RNA-seq raw counts
# using DESeq2, focusing on cases with RCB assessment and more than one cycle
# of therapy, as per the study design detailed in Methods (Statistical testing).
# It assumes the raw counts data is in a gzipped tab-separated file and metadata
# is in an Excel file.

# Load necessary libraries
# Install required packages if not already installed
packages <- c("data.table", "readxl", "DESeq2")
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

library(data.table)
library(readxl)
library(DESeq2)

# Load the data
# Assuming the file is in the current working directory or specify the full path
file_path <- "data/transneo-diagnosis-RNAseq-rawcounts.tsv.gz"
# Use fread from data.table for fast reading of gzipped tab-separated files
df <- fread(file_path, header = TRUE, sep = "\t", data.table = FALSE)

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

# Ensure that the metadata and raw counts data are aligned
# Check if the Donor.ID in metadata matches the sample names in df
samples_count <- colnames(df)[-1]
samples_meta <- metadata$Donor.ID

# Create a binary “response” column
metadata$Response <- ifelse(metadata$RCB.category == "pCR", "pCR", "RD")
# When build colData, it will carry that column through:
common   <- intersect(colnames(df)[-1], metadata$Donor.ID)
count_data <- df[, common, drop = FALSE] # nolint
col_data   <- metadata[match(common, metadata$Donor.ID), , drop = FALSE]

# Differential expression with DESeq2
# Input: un-normalized counts
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = col_data,
  design    = ~ Response
)

# pre‐filter: drop genes with very low counts across all samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

# Run DESeq
dds <- DESeq(dds)

# Check what you can test
resultsNames(dds)

res <- results(dds, contrast = c("response", "pCR", "RD"))
summary(res)
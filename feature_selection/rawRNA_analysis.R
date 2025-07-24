# Authors: Nikolai Egorov, Zhixin Mao, Mingxuan Fan, Group 10
# This R script implements a complete pipeline to identify gene signatures 
# predictive of neoadjuvant therapy response in breast cancer by loading raw 
# RNA‑seq counts and clinical metadata (RCB status, ER/HER2, chemo cycles),
# filtering for cases with >1 therapy cycle and RCB assessment, performing 
# DESeq2 differential expression (|log₂FC|>1, padj<0.05), applying a 
# variance‑stabilizing transform followed by LASSO logistic regression to select 
# the top 20 predictive genes, annotating these genes via Ensembl BioMart 
# (and cross‑referencing with IntOGen BRCA drivers), and generating 
# publication‑quality visualizations (p‑value histogram, heatmap, waterfall plot, 
# PCA scatter, and volcano plot).
# It assumes the raw counts data is in a gzipped tab-separated file and metadata
# is in an Excel file.
# References:
#   – Sammut, SJ., Crispin-Ortuzar, M., Chin, SF. et al. 
#     Multi-omic machine learning predictor of breast cancer therapy response.
#     https://www.nature.com/articles/s41586-021-04278-5
#   – GitHub repo:
#     https://github.com/cclab-brca/neoadjuvant-therapy-response-predictor
#   – IntOGen BRCA drivers:
#     https://www.intogen.org/search?cancer=BRCA

#-------------------------------------------------------------------------------
# Preparation
#-------------------------------------------------------------------------------

# 1. Install and load necessary libraries

# Install required CRAN packages if not already installed
#cran_packages <- c("data.table", "readxl", "glmnet")
#installed <- cran_packages %in% rownames(installed.packages())
#if (any(!installed)) {
#  install.packages(cran_packages[!installed])
#}
# Ensure BiocManager is installed
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}
# Install Bioconductor packages if not already installed
#bioc_packages <- c("DESeq2", "biomaRt")
#for (pkg in bioc_packages) {
#  if (!requireNamespace(pkg, quietly = TRUE)) {
#    BiocManager::install(pkg)
#  }
#}

# Load necessary libraries
library(data.table)
library(readxl)
library(DESeq2)
library(glmnet)
library(biomaRt)
library(dplyr)
library(tibble)
library(MASS)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

# 2. Load the data
# Assuming the file is in the current working directory or specify the full path
file_path <- "data/transneo-diagnosis-RNAseq-rawcounts.tsv.gz"
# Use fread from data.table for fast reading of gzipped tab-separated files
df <- fread(file_path, header = TRUE, sep = "\t", data.table = FALSE)
# Set the first column as row names
row.names(df) <- df[[1]]  
df <- df[, -1, drop = FALSE]  # Remove the first column after setting row names

# Check if the data is loaded correctly
if (is.null(df) || nrow(df) == 0) {
  stop("Data loading failed or the file is empty.")
}

# 3. Load metadata
# Assuming the metadata file is in the same directory as the raw counts file
metadata <- data.frame(read_excel("data/Supplementary_table .xlsx", sheet = 1))
# combine ER and HER2 status
metadata$ERHER2.status <- ifelse(metadata$ER.status=="POS","ER+ HER2-","ER- HER2-")
metadata$ERHER2.status <- ifelse(metadata$HER2.status=="POS","HER2+",metadata$ERHER2.status)

# 4. Filter metadata to include only cases with RCB assessment and more than one
# cycle of therapy and HER2-targeted therapy.
metadata <- metadata[!is.na(metadata$RCB.category), ]
metadata <- metadata[metadata$Chemo.cycles > 1 & metadata$aHER2.cycles > 1, ]

# Create a binary “response” column
metadata$Response <- ifelse(metadata$RCB.category == "pCR", "pCR", "RD")

# 5. update metadata - retain only samples that have RNAseq data
p <- metadata
p <- p[p$Donor.ID %in% colnames(df),]

# 6. Filter the counts to only those samples, and reorder to match p
df <- df[ , p$Donor.ID, drop = FALSE]

# (Optional) sanity-check that ordering matches
all(colnames(df) == p$Donor.ID)
# Should return TRUE

#-------------------------------------------------------------------------------
# Feature selection
#-------------------------------------------------------------------------------

# 1. Differential expression analysis with DESeq2
# When build colData, it will carry that column through:
common     <- p$Donor.ID
count_data <- df[, common, drop = FALSE]
col_data   <- p
col_data$Response <- factor(col_data$Response, levels = c("RD","pCR"))

# Input: un-normalized counts
dds <- DESeqDataSetFromMatrix(
  countData = count_data,   # raw counts
  colData   = col_data,     # metadata
  design    = ~Response     # design formula: Response as a factor
)

# Run DESeq
dds <- DESeq(dds)

# See the result
res <- results(dds)
summary(res)

# Filter result for significant genes
# Adjusted p-value < 0.05 and absolute log2 fold change > 1
de_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1 & baseMean > 10)
de_genes_test <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
# View top significant genes
head(de_genes[order(de_genes$padj), ])

# 2. Prepare data: expression matrix of significant genes
# Variance‐stabilizing transform (removes depth effects & stabilizes variance)
vsd <- vst(dds, blind = FALSE)

# Pull out the normalized matrix (genes × samples)
norm_mat <- assay(vsd)
norm_de_genes <- norm_mat[rownames(de_genes),]

# gene_mat <- t(count_data[rownames(de_genes), ])
gene_mat <- t(norm_de_genes)
response <- as.factor(col_data$Response)

# 3. Cross-validation LASSO logistic regression
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
lasso_genes <- de_genes[selected_genes , , drop=FALSE]
# Add lasso coefficients
lasso_genes$lasso_coef <- coef_vec[ rownames(lasso_genes) ]

# 4. Filter dataset, include genes that have relatively high padj and lasso coef.
# Turn our lasso'ed dataset of DESeqResults‐like object into a plain data frame
lasso_df <- as.data.frame(lasso_genes) %>%
  rownames_to_column("geneID") %>%
  dplyr::select(geneID, log2FoldChange, padj, lasso_coef)

# Now compute ranks and pick the top 20
top20 <- lasso_df %>%
  filter(!is.na(padj)) %>%               # drop any NAs just in case
  mutate(
    #rank_DE    = rank(padj,             ties.method = "first"),
    rank_LASSO = rank(-abs(lasso_coef), ties.method = "first"),
    #rank_sum   = rank_DE + rank_LASSO
  ) %>%
  arrange(rank_LASSO) %>%
  slice_head(n = 20)

feature_mat <- t(norm_de_genes[top20$geneID, , drop = FALSE ])

# 5. Write features to our feature table.

# Read both tables
#rna_features_stat  <- read.csv("data/training_set_final.csv",  stringsAsFactors = FALSE)

#rna_features_stat <- rna_features_stat[, c("Trial.ID", "STAT1.gsva", "GGI.gsva", "ESC.gsva")]
#colnames(rna_features_stat)[1] <- "Donor.ID"
#rna_features  <- as.data.frame(feature_mat) %>%
#  rownames_to_column(var = "Donor.ID")
#training_set <- data.frame(read_excel("data/training_set_new.xlsx", sheet = 1))

#merged <- training_set %>%
#  left_join(rna_features, by = "Donor.ID")

#merged <- merged %>%
#  left_join(rna_features_stat, by = "Donor.ID")

# Write merged dataset into new csv file.
#write.csv(merged,
#          file      = "data/training_set_new.csv",
#          row.names = FALSE)

#-------------------------------------------------------------------------------
# Annotation
#-------------------------------------------------------------------------------

# Extract gene ids from top20 selected genes
gene_ids <- rownames(t(feature_mat))

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
annotated_genes_count <- merge(
  data.frame(ensembl_gene_id = gene_ids, t(feature_mat)),
  annotations,
  by = "ensembl_gene_id"
)

lasso_annotated_genes <- merge(
  data.frame(ensembl_gene_id = gene_ids, top20),
  annotations,
  by = "ensembl_gene_id"
)

# View annotated results
head(annotated_genes_count)

#-------------------------------------------------------------------------------
# Plot the results
#-------------------------------------------------------------------------------
# 0. p-value distribution
ggplot(as.data.frame(res), aes(pvalue)) +
  geom_histogram(bins=50, fill="steelblue", color="white") +
  labs(
    title="P-value distribution",
    x="Raw p-value", y="Count"
  )

# 1. Heatmap
heatmap_df <- annotated_genes_count %>%
  # keep only the sample columns + gene symbol
  dplyr::select(hgnc_symbol, starts_with("T")) %>%
  pivot_longer(
    cols      = -hgnc_symbol,
    names_to  = "Donor.ID",
    values_to = "Expression"
  )

# Heatmap with counts
ggplot(heatmap_df, aes(x = Donor.ID, y = hgnc_symbol, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(
    low    = "blue",
    mid    = "white",
    high   = "red",
    midpoint = median(heatmap_df$Expression, na.rm = TRUE),
    name   = "VST\nexpression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title   = element_blank(),
    panel.grid   = element_blank()
  ) +
  labs(
    title = "Expression of Genes Across Samples",
    subtitle = "Rows = genes; columns = donors"
  ) + scale_x_discrete(
    breaks = heatmap_df$Donor.ID[seq(1, length(unique(heatmap_df$Donor.ID)), by = 5)]
  )

# 2. Waterfall plot
# Order by log₂FC
# First, drop any “novel transcript” rows with no symbol
lasso_annotated_genes <- lasso_annotated_genes %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "")

# Order by log₂FC
lasso_annotated_genes <- lasso_annotated_genes %>%
  arrange(log2FoldChange)

# Extract the unique symbols in that order
symbol_levels <- unique(lasso_annotated_genes$hgnc_symbol)

# Re-factor and add Direction
lasso_annotated_genes <- lasso_annotated_genes %>%
  mutate(
    hgnc_symbol = factor(hgnc_symbol, levels = symbol_levels),
    Direction   = ifelse(log2FoldChange > 0, "Over", "Under")
  )
ggplot(lasso_annotated_genes, aes(x = hgnc_symbol, y = log2FoldChange, fill = Direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(Over="#D73027", Under="#4575B4")) +
  labs(
    title = "Top 14 Genes: Log₂ Fold Change",
    x     = "Gene",
    y     = "log₂(Fold Change)",
    fill  = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  )

# 3. PCA 
# On whole set of DE-genes
plotPCA(vsd, intgroup = "Response") +
  ggtitle("PCA of VST‑transformed counts")

# On top 20 genes
top_genes   <- top20$geneID
vsd_sub     <- norm_mat[top_genes, , drop = FALSE]    # only those 20 rows

# Transpose and run PCA
pca_res <- prcomp(t(vsd_sub), scale. = TRUE)

# Grab the scores for samples and merge in your metadata
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("Donor.ID") %>%
  left_join(p[, c("Donor.ID","Response")], by = "Donor.ID")

# Plot PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2, color = Response)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Response), geom = "polygon", alpha = 0.2, color = NA) +
  labs(
    title    = "PCA on Top 20 Differential Genes",
    subtitle = paste0("Explained: PC1=", round(100*summary(pca_res)$importance[2,1],1),
                      "%, PC2=", round(100*summary(pca_res)$importance[2,2],1), "%"),
    x        = paste0("PC1 (", round(100*summary(pca_res)$importance[2,1],1), "%)"),
    y        = paste0("PC2 (", round(100*summary(pca_res)$importance[2,2],1), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title      = element_text(face = "bold")
  )

# 4.Volcano plot
# Prepare a data frame (if needed)
volcano_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("Gene.ID")

# Plot with EnhancedVolcano
EnhancedVolcano(
  volcano_df,
  lab           = volcano_df$Gene.ID,         # point labels
  x             = 'log2FoldChange',           # column name for x-axis
  y             = 'padj',                     # column name for y-axis
  xlim          = c(-5, 5),                   # adjust as needed
  pCutoff       = 0.05,                       # FDR threshold
  FCcutoff      = 1,                          # |log2FC| threshold
  pointSize     = 2,
  selectLab      = character(),
  col           = c('grey30', 'firebrick2', 'royalblue2', 'darkorange2'),
  legendLabels  = c('NS','|log2FC| > 1','padj < 0.05','padj < 0.05 & |log2FC| > 1'),
  legendPosition = 'right',
  title         = 'Volcano Plot: pCR vs RD Response',
  subtitle      = 'DESeq2 differential expression',
  caption       = 'Cutoffs: |log₂FC| > 1, padj < 0.05',
  drawConnectors = TRUE,                      # connector lines to labels
  max.overlaps = 20
)

#-------------------------------------------------------------------------------
# Annotation of driver genes
#-------------------------------------------------------------------------------
# Read breast cancer driver genes table
driver_genes <- fread("data/IntOGen-DriverGenes_BRCA.tsv", sep = "\t", header = TRUE)
driver_syms <- driver_genes$Symbol

# Then filter
annotated_drivers <- annotated_genes_count[
  annotated_genes_count$hgnc_symbol %in% driver_syms,
]

# Quick check
table(annotated_drivers$hgnc_symbol %in% driver_syms)  # should all be TRUE

# Result: 
# LASSO could not highlight any driver gene as significant in this dataset.

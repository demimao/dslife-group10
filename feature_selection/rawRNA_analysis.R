# Authors: Nikolai Egorov, Zhixin Mao, Mingxuan Fan from Group 10
# This script performs differential expression analysis on RNA-seq raw counts
# using DESeq2, focusing on cases with RCB assessment and more than one cycle
# of therapy, as per the study design detailed in Methods (Statistical testing).
# Furthermore, we performed association testing given features from 
# Supplementary table (sheet 3).
# It assumes the raw counts data is in a gzipped tab-separated file and metadata
# is in an Excel file.
# This project is based on the existing study: 
# https://www.nature.com/articles/s41586-021-04278-5#Fig9
# Data for analysis is provided in the github: 
# https://github.com/cclab-brca/neoadjuvant-therapy-response-predictor/tree/master

#-------------------------------------------------------------------------------
# Preparation
#-------------------------------------------------------------------------------

# 1. Install and load necessary libraries
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
metadata <- data.frame(read_excel("data/Supplementary_table.xlsx", sheet = 1))
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
# Association testing
#-------------------------------------------------------------------------------
# Load RNA data
rnadata <- data.frame(read_excel("data/Supplementary_table.xlsx", sheet = 3))

# a. GGI gsva score

rnadata <- merge(rnadata,p[,c("Donor.ID","RCB.category","pCR.RD","ER.status","HER2.status","ERHER2.status","Grade.pre.NAT")],
                 by="Donor.ID")

rnadata$RCB.category <- factor(rnadata$RCB.category,
                               levels=c("pCR","RCB-I","RCB-II","RCB-III"),
                               ordered=T)

# lm() for a continuous‐on‐continuous relationship
summary(lm(rnadata$GGI.gsva~rnadata$Grade.pre.NAT))

# GGI score is monotonically associated with RCB
# p=2e-05
# polr() for a continuous‐on‐ordered categorical relationship
m      <- polr(factor(rnadata$RCB.category,ordered = T) ~ rnadata$GGI.gsva,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
pval[1]

z <- rnadata
z$RCB.category< - factor(z$RCB.category,levels=c("RCB-III","RCB-II","RCB-I","pCR"),labels=c("III","II","I","0"))

# GGI score only associated with response in HER2- tumours, p=4e-05
wilcox.test(
  rnadata[rnadata$HER2.status=="NEG" & rnadata$pCR.RD=="pCR","GGI.gsva"],
  rnadata[rnadata$HER2.status=="NEG" & rnadata$pCR.RD=="RD","GGI.gsva"])
# HER2+ p=0.99
wilcox.test(
  rnadata[rnadata$HER2.status=="POS" & rnadata$pCR.RD=="pCR","GGI.gsva"],
  rnadata[rnadata$HER2.status=="POS" & rnadata$pCR.RD=="RD","GGI.gsva"])

wilcox.test(
  rnadata[rnadata$RCB.category=="pCR","GGI.gsva"],
  rnadata[rnadata$RCB.category=="RCB-II","GGI.gsva"])

# b. Associations with Embryonic stem cell score

# ESC score is monotonically associated with the degree of RD, p=0.0001
m      <- polr(rnadata$RCB.category ~ rnadata$ESC.gsva,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
pval[1]

table(rnadata$ERHER2.status)

# c. Associations with Taxane final score
pacli <- rnadata[, grep("Donor.ID|Taxane",colnames(rnadata))]
pacli <- merge(pacli,p,by="Donor.ID")
pacli$HER2.status <- factor(pacli$HER2.status,levels=c("NEG","POS"),labels=c("HER2-","HER2+"))

wilcox.test(pacli[pacli$HER2.status=="HER2-" & pacli$pCR.RD=="pCR","Taxane.FinalScore",],
            pacli[pacli$HER2.status=="HER2-" & pacli$pCR.RD!="pCR","Taxane.FinalScore",])

# d. Association with every other score

features <- c("GGI.gsva", "ESC.gsva", "STAT1.gsva", "CytScore.log2", "Danaher.B.cells", "Danaher.CD45",
              "Danaher.CD8.T.cells", "Danaher.Cytotoxic.cells", "Danaher.DC", "Danaher.Exhausted.CD8", 
              "Danaher.Macrophages", "Danaher.Mast.cells" , "Danaher.Neutrophils", "Danaher.NK.CD56dim.cells",
              "Danaher.NK.cells", "Danaher.T.cells", "Danaher.Th1.cells", "Danaher.Treg",
              "TIDE.Dysfunction", "TIDE.Exclusion", "TIDE.MDSC", "TIDE.CAF",
              "TIDE.TAM.M2", "TIDE.IFNG", "TIDE.CD274", "TIDE.CD8")

univ_results <- lapply(features, function(f) {
  form       <- as.formula(paste("RCB.category ~", f))
  m          <- polr(form, data=rnadata, Hess=TRUE)
  s          <- summary(m)
  # assume the predictor is in the last row of the coefficient table
  coef_tbl   <- coef(s)
  # which row corresponds to the predictor?
  # typically it’s the last one
  slope_row  <- nrow(coef_tbl)
  estimate   <- coef_tbl[slope_row, "Value"]
  std_error  <- coef_tbl[slope_row, "Std. Error"]
  # two‐sided Wald p‐value
  p_raw      <- pnorm(abs(coef_tbl[, "t value"]), lower.tail = FALSE) * 2
  
  data.frame(
    feature   = f,
    estimate  = estimate,
    std.error = std_error,
    p.value   = p_raw[1],
    stringsAsFactors = FALSE
  )
})

# Bind into one data.frame
results_df <- bind_rows(univ_results)

# Adjust for multiple testing
results_df <- results_df %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH")   # only FDR
  ) %>%
  arrange(p.adj)

# Pick candidates at FDR < 5%
candidates <- filter(results_df, p.adj < 0.05)$feature
# Inspect
print(candidates)

# Read both tables
#training_df  <- read.csv("data/training_df.csv",  stringsAsFactors = FALSE)
#training_set <- read_excel("data/training_set.csv", sheet = 1)

#selected_feats <- training_df[, c("Trial.ID", candidates), drop = FALSE]

# merged <- training_set %>%
#   left_join(selected_feats, by = "Trial.ID")

# Write merged dataset into new csv file.
# write.csv(merged,
#          file      = "data/training_set_final.csv",
#          row.names = FALSE)

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
# Adjusted p-value < 0.05 and absolute log2 fold change > 0
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

# Cross-validation LASSO logistic regression
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

# 3. Filter dataset, include genes that have relatively high padj and lasso coef.
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

# Waterfall plot
# Order by log₂FC
# 0) First, drop any “novel transcript” rows with no symbol
lasso_annotated_genes <- lasso_annotated_genes %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "")

# 1) Order by log₂FC
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
    title = "Top 20 Genes: Log₂ Fold Change",
    x     = "Gene",
    y     = "log₂(Fold Change)",
    fill  = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  )

# PCA on top 20 genes
top_genes   <- top20$geneID
vsd_sub     <- norm_mat[top_genes, , drop = FALSE]    # only those 20 rows

# 3. Transpose and run PCA
pca_res <- prcomp(t(vsd_sub), scale. = TRUE)

# 4. Grab the scores for samples and merge in your metadata
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("Donor.ID") %>%
  left_join(p[, c("Donor.ID","Response")], by = "Donor.ID")

# 5. Plot PC1 vs PC2
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

# 3. Prepare a data frame (if needed)
volcano_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("Gene.ID")

# 4. Plot with EnhancedVolcano
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

# As we can see, LASSO could not highlight any driver gene as significant in this dataset.

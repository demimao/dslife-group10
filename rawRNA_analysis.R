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
library(MASS)

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
# Load RNA data
# Includes ....
rnadata <- data.frame(read_excel("data/Supplementary_table.xlsx", sheet = 3))

# combine ER and HER2 status
metadata$ERHER2.status <- ifelse(metadata$ER.status=="POS","ER+ HER2-","ER- HER2-")
metadata$ERHER2.status <- ifelse(metadata$HER2.status=="POS","HER2+",metadata$ERHER2.status)

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

p <- metadata
# update metadata - retain only samples that have RNAseq data
p <- p[p$Donor.ID %in% colnames(df),]

#-------------------------------------------------------------------------------
# Association testing
#-------------------------------------------------------------------------------
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

# Prepare data: expression matrix of significant genes
#gene_mat <- t(count_data[rownames(sig_genes), ])
#response <- as.factor(col_data$Response)

# LASSO logistic regression
# set.seed(123)
# cvfit <- cv.glmnet(gene_mat, response, family = "binomial", alpha = 1, standardize = TRUE)

# Get non-zero coefficients (selected predictors)
#best_lambda <- cvfit$lambda.min
#coef_lasso <- coef(cvfit, s = best_lambda)
# convert to a named numeric vector and drop the intercept
#coef_vec <- as.numeric(coef_lasso)
#names(coef_vec) <- rownames(coef_lasso)
#coef_vec <- coef_vec[names(coef_vec) != "(Intercept)"]

# select only those genes with non-zero coefficient
#selected_genes <- names(coef_vec)[coef_vec != 0]

# Filter the matrix
#important_genes <- sig_genes[selected_genes , , drop=FALSE]
# Add lasso coefficients
#important_genes$lasso_coef <- coef_vec[ rownames(important_genes) ]
#head(important_genes)

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

# Assume your gene IDs are Ensembl IDs in rownames(sig_genes)
gene_ids <- rownames(sig_genes)

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
  data.frame(ensembl_gene_id = gene_ids, sig_genes),
  annotations,
  by = "ensembl_gene_id"
)

# View annotated results
head(annotated_sig_genes)

# Read breast cancer driver genes table
driver_genes <- fread("data/IntOGen-DriverGenes_BRCA.tsv", sep = "\t", header = TRUE)
driver_syms <- driver_genes$Symbol

# Then filter
annotated_drivers <- annotated_sig_genes[
  annotated_sig_genes$hgnc_symbol %in% driver_syms,
]

# Quick check
table(annotated_drivers$Symbol %in% driver_syms)  # should all be TRUE

# 1. Define “pCR” groups
res <- annotated_sig_genes
res$group_pcr <- "Not significant"
res$group_pcr[res$padj < 0.05 & res$log2FoldChange >  0] <- "Overexpressed in pCR"
res$group_pcr[res$padj < 0.05 & res$log2FoldChange <  0] <- "Underexpressed in pCR"
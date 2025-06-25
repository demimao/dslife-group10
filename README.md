# Interpretable Prediction of Breast Cancer Therapy Response Using Shapley-Based Explanations

**Authors**: Mingxuan Fan, Nikolai Egorov, Zhixin Mao  
**Project Type**: Scientific Research / Machine Learning in Bioinformatics  
**Keywords**: Breast Cancer, Multi-omics, Machine Learning, SHAP, pCR, Interpretation, Neoadjuvant Therapy 

> *Inspired by*: Sammut, S.J., Crispin-Ortuzar, M., Chin, S.F. et al. (2022).  
> *Multi-omic machine learning predictor of breast cancer therapy response.* *Nature*, 601, 623–629.  
> [https://doi.org/10.1038/s41586-021-04278-5](https://doi.org/10.1038/s41586-021-04278-5)
---

## 1. Scientific Motivation

**Research Question**  
Can we predict the response of breast cancer patients to neoadjuvant chemotherapy using multi-omics data and interpretable machine learning models?

**Motivation**  
Breast cancer treatment outcomes, especially in chemotherapy or targeted therapy, vary widely among patients. Predicting therapy response can enable tailored treatments and improve survival rates. Our project aims to develop interpretable models that integrate genomic, transcriptomic, and digital pathology features.

---

## 2. Biological Scope

- **Disease**: Breast cancer  
- **Focus**: Tumor ecosystems and therapy outcomes (pCR vs. residual disease)  
- **Tissue**: Tumor biopsies before treatment  

---

## 3. Data Overview

**Source**  
- Public GitHub repositories:
  - https://github.com/cclab-brca/neoadjuvant-therapy-response-predictor  
  - https://github.com/micrisor/NAT-ML  
- Supplementary materials from the original study  

**Summary**

| Sample Type                  | Count |
|-----------------------------|-------|
| Patients enrolled           | 180   |
| Biopsies collected          | 168   |
| WGS/WES/RNA-seq profiled    | 162–168 |
| Pathology slides available  | 166   |
| Response labeled (RCB)      | 161   |

**Data Types**
- Whole-exome sequencing (mutations)
- Copy number variation (ASCAT)
- RNA-seq expression
- Neoantigen burden (HLA typing, LOH)
- Digital pathology features (immune infiltration)

**Format**: `.tsv` tables  
**Access**: De-identified data is public. Raw data requires controlled access at EGA (EGAS00001004582).

---

## 4. Methodology

**Modeling Techniques**
- Logistic Regression  
- Random Forests  
- Support Vector Machines (SVM)  

**Interpretability**
- SHAP (Shapley Additive Explanations) for patient-level feature attribution

**Feature Engineering**
- Fisher’s Exact Test for categorical features  
- PCA for numeric data  
- Lasso for regularization and feature selection  

**Evaluation Metrics**
- AUC  
- Accuracy  
- Sensitivity  
- Specificity  

---


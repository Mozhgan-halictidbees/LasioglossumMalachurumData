# Differential Gene Expression Analysis using DESeq2

## Overview

This script performs differential gene expression (DGE) analysis using **DESeq2** on HTSeq count data. It compares two phenotypes (e.g., queen vs worker) and identifies significantly upregulated and downregulated genes based on log2 fold change and adjusted p-values.

The workflow includes:
- Importing HTSeq count files
- Building a combined gene expression matrix
- Filtering low-expression genes
- Running DESeq2 analysis
- Extracting differentially expressed genes (DEGs)
- Saving results as CSV files

---

## Requirements

### R Packages

- DESeq2

Install if needed:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

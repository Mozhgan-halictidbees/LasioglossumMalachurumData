## PCA Analysis of Fat Body and Brain Transcriptomes

## Overview

This directory  performs Principal Component Analysis (PCA) on RNA-seq gene expression data from fat body and brain tissues. The analysis uses the 500 most variable genes identified from normalized count data and applies Variance Stabilizing Transformation (VST) on all genes  prior to PCA.


The workflow generates:

* PCA plots for fat body samples
* PCA plots for brain samples
* Combined publication-quality PCA figures


---

## Required R Packages

```r
library(DESeq2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

```
---

## Input Files

### Fat Body Data

| File                          | Description                                      |
| ----------------------------- | ------------------------------------------------ |
| `raw_fatbody_mal_mozhgan.csv` | Raw count matrix (genes × samples)               |
| `ttype_fatbody.csv`           | Sample metadata containing phenotype information |

### Brain Data

| File                        | Description                                      |
| --------------------------- | ------------------------------------------------ |
| `raw_brain_mal_mozhgan.csv` | Raw count matrix (genes × samples)               |
| `ttype_brain.csv`           | Sample metadata containing phenotype information |

### Metadata Requirements

The metadata files must contain a column named:
phenotype


with sample groups such as:

queen
worker
foundress


Row names of the metadata files must match the sample names in the count matrices.

---

## Analysis Workflow

### 1. Selection of Highly Variable Genes

For each tissue:

1. Raw count data are loaded.
2. Low-expression genes are filtered:

rowSums(counts >= 10) >= 5

# We choose 5 here because our smallest group (queen) has 5 samples. For more information see Pre-filtering in https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
 

3. DESeq2 normalization is performed.
4. Variance is calculated across normalized counts.
5. The 500 genes with the highest variance are selected for PCA.

---

### 2. Principal Component Analysis

PCA is performed using:

prcomp(t(vst_counts_top500), scale. = TRUE)


where:

* Rows = samples
* Columns = genes

The first two principal components (PC1 and PC2) are plotted.

---

### 3. Sample Annotation

Here we are showing 2 egg-laying workers as A and B to differentiate them from other samples.

#### Fat Body

| Sample ID   | Label |
| ----------- | ----- |
| X109FmalW83 | A     |
| X128FmalW99 | B     |

#### Brain

| Sample ID   | Label |
| ----------- | ----- |
| X87BmalW83  | A     |
| X125BmalW99 | B     |

---

## Notes

* The PCA is based only on the 500 most variable genes.
* Gene selection is performed independently for fat body and brain datasets.
* DESeq2 normalization is applied prior to variance calculation.
* VST-transformed counts are used for PCA.

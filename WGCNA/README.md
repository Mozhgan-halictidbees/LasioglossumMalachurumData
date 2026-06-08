# WGCNA Analysis of Brain and Fat Body Transcriptomes in *Lasioglossum malachurum*

## Overview

This directory  performs weighted gene co-expression network analysis (WGCNA) on RNA-seq data from brain and fat body tissues of Lasioglossum malachurum. 

The pipeline combines DESeq2 normalization with WGCNA network construction to identify co-expressed gene modules associated with social phenotypes and biological traits.

---

## Main Steps

### 1. RNA-seq Data Processing

Raw count matrices are imported into DESeq2 and filtered to retain genes with at least 10 counts in at least 5 samples.

Data processing includes:

* DESeq2 normalization
* Variance Stabilizing Transformation (VST)
* Sample quality assessment
* Outlier detection by hierarchical clustering

---

### 2. Network Construction

For each tissue (brain and fat body):

* Soft-threshold power selection using scale-free topology criteria
* Signed co-expression network construction
* Topological Overlap Matrix (TOM) calculation
* Hierarchical clustering of genes
* Dynamic tree cutting for module detection
* Module merging based on eigengene similarity

Parameters:

| Parameter           | Value  |
| ------------------- | ------ |
| Network type        | Signed |
| Minimum module size | 30     |
| Merge cut height    | 0.25   |

---

### 3. Module–Trait Relationships

Module eigengenes are correlated with trait data using biweight midcorrelation (bicor).

Outputs include:

* Module-trait correlation matrices
* P-values
* Heatmaps of module-trait relationships

Generated figures:

* Brain module-trait heatmap
* Fat body module-trait heatmap
* Combined heatmap figure

---

### 4. Module Exploration

The workflow allows:

* Extraction of genes from individual modules
* Identification of hub genes
* Calculation of module membership (KME)
* Identification of highly connected genes (KME > 0.7)
* Detection of top hub genes within modules

---

### 5. Visualization

The pipeline generates:

* Sample clustering dendrograms
* Gene dendrograms
* Module color assignments
* Module–trait heatmaps
* Module expression heatmaps
* Combined publication-quality figures

---

## Required Input Files

### Expression Data

* `raw_brain_mal_mozhgan.csv`
* `raw_fatbody_mal_mozhgan.csv`

Rows represent genes and columns represent samples.

### Sample Metadata

* `ttype_brain.csv`
* `ttype_fatbody.csv`

### Trait Information

* `datTraits.csv`

Trait matrix used for module–trait correlation analysis.

---

## Required R Packages

```r
library(DESeq2)
library(WGCNA)
library(impute)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggrepel)
library(pheatmap)
library(grid)
library(gridExtra)
library(png)
library(sva)
```

---

## Outputs

### Network Files

* `TOM_signed.rds`
* `mal-expt1-01-dataInput_core.RData`
* `mal-expt1-NetworkConstruction-core.RData`

### Module–Trait Analysis

* `Brain_module_trait_heatmap.pdf`
* `FatBody_module_trait_heatmap.pdf`
* `Combined_Module_Trait_Heatmaps_Horizontal.pdf`

### Module Gene Lists

* `genes_<module>.csv`

### Hub Gene Results

* `hub_genes_<module>.csv`

### Module Expression Heatmaps

* `module_<module>.png`



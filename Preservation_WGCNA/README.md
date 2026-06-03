# Module Preservation Analysis

## Overview

This repository contains R scripts used to evaluate the preservation of gene co-expression modules between independent RNA-seq datasets using Weighted Gene Co-expression Network Analysis (WGCNA).

The analysis was performed separately for **fat body** and **brain** tissues. A reference co-expression network was constructed from the primary dataset and compared with an independent test dataset using WGCNA module preservation statistics.

## Workflow

### 1. Data Preprocessing

* Import raw count matrices and sample metadata.
* Filter lowly expressed genes.
* Normalize counts using DESeq2.
* Apply variance stabilizing transformation (VST).
* Inspect sample clustering and potential outliers.

### 2. Reference Network Construction

* Select a soft-thresholding power.
* Construct a signed co-expression network.
* Calculate the Topological Overlap Matrix (TOM).
* Identify gene modules using dynamic tree cutting.
* Merge highly similar modules based on eigengene correlation.
* Calculate module eigengenes and module–trait relationships.

### 3. Module Preservation Analysis

* Define the reference and test expression datasets.
* Perform module preservation analysis using WGCNA.
* Calculate preservation statistics (Zsummary) for each module.
* Assess module conservation between datasets.

### 4. Visualization

* Generate module–trait relationship heatmaps.
* Create module preservation bar plots.
* Combine brain and fat body preservation results into a single publication-quality figure.

## Required Input Files

### Reference Dataset

* Gene count matrix (`.csv`)
* Sample metadata (`.csv`)
* Trait matrix (`datTraits.csv`)
* Previously generated WGCNA objects:

  * `mal-expt1-01-dataInput_core.RData`
  * `mal-expt1-NetworkConstruction-core.RData`
  * `TOM_signed.rds`

### Test Dataset

* Gene count matrix (`.csv`)
* Sample metadata (`.csv`)

## Main Output Files

### Network Analysis

* `module_trait_heatmap.png`
* `module_trait_heatmap.pdf`
* `WGCNA_ref_test_data.RData`

### Module Preservation

* `modulePreservation_results.RData`

### Figures

* Fat body module preservation plot
* Brain module preservation plot
* `combined_plot.pdf`
* `combined_plot.png`

## Software Requirements

* R (version 4.0 or later recommended)

### R Packages

* DESeq2
* WGCNA
* ggplot2
* dplyr
* reshape2
* ggrepel
* sva
* patchwork

## Notes

* The analysis uses signed WGCNA networks.
* Preservation statistics are based on permutation testing.



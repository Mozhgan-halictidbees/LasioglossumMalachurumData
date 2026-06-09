## Overview

This repository contains the code and analysis pipelines used to investigate the molecular basis of sociality in Lasioglossum malachurum using RNA-seq data. The analyses cover the complete workflow from raw sequencing data processing to differential gene expression, gene co-expression network analysis, functional annotation, and evolutionary analyses.

The repository is organized into separate directories, each corresponding to a specific analysis. Every directory contains its own README file with detailed descriptions of the workflow, required inputs, software dependencies, and expected outputs.

## Repository Structure

### RNA_seq_Processing_Pipeline
Processing of raw RNA-seq reads, including quality control, trimming, mapping, read quantification, and generation of count matrices for downstream analyses.

### PCA_Gene_Expression_Data
Principal Component Analysis (PCA) of gene expression data is used to assess sample clustering and explore global expression patterns among samples.

### Differentially_Expressed_Gene_Analysis
Identification of differentially expressed genes between biological groups using RNA-seq count data.

### Functional_Gene_Annotation
Functional characterization of genes through Gene Ontology (GO), KEGG pathway annotation, description, and potential gene code.

###Functional_Enrichment_Analysis
Functional enrichment analysis was performed on differentially expressed genes (DEGs) using topGO to identify significantly enriched Gene Ontology (GO) terms relative to a background gene set.

### WGCNA
Weighted Gene Co-expression Network Analysis (WGCNA) is used to identify modules of co-expressed genes and evaluate their associations with biological traits of interest.

### Preservation_WGCNA
Module preservation analysis is used to assess whether co-expression modules identified in one dataset are preserved in independent datasets.

### Phylostratigraphy
Evolutionary age assignment of genes using phylostratigraphic approaches. This workflow includes taxonomic lineage reconstruction, protein homology searches, gene age assignment, overlap analyses with differentially expressed genes, and statistical testing of age distributions.

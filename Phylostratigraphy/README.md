# Phylostratigraphy Analysis

## Overview

This directory contains the scripts used to perform a phylostratigraphy analysis of L. malachurum. The workflow assigns genes to evolutionary strata based on sequence homology and taxonomic lineage information, followed by downstream analyses examining the distribution of differentially expressed genes (DEGs) across evolutionary age classes.

The focal species used for downloading the database is *Apis mellifera*.

---

## Workflow

### 1. Download Reference Proteomes

**Script:** `01_Download.R`

Downloads protein datasets (`.faa`) for the focal species and related taxa required for phylostratigraphy analyses.

### 2. Retrieve Taxonomic Lineages

**Script:** `02_taxid_to_lineage.R`

Obtains complete taxonomic lineages for all downloaded species based on NCBI taxonomy identifiers.

### 3. Assign Species to Evolutionary Strata

**Script:** `03_species_stratum_assignment.R`

Assigns species to taxonomic strata based on lineage information.

Example output:

| Stratum | Number of Species |
|----------|----------|
| Cellular organisms | 57 |
| Eukaryota | 19 |
| Bilateria | 40 |
| Insecta | 18 |
| Hymenoptera | 13 |
| Apoidea | 9 |
| Halictidae | 1 |

---

## Expanding Halictidae Sampling

The default database contained only one halictid bee (*Dufourea novaeangliae*). To improve taxonomic representation, additional halictid species were included. The *Dufourea novaeangliae* proteome was replaced with a more complete protein set.

Additional species:

| Species | TaxID |
|----------|----------|
| Nomia melanderi | 2448451 |
| Agapostemon virescens | 115084 |
| Dufourea novaeangliae | 178035 |
| Halictus quadricinctus | 115107 |
| Lasioglossum albipes | 88501 |
| Lasioglossum leucozonium | 88532 |
| Lasioglossum figueresi | 160208 |
| Lasioglossum pauxillum | 88516 |
| Lasioglossum zephyrum | 88500 |
| Megalopta genalis | 115081 |

### 4. Download and Process Additional Proteomes

**Script:** `04_add_related_stratum.sh`

Downloads proteomes for the additional halictid species, extracts the longest isoform for each gene, and saves the resulting protein datasets as `.faa` files named according to their TaxID.

### 5. Similarity Searches

**Script:** `05_blast.R`

Performs sequence similarity searches required for gene age assignment.

### 6. Retrieve Taxonomic Lineages (Updated)

**Script:** `06_taxid_to_lineage.R`

Updates lineage information after inclusion of the additional species.

### 7. Assign Species to Evolutionary Strata (Updated)

**Script:** `07_species_stratum_assignment.R`

Updates species-stratum assignments after database expansion.

### 8. Gene Age Assignment

**Script:** `08_phylostratigraphy_gene_age_assignment.R`

Assigns each gene to its oldest detectable evolutionary stratum.

### 9. Export Genes by Stratum

**Script:** `09_save_genes_by_stratum.R`

Creates separate gene lists for each evolutionary stratum.

### 10. DEG Overlap Analysis

**Script:** `10_common_DEGs.R`

Determines the overlap between genes assigned to each evolutionary stratum and:

- Brain DEGs
- Fat body DEGs

### 11. Visualization

**Script:** `11_plotting.R`

Generates figures summarizing gene age distributions and DEG enrichment patterns.

### 12. Statistical Testing

**Script:** `12_chisquare_test.R`

Performs chi-square tests to evaluate whether DEG distributions differ significantly among evolutionary strata.

---

## Required R Packages

```r
library(DESeq2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
```

## External Software

- BLAST+
- Bash
- R (version 4.0 or later)

---

## Outputs

The workflow produces:

- Species-to-stratum assignments
- Gene age assignments
- Gene lists for each phylostratum
- DEG overlap summaries
- Chi-square test results

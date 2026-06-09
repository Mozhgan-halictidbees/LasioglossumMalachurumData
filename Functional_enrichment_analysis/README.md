# GO Enrichment Analysis with topGO

## Description

This script performs Gene Ontology (GO) enrichment analysis using the R package topGO. It identifies GO terms that are significantly overrepresented among differentially expressed genes (DEGs) compared to a background gene set.

The analysis uses Fisher's exact test with the "weight01" and "parentchild" algorithms implemented in topGO.

## Requirements

R packages:

* topGO
* ggplot2

## Input Files

1. DEG file
   A text file containing one differentially expressed gene (DEG) ID per line.

2. Background gene file
   A text file containing all genes used as the reference/background set.

3. Gene-to-GO mapping file
   A custom mapping file linking gene IDs to GO terms.

## Parameters

ontology:
GO ontology to analyze:

* BP = Biological Process
* MF = Molecular Function
* CC = Cellular Component

pval_threshold:
Significance threshold for enriched GO terms (default: 0.05).

## Output Files

GO_enrichment_results.csv
Complete GO enrichment results table. The table shows enriched GO Biological Process terms with their identifiers (GO.ID), full descriptions (Full Term), number of associated significant genes (Significant), and Fisher exact test p-values (Fisher). 

GO_enrichment_barplot.pdf
Barplot of significant GO terms.

GO_enrichment_barplot_top10.pdf
Barplot showing the most significant enriched GO terms.


## Workflow

1. Load DEG and background gene lists.
2. Create a binary gene list for enrichment testing.
3. Load gene-to-GO annotations.
4. Construct a topGO data object.
5. Run enrichment analysis using Fisher's exact test.
6. Extract and save enriched GO terms.
7. Generate barplots of significant GO terms.

## Notes

The script currently analyzes the Molecular Function (MF) ontology. To analyze Biological Process (BP) or Cellular Component (CC), modify the ontology parameter accordingly.


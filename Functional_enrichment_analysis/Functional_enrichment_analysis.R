# Load required libraries
library(topGO)
library(ggplot2)

# Input files and parameters
deg_file <- "C:/Users/....../..../DEG_queen_worker_mal_mozhgan_fatbody.txt"  # File containing a list of DEGs
bg_file <- "C:/Users/....../Csec_universe_reference_topGO.txt"  # File containing the background gene list
gene2go_file <- "C:/Users...../gene2go_bees.txt"  # Custom Gene-to-GO mapping file
output_folder <- "C:/Users/......./topGO2"  # Directory to save results
pval_threshold <- 0.05  # Significance threshold for p-value
ontology <- "MF"  # Change to "MF" or "CC" if you want to analyze other ontologies

# Load gene lists
deg_genes <- scan(deg_file, what = "", quiet = TRUE)
bg_genes <- scan(bg_file, what = "", quiet = TRUE)

# Create geneList: 1 for DEGs, 0 for background genes
all_genes <- unique(c(bg_genes, deg_genes))
gene_list <- factor(as.integer(all_genes %in% deg_genes))
names(gene_list) <- all_genes

# Load the gene-to-GO mapping
gene2GO <- readMappings(file = gene2go_file)

# Create topGOdata object
GOdata <- new("topGOdata",
              ontology = ontology,
              allGenes = gene_list,
              geneSelectionFun = function(x) x == 1,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)

# Run enrichment analysis
result_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "Fisher")
result_parentchild <- runTest(GOdata, algorithm = "parentchild", statistic = "Fisher")

# Get significant GO terms
significant_terms <- GenTable(GOdata,
                              weight01 = result_weight01,
                              parentchild = result_parentchild,
                              orderBy = "weight01",
                              topNodes = length(usedGO(GOdata)))

# Save results
write.csv(significant_terms, file = file.path(output_folder, "GO_enrichment_results.csv"), row.names = FALSE)

# Generate bar chart of significant GO terms
top_terms <- significant_terms[significant_terms$weight01 < pval_threshold, ]
if (nrow(top_terms) > 0) {
  ggplot(top_terms, aes(x = reorder(Term, -weight01), y = -log10(weight01))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "GO Enrichment Analysis", x = "GO Term", y = "-log10(p-value)") +
    theme_minimal() +
    ggsave(file.path(output_folder, "GO_enrichment_barplot.pdf"))
} else {
  message("No significant GO terms found.")
}



str(significant_terms)
significant_terms$weight01 <- as.numeric(significant_terms$weight01)

top_terms <- significant_terms[significant_terms$weight01 < pval_threshold, ]
top_terms <- significant_terms[significant_terms$weight01 < pval_threshold, ]

if (nrow(top_terms) > 0) {
  # Create the plot object
  p <- ggplot(top_terms, aes(x = reorder(Term, -weight01), y = -log10(weight01))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "GO Enrichment Analysis", x = "GO Term", y = "-log10(p-value)") +
    theme_minimal()
  
  # Save the plot using ggsave
  ggsave(file.path(output_folder, "GO_enrichment_barplot.pdf"), plot = p, width = 8, height = 6)
} else {
  message("No significant GO terms found.")
}


top_terms <- significant_terms[significant_terms$weight01 < pval_threshold, ]

if (nrow(top_terms) > 0) {
  # Create the plot object
  p <- ggplot(top_terms, aes(x = reorder(Term, -weight01), y = -log10(weight01))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "GO Enrichment Analysis", x = "GO Term", y = "-log10(p-value)") +
    theme_minimal()
  
  # Print the plot to view it in the R session
  print(p)
} else {
  message("No significant GO terms found.")
}


# Limit the number of top terms to display
top_n <- 25
top_terms <- significant_terms[significant_terms$weight01 < pval_threshold, ]
top_terms <- head(top_terms[order(top_terms$weight01), ], top_n)

if (nrow(top_terms) > 0) {
  # Create the plot object
  p <- ggplot(top_terms, aes(x = reorder(Term, -weight01), y = -log10(weight01))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Top GO Enrichment Terms", x = "GO Term", y = "-log10(p-value)") +
    theme_minimal()
  
  # Print the plot to view it in the R session
  print(p)
  
  # Save the plot to a file
  ggsave(file.path(output_folder, "GO_enrichment_barplot_top10.pdf"), plot = p, width = 8, height = 6)
} else {
  message("No significant GO terms found.")
}



# Generate a file for REVIGO input
revigo_file <- file.path(output_folder, "go_ids_pvalues.csv")
write.table(significant_terms[, c("GO.ID", "weight01")],
            file = revigo_file,
            sep = ",", row.names = FALSE, col.names = c("GO.ID", "p-value"))

message("GO enrichment analysis completed. Results saved to: ", output_folder)

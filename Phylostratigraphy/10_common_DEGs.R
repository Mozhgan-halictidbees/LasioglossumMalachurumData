# =========================================================
# Count overlap between phylostratigraphy gene lists
# and DEG lists (brain and fatbody)
# =========================================================

# Set working directory
setwd("C:/Users/...../Desktop/new_phylostratigraphy/stratum_gene_lists")

# ---------------------------------------------------------
# Files
# ---------------------------------------------------------

stratum_files <- c(
  "Cellular_organisms_genes.csv",
  "Eukaryota_genes.csv",
  "Bilateria_genes.csv",
  "Insecta_genes.csv",
  "Hymenoptera_genes.csv",
  "Apoidea_genes.csv",
  "Halictidae_genes.csv",
  "Lasioglossum_genes.csv"
)

# DEG files
brain_file <- "DEG_brain.csv"
fatbody_file <- "DEG_fatbody.csv"

# ---------------------------------------------------------
# Function to clean gene IDs
# Removes -RA, -RB, -RC, etc.
# ---------------------------------------------------------

clean_genes <- function(x) {
  x <- trimws(x)
  x <- sub("-R[A-Z]+$", "", x)
  return(unique(x))
}

# ---------------------------------------------------------
# Read DEG files
# (files have NO column names)
# ---------------------------------------------------------

deg_brain <- read.csv(brain_file, header = FALSE, stringsAsFactors = FALSE)
deg_fatbody <- read.csv(fatbody_file, header = FALSE, stringsAsFactors = FALSE)

deg_brain_genes <- clean_genes(deg_brain$V1)
deg_fatbody_genes <- clean_genes(deg_fatbody$V1)

# ---------------------------------------------------------
# Create results dataframe
# ---------------------------------------------------------

results <- data.frame(
  Group = character(),
  Total_Genes = numeric(),
  Common_with_Brain = numeric(),
  Common_with_Fatbody = numeric(),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------
# Process each stratum file
# ---------------------------------------------------------

for (file in stratum_files) {
  
  # Read file
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Extract gene column
  genes <- clean_genes(df$gene)
  
  # Count overlaps
  brain_overlap <- length(intersect(genes, deg_brain_genes))
  fatbody_overlap <- length(intersect(genes, deg_fatbody_genes))
  
  # Add results
  results <- rbind(
    results,
    data.frame(
      Group = gsub("_genes.csv", "", file),
      Total_Genes = length(genes),
      Common_with_Brain = brain_overlap,
      Common_with_Fatbody = fatbody_overlap,
      stringsAsFactors = FALSE
    )
  )
}

# ---------------------------------------------------------
# Print results
# ---------------------------------------------------------

print(results)

# ---------------------------------------------------------
# Save results
# ---------------------------------------------------------

write.csv(results,
          "DEG_overlap_summary.csv",
          row.names = FALSE)

cat("\nAnalysis completed.\n")
cat("Results saved as: DEG_overlap_summary.csv\n")

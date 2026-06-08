# ============================================
# LOAD LIBRARIES
# ============================================

library(readr)
library(dplyr)

# ============================================
# INPUT FILE
# ============================================

gene_age_file <- "/media/..../new_phylostratigraphy/phylostratigraphy_gene_age_assignment/gene_ages.csv"

# ============================================
# READ GENE AGE TABLE
# ============================================

gene_age_final <- read_csv(gene_age_file)

# ============================================
# CREATE OUTPUT DIRECTORY
# ============================================

output_dir <- "/media/...../new_phylostratigraphy/phylostratigraphy_gene_age_assignment/stratum_gene_lists"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# CHECK COLUMN NAMES
# ============================================

print(colnames(gene_age_final))

# ============================================
# SAVE GENE LISTS FOR EACH STRATUM
# ============================================

strata_levels <- unique(gene_age_final$assigned_stratum)

for (s in strata_levels) {
  
  genes <- gene_age_final %>%
    filter(assigned_stratum == s) %>%
    select(gene)
  
  # Replace spaces and special characters
  clean_name <- gsub("[^A-Za-z0-9_]", "_", s)
  
  outname <- paste0(clean_name, "_genes.csv")
  
  write_csv(
    genes,
    file.path(output_dir, outname)
  )
  
  cat("Saved:", outname, "\n")
}

cat("\nAll gene lists saved successfully.\n")


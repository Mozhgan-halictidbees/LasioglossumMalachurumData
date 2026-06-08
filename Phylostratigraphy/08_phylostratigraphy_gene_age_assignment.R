
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# =========================================================
# PATHS
# =========================================================

blast_dir <- "/media/...../new_phylostratigraphy/faa"

taxonomy_file <- "/media/.../new_phylostratigraphy/taxid_to_lineage/taxonomy_lineages.csv"

output_dir <- "/media/...../new_phylostratigraphy/phylostratigraphy_gene_age_assignment"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# STRATA ORDER
# oldest --> youngest
# =========================================================

strata_order <- c(
  "Cellular organisms",
  "Eukaryota",
  "Bilateria",
  "Insecta",
  "Hymenoptera",
  "Apoidea",
  "Halictidae",
  "Lasioglossum"
)

# =========================================================
# LOAD TAXONOMY
# =========================================================

taxonomy <- read_csv(
  taxonomy_file,
  col_names = c("taxid", "rank", "taxon"),
  show_col_types = FALSE
)

# =========================================================
# MAP LINEAGES TO CUSTOM STRATA
# =========================================================

taxonomy <- taxonomy %>%
  mutate(
    stratum = case_when(
      
      taxon == "Lasioglossum" ~ "Lasioglossum",
      
      taxon == "Halictidae" ~ "Halictidae",
      
      taxon == "Apoidea" ~ "Apoidea",
      
      taxon == "Hymenoptera" ~ "Hymenoptera",
      
      taxon == "Insecta" ~ "Insecta",
      
      taxon == "Bilateria" ~ "Bilateria",
      
      taxon == "Eukaryota" ~ "Eukaryota",
      
      taxon == "cellular organisms" ~ "Cellular organisms",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stratum))

# =========================================================
# ORDER STRATA
# =========================================================

taxonomy$stratum <- factor(
  taxonomy$stratum,
  levels = strata_order,
  ordered = TRUE
)

# =========================================================
# FOR EACH SPECIES:
# KEEP MOST SPECIFIC STRATUM
# =========================================================

species_strata <- taxonomy %>%
  arrange(desc(stratum)) %>%
  group_by(taxid) %>%
  slice(1) %>%
  ungroup()

# =========================================================
# LOAD BLAST FILES
# =========================================================

blast_files <- list.files(
  blast_dir,
  pattern = "\\.tab$",
  full.names = TRUE
)

read_blast <- function(f) {
  
  species_taxid <- str_remove(basename(f), ".tab")
  
  df <- read_tsv(
    f,
    col_names = c(
      "qseqid",
      "sseqid",
      "qstart",
      "qend",
      "sstart",
      "send",
      "evalue",
      "bitscore",
      "staxids"
    ),
    show_col_types = FALSE
  )
  
  df$species_taxid <- species_taxid
  
  return(df)
}

blast_all <- map_dfr(blast_files, read_blast)

# =========================================================
# FILTER BLAST
# =========================================================

blast_filtered <- blast_all %>%
  filter(evalue < 1e-5)

# =========================================================
# MERGE SPECIES STRATA
# =========================================================

blast_tax <- blast_filtered %>%
  left_join(
    species_strata,
    by = c("species_taxid" = "taxid")
  ) %>%
  filter(!is.na(stratum))

# =========================================================
# FOR EACH GENE:
# CHOOSE OLDEST DETECTABLE STRATUM
# =========================================================

gene_age <- blast_tax %>%
  arrange(stratum, evalue, desc(bitscore)) %>%
  group_by(qseqid) %>%
  slice(1) %>%
  ungroup()

# =========================================================
# FINAL TABLE
# =========================================================

gene_age_final <- gene_age %>%
  select(
    gene = qseqid,
    assigned_stratum = stratum,
    evalue,
    bitscore,
    species_taxid
  )

# =========================================================
# SAVE
# =========================================================

write_csv(
  gene_age_final,
  file.path(output_dir, "gene_ages.csv")
)

summary_table <- gene_age_final %>%
  count(assigned_stratum)

write_csv(
  summary_table,
  file.path(output_dir, "stratum_counts.csv")
)

# =========================================================
# PRINT
# =========================================================

cat("\n====================================\n")
cat("DONE\n")
cat("====================================\n")

print(summary_table)

cat("\nResults saved to:\n")
cat(output_dir, "\n")


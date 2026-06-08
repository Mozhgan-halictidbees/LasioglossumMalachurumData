
library(dplyr)
library(readr)

# ==========================================
# LOAD TAXONOMY
# ==========================================

taxonomy <- read_csv(
  "/media/mozhgan/..../taxid_to_lineage/taxonomy_lineages.csv",
  col_names = c("taxid", "rank", "taxon"),
  show_col_types = FALSE
)

# ==========================================
# DEFINE STRATA
# oldest --> youngest
# ==========================================

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

# ==========================================
# MAP TAXONOMY TO STRATA
# ==========================================

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

# ==========================================
# ORDER STRATA
# ==========================================

taxonomy$stratum <- factor(
  taxonomy$stratum,
  levels = strata_order,
  ordered = TRUE
)

# ==========================================
# KEEP MOST SPECIFIC STRATUM
# FOR EACH SPECIES
# ==========================================

species_assignment <- taxonomy %>%
  arrange(desc(stratum)) %>%
  group_by(taxid) %>%
  slice(1) %>%
  ungroup()

# ==========================================
# COUNT SPECIES
# ==========================================

species_counts <- species_assignment %>%
  count(stratum)

print(species_counts)

# ==========================================
# SAVE RESULTS
# ==========================================

write_csv(
  species_assignment,
  "/media/mozhgan/....../species_stratum_assignment/species_assignment.csv"
)

write_csv(
  species_counts,
  "/media/mozhgan/....../species_stratum_assignment/species_counts.csv"
)

cat("\nDONE\n")


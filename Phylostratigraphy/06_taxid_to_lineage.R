
library(taxizedb)
library(dplyr)
library(ape)

cat("\n=== PREPARING TAXONOMY TREE ===\n")

# --------------------------------------------------
# WORK & OUTPUT DIRECTORY
# --------------------------------------------------

data_dir <- "/media/...../new_phylostratigraphy/faa"
out_dir  <- "/media/....../new_phylostratigraphy/taxid_to_lineage"

# create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

setwd(data_dir)

# --------------------------------------------------
# GET TAXA FROM FAA FILES
# --------------------------------------------------

faa_files <- list.files(pattern="\\.faa$")

taxids <- gsub("\\.faa$", "", faa_files)

# keep only numeric taxids
taxids <- taxids[grepl("^[0-9]+$", taxids)]

cat("Taxa found:", length(taxids), "\n")

# --------------------------------------------------
# RETRIEVE TAXONOMY
# --------------------------------------------------

cat("\nRetrieving taxonomy from NCBI...\n")

lineages <- classification(
  taxids,
  db = "ncbi"
)

# --------------------------------------------------
# SAVE RAW TAXONOMY
# --------------------------------------------------

saveRDS(
  lineages,
  file.path(out_dir, "taxonomy_lineages.rds")
)

# --------------------------------------------------
# CONVERT TO TABLE
# --------------------------------------------------

cat("\nConverting lineage tables...\n")

lineage_df <- bind_rows(lapply(names(lineages), function(x){
  
  lin <- lineages[[x]]
  
  if(is.null(lin)) return(NULL)
  
  data.frame(
    staxid = x,
    rank = lin$rank,
    taxon = lin$name,
    stringsAsFactors = FALSE
  )
  
}))

# --------------------------------------------------
# SAVE TABLE
# --------------------------------------------------

write.csv(
  lineage_df,
  file.path(out_dir, "taxonomy_lineages.csv"),
  row.names = FALSE
)

# --------------------------------------------------
# SUMMARY
# --------------------------------------------------

cat("\n=== DONE ===\n")

cat("Output directory:", out_dir, "\n")
cat("Taxa processed:", length(unique(lineage_df$staxid)), "\n")

print(head(lineage_df))



library(phylostratr)
library(dplyr)
library(ape)

cat("\n=== DOWNLOAD ONLY PIPELINE (HONEY BEE FOCAL) ===\n")

# ==================================================
# WORK DIRECTORY
# ==================================================

data_dir <- "/media/...."

setwd(data_dir)

# ==================================================
# CLEAN OLD CACHE
# ==================================================

cat("\n[0] Cleaning old files...\n")

unlink("cache", recursive = TRUE, force = TRUE)
unlink("faa", recursive = TRUE, force = TRUE)

dir.create("cache", showWarnings = FALSE)
dir.create("faa", showWarnings = FALSE)

# ==================================================
# STEP 1: BUILD INITIAL TREE
# ==================================================

cat("\n[1] Building tree...\n")

# 7460 = Apis mellifera
strata <- uniprot_strata(
  "7460",
  from = 2
)

cat("Initial taxa:",
    length(strata@tree$tip.label), "\n")

# ==================================================
# STEP 2: PRUNE EUKARYOTIC TREE
# ==================================================

cat("\n[2] Pruning eukaryotic tree...\n")

strata <- strata %>%
  strata_apply(
    f = diverse_subtree,
    n = 5,
    weights = uniprot_weight_by_ref()
  )

cat("Taxa after pruning:",
    length(strata@tree$tip.label), "\n")

# ==================================================
# STEP 3: ADD RECOMMENDED PROKARYOTES
# ==================================================

cat("\n[3] Adding recommended prokaryotes...\n")

strata <- use_recommended_prokaryotes(strata)

cat("Taxa after adding prokaryotes:",
    length(strata@tree$tip.label), "\n")

# ==================================================
# STEP 4: ADD IMPORTANT MANUAL TAXA
# ==================================================

cat("\n[4] Adding manual taxa...\n")

strata <- add_taxa(
  strata,
  c(
    "7460",   # Apis mellifera
    "7425",   # Bombus terrestris
    "7227",   # Drosophila melanogaster
    "9606",   # Human
    "3702"    # Arabidopsis
  )
)

cat("Taxa after manual additions:",
    length(strata@tree$tip.label), "\n")

# ==================================================
# STEP 5: DOWNLOAD PROTEOMES SAFELY
# ==================================================

cat("\n[5] Downloading UniProt proteomes...\n")

all_taxa <- strata@tree$tip.label

successful_taxa <- c()
failed_taxa <- c()

for(taxid in all_taxa){
  
  cat("\n----------------------------------\n")
  cat("Downloading taxid:", taxid, "\n")
  
  temp_strata <- strata
  
  temp_strata@tree <- keep.tip(
    temp_strata@tree,
    taxid
  )
  
  tryCatch(
    
    {
      
      temp_strata <- uniprot_fill_strata(
        temp_strata,
        dir = file.path(data_dir, "faa")
      )
      
      faa_file <- file.path(
        data_dir,
        "faa",
        paste0(taxid, ".faa")
      )
      
      if(file.exists(faa_file)){
        
        successful_taxa <- c(
          successful_taxa,
          taxid
        )
        
        cat("SUCCESS:", taxid, "\n")
        
      } else {
        
        failed_taxa <- c(
          failed_taxa,
          taxid
        )
        
        cat("FAILED:", taxid, "\n")
        cat("Reason: faa file not created\n")
      }
      
    },
    
    error = function(e){
      
      failed_taxa <<- c(
        failed_taxa,
        taxid
      )
      
      cat("FAILED:", taxid, "\n")
      cat("Reason:", e$message, "\n")
    }
  )
}

# ==================================================
# STEP 6: SUMMARY
# ==================================================

cat("\n====================================\n")
cat("DOWNLOAD COMPLETE\n")
cat("====================================\n")

cat("\nSuccessful downloads:",
    length(successful_taxa), "\n")

cat("Failed downloads:",
    length(failed_taxa), "\n")

# --------------------------------------------------
# SAVE SUCCESSFUL TAXA
# --------------------------------------------------

write.table(
  successful_taxa,
  file = "successful_taxa.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# --------------------------------------------------
# SAVE FAILED TAXA
# --------------------------------------------------

write.table(
  failed_taxa,
  file = "failed_taxa.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# --------------------------------------------------
# PRINT FAILED TAXA
# --------------------------------------------------

cat("\nFAILED TAXA:\n")

if(length(failed_taxa) > 0){
  
  print(failed_taxa)
  
} else {
  
  cat("None\n")
}

# --------------------------------------------------
# COUNT FAA FILES
# --------------------------------------------------

faa_files <- list.files(
  "faa",
  pattern = "\\.faa$"
)

cat("\nFAA files available:",
    length(faa_files), "\n")

# ==================================================
# SAVE FINAL OBJECT
# ==================================================

saveRDS(
  strata,
  "strata_downloaded_only.rds"
)

cat("\nSaved: strata_downloaded_only.rds\n")


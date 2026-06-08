
library(dplyr)

cat("\n=== 88512 PHYLOSTRATIGRAPHY BLAST ===\n")

# --------------------------------------------------
# WORK DIRECTORY
# --------------------------------------------------

data_dir <- "/media/....../new_phylostratigraphy/faa"

setwd(data_dir)

# --------------------------------------------------
# QUERY FASTA
# --------------------------------------------------

query_fasta <- file.path(
  data_dir,
  "88512.faa"
)

if(!file.exists(query_fasta)){
  stop("88512.faa not found")
}

cat("\nUsing query:\n")
cat(query_fasta, "\n")

# --------------------------------------------------
# FIND TARGET TAXA FASTAS
# --------------------------------------------------

target_fastas <- list.files(
  data_dir,
  pattern="\\.faa$",
  full.names=TRUE
)

# remove self blast
target_fastas <- target_fastas[
  basename(target_fastas) != "88512.faa"
]

cat("\nTarget taxa detected:", length(target_fastas), "\n")

# --------------------------------------------------
# MAKE BLAST DATABASE DIRECTORY
# --------------------------------------------------

blastdb_dir <- file.path(data_dir, "blastdb")

dir.create(blastdb_dir, showWarnings = FALSE)

# --------------------------------------------------
# LOG FILE
# --------------------------------------------------

log_file <- file.path(data_dir, "blast.log")

cat(
  paste0(
    "\n=== BLAST RUN STARTED: ",
    Sys.time(),
    " ===\n"
  ),
  file = log_file,
  append = TRUE
)

# --------------------------------------------------
# RUN BLAST
# --------------------------------------------------

for(faa in target_fastas){
  
  taxid <- basename(faa)
  taxid <- sub("\\.faa$", "", taxid)
  
  cat("\n====================================\n")
  cat("Processing taxon:", taxid, "\n")
  cat("====================================\n")
  
  write(
    paste0("Processing ", taxid),
    file = log_file,
    append = TRUE
  )
  
  # --------------------------------------------------
  # BUILD DATABASE
  # --------------------------------------------------
  
  db_path <- file.path(
    blastdb_dir,
    taxid
  )
  
  cmd_db <- paste(
    "makeblastdb",
    "-in", shQuote(faa),
    "-dbtype prot",
    "-out", shQuote(db_path)
  )
  
  cat("\nBuilding database...\n")
  
  db_status <- system(cmd_db)
  
  if(db_status != 0){
    
    cat("ERROR building database for", taxid, "\n")
    
    write(
      paste0("ERROR building database for ", taxid),
      file = log_file,
      append = TRUE
    )
    
    next
  }
  
  # --------------------------------------------------
  # OUTPUT FILE
  # --------------------------------------------------
  
  out_file <- file.path(
    data_dir,
    paste0(taxid, ".tab")
  )
  
  # --------------------------------------------------
  # BLASTP
  # --------------------------------------------------
  
  cmd_blast <- paste(
    "blastp",
    "-query", shQuote(query_fasta),
    "-db", shQuote(db_path),
    "-out", shQuote(out_file),
    "-evalue 1e-5",
    "-outfmt '6 qseqid sseqid qstart qend sstart send evalue bitscore'",
    "-num_threads 8",
    "-seg no"
  )
  
  cat("Running BLASTP...\n")
  
  blast_status <- system(cmd_blast)
  
  if(blast_status != 0){
    
    cat("ERROR running blast for", taxid, "\n")
    
    write(
      paste0("ERROR running blast for ", taxid),
      file = log_file,
      append = TRUE
    )
    
    next
  }
  
  cat("Finished:", taxid, "\n")
  
  write(
    paste0("Finished ", taxid),
    file = log_file,
    append = TRUE
  )
  
}

cat("\n=== ALL BLASTS FINISHED ===\n")

write(
  paste0(
    "\n=== BLAST RUN FINISHED: ",
    Sys.time(),
    " ===\n"
  ),
  file = log_file,
  append = TRUE
  

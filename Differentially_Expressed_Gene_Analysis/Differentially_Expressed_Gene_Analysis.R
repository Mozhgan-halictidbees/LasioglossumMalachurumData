library(DESeq2)

# =====================================
# Read HTSeq count files
# =====================================

count_files <- dir(
  "C:/Users/",
  full.names = TRUE
)

count_list <- lapply(count_files, function(file) {
  read.table(file)[, 2]
})

combined_counts <- do.call(cbind, count_list)

# Gene names
rownames(combined_counts) <- read.table(count_files[1])[, 1]

# Sample names
colnames(combined_counts) <- gsub(
  ".count",
  "",
  dir("C:/Users/")
)

# Remove HTSeq summary rows
combined_counts <- combined_counts[
  !grepl("^__", rownames(combined_counts)),
]

# =====================================
# Read metadata
# =====================================

ttype <- read.csv(
  "C:/Users/",
  row.names = 1
)

ttype$phenotype <- factor(ttype$phenotype)

# =====================================
# Create DESeq2 object
# =====================================

dds <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData = ttype,
  design = ~ phenotype
)

# Keep genes with at least 10 counts in at least 5 samples
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Optional VST transformation
dds_vst <- vst(dds, blind = FALSE)

# =====================================
# Differential expression
# Queen vs Worker
# =====================================

res <- results(
  dds,
  contrast = c("phenotype", "queen", "worker")
)

res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df <- na.omit(res_df)

# =====================================
# Extract DEGs
# =====================================

upregulated <- res_df[
  res_df$log2FoldChange > 1 &
    res_df$padj < 0.05,
]

downregulated <- res_df[
  res_df$log2FoldChange < -1 &
    res_df$padj < 0.05,
]

# =====================================
# Save results
# =====================================

write.csv(
  upregulated,
  "C:/Users/",
  row.names = FALSE
)

write.csv(
  downregulated,
  "C:/Users/",
  row.names = FALSE
)

cat("Analysis completed.\n")
cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")

#########################################
# PCA Analysis for Fat Body and Brain

#########################################

library(DESeq2)
library(ggplot2)
library(dplyr)
library(patchwork)

# ======================================
# === 1. FAT BODY
# ======================================

# ---------------------------
# Step 1: Select top 500 variable genes 
# ---------------------------
setwd("C:/Users/")  #path to the excluding egg-laying workers dataset

combined_counts_fatbody <- read.csv("./raw_fatbody_mal_mozhgan.csv", header = TRUE, row.names = 1)
phenotype_data_fatbody <- read.csv("./ttype_fatbody.csv", header = TRUE, row.names = 1)

dds_fatbody <- DESeqDataSetFromMatrix(countData = combined_counts_fatbody,
                                      colData = phenotype_data_fatbody,
                                      design = ~ phenotype)

keep <- rowSums(counts(dds_fatbody) >= 10) >= 5
dds_fatbody <- dds_fatbody[keep,]

dds_fatbody <- DESeq(dds_fatbody)
norm_counts_fatbody <- counts(dds_fatbody, normalized = TRUE)
gene_variances_fatbody <- apply(norm_counts_fatbody, 1, var)
top500_genes_fatbody <- names(sort(gene_variances_fatbody, decreasing = TRUE))[1:500]
cat("Top 500 genes selected for PCA (Fat body)\n")

# ---------------------------
# Step 2: Do VST + PCA 
# ---------------------------
setwd("C:/Users/") # path to the complete dataset (including the egg-laying workers)


combined_counts_fatbody <- read.csv("./raw_fatbody_mal_mozhgan.csv", header = TRUE, row.names = 1)
phenotype_data_fatbody <- read.csv("./ttype_fatbody.csv", header = TRUE, row.names = 1)

dds_fatbody <- DESeqDataSetFromMatrix(countData = combined_counts_fatbody,
                                      colData = phenotype_data_fatbody,
                                      design = ~ phenotype)

keep <- rowSums(counts(dds_fatbody) >= 10) >= 5
dds_fatbody <- dds_fatbody[keep,]

dds_fatbody <- DESeq(dds_fatbody)

dds_vst_fatbody <- vst(dds_fatbody, blind = FALSE)
vst_counts_fatbody <- assay(dds_vst_fatbody)
vst_counts_top500_fatbody <- vst_counts_fatbody[top500_genes_fatbody, ]

pca_vst_fatbody <- prcomp(t(vst_counts_top500_fatbody), scale. = TRUE)

pca_vst_df_fatbody <- data.frame(pca_vst_fatbody$x,
                                 Phenotype = phenotype_data_fatbody$phenotype,
                                 Sample = rownames(phenotype_data_fatbody))

pca_vst_df_fatbody$Label <- ifelse(pca_vst_df_fatbody$Sample == "X109FmalW83", "A",
                                   ifelse(pca_vst_df_fatbody$Sample == "X128FmalW99", "B", ""))

p_pca_vst_fatbody <- ggplot(pca_vst_df_fatbody, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), vjust = -1.2, size = 5, color = "black", fontface = "bold") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("queen" = "dodgerblue",
                                "worker" = "#ff1493",
                                "foundress" = "#969696")) +
  labs(title = "Fat body transcriptome",
       x = paste0("PC1 (", round(summary(pca_vst_fatbody)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_vst_fatbody)$importance[2,2] * 100, 1), "%)"))
print(p_pca_vst_fatbody)

ggsave("PCA_top500_vst_fatbody.jpg", plot = p_pca_vst_fatbody, width = 7, height = 5, dpi = 300)
ggsave("PCA_top500_vst_fatbody.pdf", plot = p_pca_vst_fatbody, width = 7, height = 5)


# ======================================
# === 2. BRAIN
# ======================================

# ---------------------------
# Step 1: Select top 500 variable genes 
# ---------------------------
setwd("C:/Users/") #path to the excluding egg-laying workers dataset

combined_counts_brain <- read.csv("./raw_brain_mal_mozhgan.csv", header = TRUE, row.names = 1)
phenotype_data_brain <- read.csv("./ttype_brain.csv", header = TRUE, row.names = 1)

dds_brain <- DESeqDataSetFromMatrix(countData = combined_counts_brain,
                                    colData = phenotype_data_brain,
                                    design = ~ phenotype)

keep <- rowSums(counts(dds_brain) >= 10) >= 5
dds_brain <- dds_brain[keep,]

dds_brain <- DESeq(dds_brain)
norm_counts_brain <- counts(dds_brain, normalized = TRUE)
gene_variances_brain <- apply(norm_counts_brain, 1, var)
top500_genes_brain <- names(sort(gene_variances_brain, decreasing = TRUE))[1:500]
cat("Top 500 genes selected for PCA (Brain)\n")

# ---------------------------
# Step 2: Do VST + PCA 
# ---------------------------
setwd("C:/Users/") ## path to the complete dataset (including the egg-laying workers)

combined_counts_brain <- read.csv("./raw_brain_mal_mozhgan.csv", header = TRUE, row.names = 1)
phenotype_data_brain <- read.csv("./ttype_brain.csv", header = TRUE, row.names = 1)

dds_brain <- DESeqDataSetFromMatrix(countData = combined_counts_brain,
                                    colData = phenotype_data_brain,
                                    design = ~ phenotype)

keep <- rowSums(counts(dds_brain) >= 10) >= 5
dds_brain <- dds_brain[keep,]

dds_brain <- DESeq(dds_brain)

dds_vst_brain <- vst(dds_brain, blind = FALSE)
vst_counts_brain <- assay(dds_vst_brain)
vst_counts_top500_brain <- vst_counts_brain[top500_genes_brain, ]

pca_vst_brain <- prcomp(t(vst_counts_top500_brain), scale. = TRUE)

pca_vst_df_brain <- data.frame(pca_vst_brain$x,
                               Phenotype = phenotype_data_brain$phenotype,
                               Sample = rownames(phenotype_data_brain))

pca_vst_df_brain$Label <- ifelse(pca_vst_df_brain$Sample == "X87BmalW83", "A",
                                 ifelse(pca_vst_df_brain$Sample == "X125BmalW99", "B", ""))

p_pca_vst_brain <- ggplot(pca_vst_df_brain, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), vjust = -1.2, size = 5, color = "black", fontface = "bold") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("queen" = "dodgerblue",
                                "worker" = "#ff1493",
                                "foundress" = "#969696")) +
  labs(title = "Brain transcriptome",
       x = paste0("PC1 (", round(summary(pca_vst_brain)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_vst_brain)$importance[2,2] * 100, 1), "%)"))

ggsave("PCA_top500_vst_brain.jpg", plot = p_pca_vst_brain, width = 7, height = 5, dpi = 300)
ggsave("PCA_top500_vst_brain.pdf", plot = p_pca_vst_brain, width = 7, height = 5)


# ======================================
# === 3. Combine Plots
# ======================================
p_pca_vst_fatbody <- p_pca_vst_fatbody + theme(legend.position = "bottom")
p_pca_vst_brain <- p_pca_vst_brain + theme(legend.position = "none")

combined_pca <- p_pca_vst_fatbody + p_pca_vst_brain +
  plot_annotation(title = "PCA plots of Fat Body and Brain Transcriptomes") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(combined_pca)

ggsave("Combined_PCA_vst.jpg", plot = combined_pca, width = 12, height = 6, dpi = 300)
ggsave("Combined_PCA_vst.pdf", plot = combined_pca, width = 12, height = 6)


###################publication
# ======================================
# === 3. Combine Plots (side-by-side, legend centered)
# ======================================

library(patchwork)

# Remove legend from brain, keep it in fat body
p_pca_vst_fatbody <- p_pca_vst_fatbody + theme(legend.position = "right",
                                               legend.title = element_blank(),
                                               legend.text = element_text(size = 12))
p_pca_vst_brain <- p_pca_vst_brain + theme(legend.position = "none")

# Create an empty plot for spacing (optional)
empty <- ggplot() + theme_void()

# Arrange plots: fatbody | legend | brain
combined_pca <- (p_pca_vst_fatbody + plot_spacer() + p_pca_vst_brain) + 
  plot_layout(ncol = 3, widths = c(4, 1, 4), guides = "collect") & 
  theme(legend.position = "right")  # legend will appear in middle column

# Remove main title
combined_pca <- combined_pca & theme(plot.title = element_blank())
print(combined_pca)
# Save
ggsave("Combined_PCA_vst.jpg", plot = combined_pca, width = 14, height = 6, dpi = 300)
ggsave("Combined_PCA_vst.pdf", plot = combined_pca, width = 14, height = 6)










library(patchwork)
library(cowplot)  # for get_legend

# Extract the legend from the fatbody plot
legend <- get_legend(
  p_pca_vst_fatbody + theme(legend.position = "right",
                            legend.title = element_blank(),
                            legend.text = element_text(size = 12))
)

# Remove legends from both plots
p_pca_vst_fatbody_clean <- p_pca_vst_fatbody + theme(legend.position = "none") +
  labs(title = "a) Fat bodies")
p_pca_vst_brain_clean <- p_pca_vst_brain + theme(legend.position = "none") +
  labs(title = "b) Brains")

# Combine plots: fatbody | legend | brain
combined_pca <- plot_grid(
  p_pca_vst_fatbody_clean,
  legend,
  p_pca_vst_brain_clean,
  ncol = 3,
  rel_widths = c(4, 1, 4),
  align = "v"
)
print(combined_pca)
# Save
ggsave("Combined_PCA_vst.jpg", plot = combined_pca, width = 14, height = 6, dpi = 300)
ggsave("Combined_PCA_vst.pdf", plot = combined_pca, width = 14, height = 6)

# Center titles
p_pca_vst_fatbody_clean <- p_pca_vst_fatbody_clean + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

p_pca_vst_brain_clean <- p_pca_vst_brain_clean + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
combined_pca <- plot_grid(
  p_pca_vst_fatbody_clean,
  legend,
  p_pca_vst_brain_clean,
  ncol = 3,
  rel_widths = c(4, 1, 4),
  align = "v"
)

ggsave("Combined_PCA_vst.jpg", plot = combined_pca, width = 14, height = 6, dpi = 300)
ggsave("Combined_PCA_vst.pdf", plot = combined_pca, width = 14, height = 6)

# ======================================
# === Save Final Combined PCA (High Resolution)
# ======================================

# Define output filenames
pdf_file <- "Combined_PCA_vst_highres.pdf"
png_file <- "Combined_PCA_vst_highres.png"

# Save as high-resolution PDF (vector format, no pixelation)
ggsave(filename = pdf_file,
       plot = combined_pca,
       width = 14,
       height = 6,
       units = "in",
       device = cairo_pdf,  # ensures text is embedded properly
       dpi = 600)           # not necessary for PDF, but harmless

# Save as high-resolution PNG (for digital display)
ggsave(filename = png_file,
       plot = combined_pca,
       width = 14,
       height = 6,
       units = "in",
       dpi = 600,
       type = "cairo-png")  # improves text rendering
getwd()


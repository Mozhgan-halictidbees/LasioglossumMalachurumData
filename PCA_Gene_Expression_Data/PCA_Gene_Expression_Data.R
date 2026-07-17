#########################################
# PCA Analysis for Fat Body and Brain
#########################################

library(DESeq2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

setwd("C:/Users/Mozhgan/Desktop/malachurum")

# ==========================================================
# Function to perform DESeq2 normalization, select top500,
# perform VST on ALL genes, and PCA on the saved top500 genes
# ==========================================================

run_pca <- function(count_file,
                    phenotype_file,
                    sampleA,
                    sampleB,
                    title,
                    output_prefix){
  
  #-----------------------------
  # Read data
  #-----------------------------
  counts <- read.csv(count_file,
                     header = TRUE,
                     row.names = 1)
  
  pheno <- read.csv(phenotype_file,
                    header = TRUE,
                    row.names = 1)
  
  #-----------------------------
  # Create DESeq object
  #-----------------------------
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = pheno,
    design = ~ phenotype
  )
  
  #-----------------------------
  # Filter low-count genes
  #-----------------------------
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds <- dds[keep, ]
  
  #-----------------------------
  # DESeq normalization
  #-----------------------------
  dds <- DESeq(dds)
  
  #-----------------------------
  # Normalized counts (ALL genes)
  #-----------------------------
  norm_counts <- counts(dds,
                        normalized = TRUE)
  
  #-----------------------------
  # Select top500 variable genes
  #-----------------------------
  gene_variances <- apply(norm_counts,1,var)
  
  top500_genes <-
    names(sort(gene_variances,
               decreasing = TRUE))[1:500]
  
  cat("\n")
  cat("Selected top 500 genes for",title,"\n")
  
  #-----------------------------
  # Save gene list
  #-----------------------------
  write.table(
    top500_genes,
    paste0(output_prefix,"_Top500Genes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  #-----------------------------
  # VST on ALL genes
  #-----------------------------
  dds_vst <- vst(dds,
                 blind = FALSE)
  
  vst_counts <- assay(dds_vst)
  
  #-----------------------------
  # Keep only saved top500 genes
  #-----------------------------
  vst_top500 <- vst_counts[top500_genes, ]
  
  #-----------------------------
  # PCA
  #-----------------------------
  pca <- prcomp(
    t(vst_top500),
    scale. = TRUE
  )
  
  pca_df <- data.frame(
    pca$x,
    Phenotype = pheno$phenotype,
    Sample = rownames(pheno)
  )
  
  pca_df$Label <-
    ifelse(
      pca_df$Sample == sampleA,
      "A",
      ifelse(
        pca_df$Sample == sampleB,
        "B",
        ""
      )
    )
  
  #-----------------------------
  # PCA plot
  #-----------------------------
  p <- ggplot(
    pca_df,
    aes(PC1,
        PC2,
        color = Phenotype)
  ) +
    geom_point(size = 4) +
    geom_text(
      aes(label = Label),
      vjust = -1.2,
      size = 5,
      color = "black",
      fontface = "bold"
    ) +
    theme_minimal(base_size = 14) +
    scale_color_manual(
      values = c(
        queen = "dodgerblue",
        worker = "#ff1493",
        foundress = "#969696"
      )
    ) +
    labs(
      title = title,
      x = paste0(
        "PC1 (",
        round(
          summary(pca)$importance[2,1] * 100,
          1
        ),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(
          summary(pca)$importance[2,2] * 100,
          1
        ),
        "%)"
      )
    )
  
  ggsave(
    paste0(output_prefix,"_PCA.jpg"),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    paste0(output_prefix,"_PCA.pdf"),
    plot = p,
    width = 7,
    height = 5
  )
  
  return(p)
  
}

###########################################################
# FAT BODY
###########################################################

p_pca_vst_fatbody <- run_pca(
  
  count_file =
    "raw_fatbody_mal_mozhgan.csv",
  
  phenotype_file =
    "ttype_fatbody.csv",
  
  sampleA =
    "X109FmalW83",
  
  sampleB =
    "X128FmalW99",
  
  title =
    "Fat body transcriptome",
  
  output_prefix =
    "FatBody"
  
)

print(p_pca_vst_fatbody)

###########################################################
# BRAIN
###########################################################

p_pca_vst_brain <- run_pca(
  
  count_file =
    "raw_brain_mal_mozhgan.csv",
  
  phenotype_file =
    "ttype_brain.csv",
  
  sampleA =
    "X87BmalW83",
  
  sampleB =
    "X125BmalW99",
  
  title =
    "Brain transcriptome",
  
  output_prefix =
    "Brain"
  
)

print(p_pca_vst_brain)
#################################



############################################################
# PUBLICATION FIGURE
############################################################

library(cowplot)
library(patchwork)

# ----------------------------------------------------------
# Extract legend from fat body plot
# ----------------------------------------------------------

legend <- get_legend(
  
  p_pca_vst_fatbody +
    
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
  
)

# ----------------------------------------------------------
# Remove legends and add panel titles
# ----------------------------------------------------------

p_pca_vst_fatbody_clean <-
  
  p_pca_vst_fatbody +
  
  labs(title = "a) Fat body") +
  
  theme(
    
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    )
    
  )


p_pca_vst_brain_clean <-
  
  p_pca_vst_brain +
  
  labs(title = "b) Brain") +
  
  theme(
    
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    )
    
  )

# ----------------------------------------------------------
# Combine plots
# ----------------------------------------------------------

combined_pca <- plot_grid(
  
  p_pca_vst_fatbody_clean,
  
  legend,
  
  p_pca_vst_brain_clean,
  
  ncol = 3,
  
  rel_widths = c(4,1,4),
  
  align = "v"
  
)

print(combined_pca)

############################################################
# SAVE FIGURE
############################################################

ggsave(
  
  filename = "Combined_PCA_vst_highres.pdf",
  
  plot = combined_pca,
  
  width = 14,
  
  height = 6,
  
  units = "in",
  
  device = cairo_pdf,
  
  dpi = 600
  
)

ggsave(
  
  filename = "Combined_PCA_vst_highres.png",
  
  plot = combined_pca,
  
  width = 14,
  
  height = 6,
  
  units = "in",
  
  dpi = 600,
  
  type = "cairo-png"
  
)

ggsave(
  
  filename = "Combined_PCA_vst_highres.jpg",
  
  plot = combined_pca,
  
  width = 14,
  
  height = 6,
  
  units = "in",
  
  dpi = 600
  
)



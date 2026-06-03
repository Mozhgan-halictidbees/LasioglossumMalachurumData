#####module preservation, my data reference, Andreias data test
# Load required libraries
library(DESeq2)
library(WGCNA)
#lets prepare the test:
# ---------------------------
# Read count data
# ---------------------------
setwd("C:/Users/Mozhgan/Desktop/preservation_andreia/fatbody_queen_and_worker")
combined_counts <- read.csv("./merged_output.csv", h = TRUE, row.names = 1)
phenotype_data <- read.csv("./ttype_fatbody.csv", h = TRUE, row.names = 1)
dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                              colData = phenotype_data,
                              design = ~ phenotype)
smallestGroupSize <- 9
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


# Perform DESeq normalization
dds <- DESeq(dds)

# Perform variance stabilizing transformation
dds_vst <- vst(dds, blind = FALSE)


# Get the transformed expression values
datExpr <- assay(dds_vst)

# Check the class of vst_expression
class(datExpr)

class (datExpr)

datExpr <- as.data.frame(t(datExpr))
datExpr_test <- datExpr

##########
#mozhgan: to see the outliers
# ---------------------------
# Check for missing values
# ---------------------------
gsg <- goodSamplesGenes(datExpr_test, verbose = 3)
if (!gsg$allOK) {
  # If there are problematic genes or samples, remove them
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr_test)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr_test)[!gsg$goodSamples], collapse = ", ")))
  datExpr_test <- datExpr_test[gsg$goodSamples, gsg$goodGenes]
}

# ---------------------------
# Sample clustering to detect outliers
# ---------------------------
sampleTree <- hclust(dist(datExpr_test), method = "average")

# ---------------------------
# Plot the sample dendrogram
# ---------------------------
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#########################
# Module Preservation: ref (reference) vs test (test)
#########################

# Load required libraries
library(DESeq2)
library(WGCNA)

# Allow multi-threading
enableWGCNAThreads()




####################################
# Set working directory
##now preparing the reference network
setwd("C:/Users/Mozhgan/Desktop/WGCNA_MOM/eye/malachurum_fatbody/mal_mozhgan_excluding_2_egg_laying_workers")

# Load packages
library(DESeq2)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggrepel)
library(sva)  # For ComBat

# ---------------------------
# Read count and phenotype data
# ---------------------------
combined_counts <- read.csv("./raw_fatbody_mal_mozhgan.csv", h = TRUE, row.names = 1)
phenotype_data <- read.csv("./ttype_fatbody.csv", h = TRUE, row.names = 1)
dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                              colData = phenotype_data,
                              design = ~ phenotype)
smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


# Perform DESeq normalization
dds <- DESeq(dds)

# Perform variance stabilizing transformation
dds_vst <- vst(dds, blind = FALSE)


# Get the transformed expression values
datExpr <- assay(dds_vst)

# Check the class of vst_expression
class(datExpr)

class (datExpr)

datExpr <- as.data.frame(t(datExpr))
datExpr_ref <- datExpr
load("./mal-expt1-01-dataInput_core.RData")
load("./mal-expt1-NetworkConstruction-core.RData")
TOM <- readRDS("./TOM_signed.rds")

sampleTree <- hclust(dist(datExpr_ref), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="")


#####################
powers = c(1:10, seq(12, 30, 2))
sft = pickSoftThreshold(datExpr_ref, powerVector = powers, verbose = 5)

# Visualize scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale-free topology R^2",
     type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.85, col="red")



softPower <- 17  # replace with your chosen power
adjacency <- adjacency(datExpr_ref, power = softPower, type = "signed")
#TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# Cluster genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamic tree cut
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Merge similar modules
MEList <- moduleEigengenes(datExpr_ref, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

MEDissThres <- 0.25  # merge modules with correlation > 0.75
merge <- mergeCloseModules(datExpr_ref, dynamicColors, cutHeight = MEDissThres, verbose = 3)
moduleColors <- merge$colors
MEs <- merge$newMEs


#####visulaization
##### Visualization
# -------------------------------
# Step 1: Prepare traits for correlation
# -------------------------------
datTraits_ref = read.csv("./datTraits.csv", h=T, row.names=1)
rownames(datTraits) <- datTraits$SampleID



# -------------------------------
# Step 2: Correlate modules with traits
# -------------------------------
library(WGCNA)
nSamples <- nrow(datExpr_ref)

# Calculate robust correlations between module eigengenes and traits
moduleTraitCor <- bicor(MEs, datTraits_ref, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# -------------------------------
# Step 3: Format text for the heatmap
# -------------------------------
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

# -------------------------------
# Step 4: Plot the heatmap with improved layout
# -------------------------------
# Set plot window size
sizeGrWindow(width = 12, height = 8)  # wider and taller

# Adjust margins: c(bottom, left, top, right)
par(mar = c(6, 8, 4, 2))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits_ref),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1,1),
  main = "Module-trait relationships",
  xLabelsAngle = 0  # makes labels horizontal
)
##### Save as PNG
png(filename = "module_trait_heatmap.png", width = 1600, height = 1200, res = 150)  # adjust size & resolution
par(mar = c(6, 8, 4, 2))  # margins
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits_ref),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1,1),
  main = "Module-trait relationships",
  xLabelsAngle = 0
)
dev.off()  # close PNG device

##### Save as PDF
pdf(file = "module_trait_heatmap.pdf", width = 12, height = 8)  # same size as sizeGrWindow
par(mar = c(6, 8, 4, 2))  # margins
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits_ref),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1,1),
  main = "Module-trait relationships",
  xLabelsAngle = 0
)
dev.off()  # close PDF device


####saving
save(datExpr_ref, datExpr_test, moduleColors, MEs, geneTree, file="WGCNA_ref_test_data.RData")



##############performing preservation
multiExpr <- list(
  Reference = list(data = datExpr_ref),
  Test      = list(data = datExpr_test)
)

colorList <- list(Reference = moduleColors)

# Module preservation
mp <- modulePreservation(multiExpr, colorList,
                         referenceNetworks = 1,
                         nPermutations = 200,  # increase to 1000 for robust stats
                         randomSeed = 123,
                         networkType = "signed",
                         verbose = 3)
setwd("C:/Users/Mozhgan/Desktop/preservation_andreia/fatbody_queen_and_worker")
save(mp, file = "modulePreservation_results.RData")


# View preservation statistics
ref_preservation <- mp$preservation$Z$ref.Reference
ref_preservation




pres_df <- ref_preservation$inColumnsAlsoPresentIn.Test
pres_df$module <- rownames(pres_df)


###########################################

p_fatbody <- ggplot(pres_df, aes(x = reorder(module, Zsummary.pres), y = Zsummary.pres, fill = Zsummary.pres)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "forestgreen", midpoint = 10) +
  geom_hline(yintercept = c(2,10), linetype = "dashed", color = "black") +
  labs(x = "Module", y = "Zsummary (Preservation)", title = "Module Preservation in Fat Body") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
p_fatbody

###############################################
####now do the same for brain
#####module preservation, my data reference, Andreias data test
# Load required libraries
library(DESeq2)
library(WGCNA)
#lets prepare the test:
# ---------------------------
# Read count data
# ---------------------------
setwd("C:/Users/Mozhgan/Desktop/preservation_andreia/brain_queen_and_worker")
getwd()
combined_counts <- read.csv("./merged_output.csv", h = TRUE, row.names = 1)
phenotype_data <- read.csv("./ttype_brain.csv", h = TRUE, row.names = 1)
dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                              colData = phenotype_data,
                              design = ~ phenotype)
smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


# Perform DESeq normalization
dds <- DESeq(dds)

# Perform variance stabilizing transformation
dds_vst <- vst(dds, blind = FALSE)


# Get the transformed expression values
datExpr <- assay(dds_vst)

# Check the class of vst_expression
class(datExpr)

class (datExpr)

datExpr <- as.data.frame(t(datExpr))
datExpr_test <- datExpr

##########
#mozhgan: to see the outliers
# ---------------------------
# Check for missing values
# ---------------------------
gsg <- goodSamplesGenes(datExpr_test, verbose = 3)
if (!gsg$allOK) {
  # If there are problematic genes or samples, remove them
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr_test)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr_test)[!gsg$goodSamples], collapse = ", ")))
  datExpr_test <- datExpr_test[gsg$goodSamples, gsg$goodGenes]
}

# ---------------------------
# Sample clustering to detect outliers
# ---------------------------
sampleTree <- hclust(dist(datExpr_test), method = "average")

# ---------------------------
# Plot the sample dendrogram
# ---------------------------
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#########################
# Module Preservation: ref (reference) vs test (test)
#########################

# Load required libraries
library(DESeq2)
library(WGCNA)

# Allow multi-threading
enableWGCNAThreads()




####################################
# Set working directory
##now preparing the reference network
setwd("C:/Users/Mozhgan/Desktop/WGCNA_MOM/eye/malachurum_brain/mal_mozhgan_excluding_egg_laying_workers_brain")

# Load packages
library(DESeq2)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggrepel)
library(sva)  # For ComBat

# ---------------------------
# Read count and phenotype data
# ---------------------------
combined_counts <- read.csv("./raw_brain_mal_mozhgan.csv", h = TRUE, row.names = 1)
phenotype_data <- read.csv("./ttype_brain.csv", h = TRUE, row.names = 1)
dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                              colData = phenotype_data,
                              design = ~ phenotype)
smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


# Perform DESeq normalization
dds <- DESeq(dds)

# Perform variance stabilizing transformation
dds_vst <- vst(dds, blind = FALSE)


# Get the transformed expression values
datExpr <- assay(dds_vst)

# Check the class of vst_expression
class(datExpr)

class (datExpr)

datExpr <- as.data.frame(t(datExpr))
datExpr_ref <- datExpr
load("./mal-expt1-01-dataInput_core.RData")
load("./mal-expt1-NetworkConstruction-core.RData")
TOM <- readRDS("./TOM_signed.rds")

sampleTree <- hclust(dist(datExpr_ref), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="")


#####################
powers = c(1:10, seq(12, 30, 2))
sft = pickSoftThreshold(datExpr_ref, powerVector = powers, verbose = 5)

# Visualize scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale-free topology R^2",
     type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.85, col="red")



softPower <- 16  # replace with your chosen power
adjacency <- adjacency(datExpr_ref, power = softPower, type = "signed")
#TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# Cluster genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamic tree cut
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Merge similar modules
MEList <- moduleEigengenes(datExpr_ref, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

MEDissThres <- 0.25  # merge modules with correlation > 0.75
merge <- mergeCloseModules(datExpr_ref, dynamicColors, cutHeight = MEDissThres, verbose = 3)
moduleColors <- merge$colors
MEs <- merge$newMEs


#####visulaization
##### Visualization
# -------------------------------
# Step 1: Prepare traits for correlation
# -------------------------------
datTraits_ref = read.csv("./datTraits.csv", h=T, row.names=1)
rownames(datTraits) <- datTraits$SampleID



# -------------------------------
# Step 2: Correlate modules with traits
# -------------------------------
library(WGCNA)
nSamples <- nrow(datExpr_ref)

# Calculate robust correlations between module eigengenes and traits
moduleTraitCor <- bicor(MEs, datTraits_ref, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# -------------------------------
# Step 3: Format text for the heatmap
# -------------------------------
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

# -------------------------------
# Step 4: Plot the heatmap with improved layout
# -------------------------------
# Set plot window size
sizeGrWindow(width = 12, height = 8)  # wider and taller

# Adjust margins: c(bottom, left, top, right)
par(mar = c(6, 8, 4, 2))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits_ref),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1,1),
  main = "Module-trait relationships",
  xLabelsAngle = 0  # makes labels horizontal
)
##### Save as PNG
png(filename = "module_trait_heatmap.png", width = 1600, height = 1200, res = 150)  # adjust size & resolution
par(mar = c(6, 8, 4, 2))  # margins
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits_ref),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1,1),
  main = "Module-trait relationships",
  xLabelsAngle = 0
)
dev.off()  # close PNG device

##### Save as PDF
pdf(file = "module_trait_heatmap.pdf", width = 12, height = 8)  # same size as sizeGrWindow
par(mar = c(6, 8, 4, 2))  # margins
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits_ref),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1,1),
  main = "Module-trait relationships",
  xLabelsAngle = 0
)
dev.off()  # close PDF device


####saving
save(datExpr_ref, datExpr_test, moduleColors, MEs, geneTree, file="WGCNA_ref_test_data.RData")



##############performing preservation
multiExpr <- list(
  Reference = list(data = datExpr_ref),
  Test      = list(data = datExpr_test)
)

colorList <- list(Reference = moduleColors)

# Module preservation
mp <- modulePreservation(multiExpr, colorList,
                         referenceNetworks = 1,
                         nPermutations = 200,  # increase to 1000 for robust stats
                         randomSeed = 123,
                         networkType = "signed",
                         verbose = 3)
setwd("C:/Users/Mozhgan/Desktop/preservation_andreia/brain_queen_and_worker")
save(mp, file = "modulePreservation_results.RData")


# View preservation statistics
ref_preservation <- mp$preservation$Z$ref.Reference
ref_preservation




pres_df <- ref_preservation$inColumnsAlsoPresentIn.Test
pres_df$module <- rownames(pres_df)


###########################################

p_brain <- ggplot(pres_df, aes(x = reorder(module, Zsummary.pres), y = Zsummary.pres, fill = Zsummary.pres)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "forestgreen", midpoint = 10) +
  geom_hline(yintercept = c(2,10), linetype = "dashed", color = "black") +
  labs(x = "Module", y = "Zsummary (Preservation)", title = "Module Preservation in Brain") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
p_brain


###
#combining
library(patchwork)
combined_plot <- (p_fatbody | p_brain)

# Display
combined_plot
setwd("C:/Users/Mozhgan/Desktop/preservation_andreia")
pdf("combined_plot.pdf", width = 10, height = 6, useDingbats = FALSE)
combined_plot
dev.off()
ggsave("combined_plot.png", 
       plot = combined_plot, 
       width = 10, height = 6, units = "in",
       dpi = 600)   # 600 dpi = publication quality


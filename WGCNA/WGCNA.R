

####################################
# Set working directory
setwd("C:/Users/Mozhgan/....")

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

datExpr_brain <- as.data.frame(t(datExpr))

load("./mal-expt1-01-dataInput_core.RData")
load("./mal-expt1-NetworkConstruction-core.RData")
TOM <- readRDS("./TOM_signed.rds")

#code1:
#install.packages("impute")
#install.packages("WGCNA")
#install.packages("WGCNA", dependencies = TRUE)
#BiocManager::install("impute")
library(impute)

setwd(workingDir)

library(WGCNA)
#prevent automatic conversion of character strings into factors when reading or manipulating data
options(stringsAsFactors = FALSE)


# cluster the samples to look for outliers
sampleTree = hclust(dist(datExpr), method = "average")
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# loading trait data
datTraits = read.csv("./datTraits.csv", h=T, row.names=1)
datTraits_brain <- data.frame(lapply(datTraits, as.numeric), row.names = rownames(datTraits))
dim(datTraits_brain)

names(datTraits_brain)


#collectGarbage()
gc()


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")



# Convert phenotype to factor, then to numeric (mozhgan)
#datTraits$phenotype <- as.numeric(factor(datTraits$phenotype, levels = c("queen", "foundress", "worker")))
# mozhgan: Convert non-numeric traits to numeric

#datTraits$Group_numeric <- as.numeric(factor(datTraits$Group, levels = c("Worker", "Queen")))

ncol(datExpr)
nrow(datTraits)


# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.

png("sample_dendro_trait_heatmap.png", width = 1200, height = 800, res = 150)
plotDendroAndColors(
  sampleTree2, 
  traitColors, 
  groupLabels = names(datTraits), 
  main = "Sample dendrogram and trait heatmap of brain"
)
dev.off()
save(datExpr, datTraits, file = "mal-expt1-01-dataInput_core.RData")





#code2:
#wgcna2


library(WGCNA)
options(stringsAsFactors = FALSE)

# Load the data saved in the first part
lnames = load(file = "mal-expt1-01-dataInput_core.RData")
lnames

# Choose a set of soft-thresholding powers #(for scale free topology)
powers = c(c(1:10), seq(from = 12, to=30, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.85
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# calculating co-expression similarity and adjacency (16 instead of 9 by mozhgan)
softPower = 16
adjacency = adjacency(datExpr, power = softPower, type = "signed")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed")

# Save the TOM matrix to the working directory
saveRDS(TOM, file = "TOM_signed.rds")

dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram): 
# branches of the dendogram group together densely interconnected, highly coexpressed genes
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
## The height of the branches in the dendrogram represents the level of dissimilarity at which the genes are grouped.


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
# this function returned 34 modules labeled 1-34 largest to smallest.
# label 0 is reserved for unassigned genes



# we now plot the module assignment under the gene dendograms
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col="red")

# merging modules whose expression profiles are very similar
# here choose correlation of 0.75 to merge = height cut of 0.25
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs_brain = mergedMEs
MEs_brain
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "mal-expt1-NetworkConstruction-core.RData")










# Define numbers of genes and samples
nGenes = ncol(datExpr_brain)
nSamples = nrow(datExpr_brain)

# Recalculate MEs with color labels  (eigengene=pc1)
MEs0 = moduleEigengenes(datExpr_brain, moduleColors)$eigengenes
MEs_brain = orderMEs(MEs0)

moduleTraitCor_brain = bicor(MEs_brain, datTraits_brain)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor_brain, nSamples)

# create text matrix with correlation and p-values
textMatrix_brain = paste(signif(moduleTraitCor_brain, 2), "\n(",
                         signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix_brain) = dim(moduleTraitCor_brain)

# Open JPEG device
# Create labels without the "ME" prefix
#moduleNames <- gsub("^ME", "", names(MEs))
# Create clean module names without "ME"
moduleNames <- gsub("^ME", "", names(MEs))


library(WGCNA)
setwd("C:/Users/Mozhgan/Desktop/...")
# Open PDF device
pdf("Brain_module_trait_heatmap.pdf", width = 10, height = 8)

# Adjust margins
par(mar = c(6, 10, 4, 2))  

brain_heatmap <- labeledHeatmap(
  Matrix = moduleTraitCor_brain,
  xLabels = names(datTraits_brain),
  yLabels = names(MEs_brain),
  ySymbols = gsub("^ME", "", names(MEs_brain)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_brain,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Brain module-trait relationships"
)

brain_heatmap
# Close device and save
dev.off()




####################################
# Set working directory
setwd("C:/Users/Mozhgan/.....")

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

datExpr_fatbody <- as.data.frame(t(datExpr))

load("./mal-expt1-01-dataInput_core.RData")
load("./mal-expt1-NetworkConstruction-core.RData")
TOM <- readRDS("./TOM_signed.rds")

#code1:
#install.packages("impute")
#install.packages("WGCNA")
#install.packages("WGCNA", dependencies = TRUE)
#BiocManager::install("impute")
library(impute)

setwd(workingDir)

library(WGCNA)
#prevent automatic conversion of character strings into factors when reading or manipulating data
options(stringsAsFactors = FALSE)


# cluster the samples to look for outliers
sampleTree = hclust(dist(datExpr), method = "average")
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# loading trait data
datTraits = read.csv("./datTraits.csv", h=T, row.names=1)
datTraits_fatbody <- data.frame(lapply(datTraits, as.numeric), row.names = rownames(datTraits))
dim(datTraits_fatbody)

names(datTraits_fatbody)


#collectGarbage()
gc()


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")



# Convert phenotype to factor, then to numeric (mozhgan)
#datTraits$phenotype <- as.numeric(factor(datTraits$phenotype, levels = c("queen", "foundress", "worker")))
# mozhgan: Convert non-numeric traits to numeric

#datTraits$Group_numeric <- as.numeric(factor(datTraits$Group, levels = c("Worker", "Queen")))

ncol(datExpr)
nrow(datTraits)


# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.

png("sample_dendro_trait_heatmap.png", width = 1200, height = 800, res = 150)
plotDendroAndColors(
  sampleTree2, 
  traitColors, 
  groupLabels = names(datTraits), 
  main = "Sample dendrogram and trait heatmap of fatbody"
)
dev.off()
save(datExpr, datTraits, file = "mal-expt1-01-dataInput_core.RData")





#code2:
#wgcna2


library(WGCNA)
options(stringsAsFactors = FALSE)

# Load the data saved in the first part
lnames = load(file = "mal-expt1-01-dataInput_core.RData")
lnames

# Choose a set of soft-thresholding powers #(for scale free topology)
powers = c(c(1:10), seq(from = 12, to=30, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.85
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# calculating co-expression similarity and adjacency (16 instead of 9 by mozhgan)
softPower = 16
adjacency = adjacency(datExpr, power = softPower, type = "signed")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed")

# Save the TOM matrix to the working directory
saveRDS(TOM, file = "TOM_signed.rds")

dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram): 
# branches of the dendogram group together densely interconnected, highly coexpressed genes
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
## The height of the branches in the dendrogram represents the level of dissimilarity at which the genes are grouped.


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
# this function returned 34 modules labeled 1-34 largest to smallest.
# label 0 is reserved for unassigned genes



# we now plot the module assignment under the gene dendograms
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col="red")

# merging modules whose expression profiles are very similar
# here choose correlation of 0.75 to merge = height cut of 0.25
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs_fatbody = mergedMEs
MEs_fatbody
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "mal-expt1-NetworkConstruction-core.RData")










# Define numbers of genes and samples
nGenes = ncol(datExpr_brain)
nSamples = nrow(datExpr_brain)

# Recalculate MEs with color labels  (eigengene=pc1)
MEs0 = moduleEigengenes(datExpr_fatbody, moduleColors)$eigengenes
MEs_fatbody = orderMEs(MEs0)

moduleTraitCor_fatbody = bicor(MEs_fatbody, datTraits_fatbody)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor_fatbody, nSamples)

# create text matrix with correlation and p-values
textMatrix_fatbody = paste(signif(moduleTraitCor_fatbody, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix_fatbody) = dim(moduleTraitCor_fatbody)

# Open JPEG device
# Create labels without the "ME" prefix
#moduleNames <- gsub("^ME", "", names(MEs))
# Create clean module names without "ME"
moduleNames <- gsub("^ME", "", names(MEs))


library(WGCNA)
setwd("C:/Users/Mozhgan/Desktop")
pdf("Fatbody_module_trait_heatmap.pdf", width = 10, height = 8)
# Fat body heatmap
fatbody_heatmap <- labeledHeatmap(
  Matrix = moduleTraitCor_fatbody,
  xLabels = names(datTraits_fatbody),
  yLabels = names(MEs_fatbody),
  ySymbols = gsub("^ME", "", names(MEs_fatbody)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_fatbody,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Fat body module-trait relationships"
)




library(gridExtra)

# Open a high-resolution PDF or PNG
pdf("Module_Trait_Heatmaps_Combined.pdf", width = 12, height = 14)

grid.arrange(
  grobs = list(brain_heatmap$gtable, fatbody_heatmap$gtable),
  ncol = 1,  # stack vertically
  heights = c(1, 1)  # adjust relative height if needed
)

dev.off()
getwd()



png("Module_Trait_Heatmaps_Combined.png", width = 3000, height = 4000, res = 300)

grid.arrange(
  grobs = list(brain_heatmap$gtable, fatbody_heatmap$gtable),
  ncol = 1,
  heights = c(1, 1)
)

dev.off()







library(gridExtra)
library(grid)

# Capture labeledHeatmap output as grob
brain_grob <- grid.grabExpr({
  labeledHeatmap(
    Matrix = moduleTraitCor_brain,
    xLabels = names(datTraits_brain),
    yLabels = names(MEs_brain),
    ySymbols = gsub("^ME", "", names(MEs_brain)),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix_brain,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1, 1),
    main = "Brain module-trait relationships"
  )
})

fatbody_grob <- grid.grabExpr({
  labeledHeatmap(
    Matrix = moduleTraitCor_fatbody,
    xLabels = names(datTraits_fatbody),
    yLabels = names(MEs_fatbody),
    ySymbols = gsub("^ME", "", names(MEs_fatbody)),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix_fatbody,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1, 1),
    main = "Fat body module-trait relationships"
  )
})

# Save combined figure
pdf("Module_Trait_Heatmaps_Combined.pdf", width = 12, height = 14)
grid.arrange(brain_grob, fatbody_grob, ncol = 1)
dev.off()



pdf("Module_Trait_Heatmaps_Combined.pdf", width = 12, height = 14)

# Split the plotting area into 2 rows
par(mfrow = c(2, 1), mar = c(6, 10, 4, 2))  # adjust margins

# Brain heatmap
labeledHeatmap(
  Matrix = moduleTraitCor_brain,
  xLabels = names(datTraits_brain),
  yLabels = names(MEs_brain),
  ySymbols = gsub("^ME", "", names(MEs_brain)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_brain,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Brain module-trait relationships"
)

# Fat body heatmap
labeledHeatmap(
  Matrix = moduleTraitCor_fatbody,
  xLabels = names(datTraits_fatbody),
  yLabels = names(MEs_fatbody),
  ySymbols = gsub("^ME", "", names(MEs_fatbody)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_fatbody,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Fat body module-trait relationships"
)

# Open PDF device
pdf("Brain_module_trait_heatmap.pdf", width = 10, height = 8)

# Adjust margins
par(mar = c(6, 10, 4, 2))  

# Plot the heatmap
labeledHeatmap(
  Matrix = moduleTraitCor_brain,
  xLabels = names(datTraits_brain),
  yLabels = names(MEs_brain),
  ySymbols = moduleNames,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_brain,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Brain module-trait relationships"
)

# Close device and save
dev.off()




png("Brain_module_trait_heatmap.png", width = 900, height = 1000, res = 150)

par(mar = c(6, 10, 4, 2))  

labeledHeatmap(
  Matrix = moduleTraitCor_brain,
  xLabels = names(datTraits_brain),
  yLabels = names(MEs_brain),
  ySymbols = moduleNames,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_brain,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Brain module-trait relationships"
)

dev.off()



# Load WGCNA
library(WGCNA)

# -------------------------------
# Brain module-trait heatmap
# -------------------------------

# Open PDF device
pdf("Brain_module_trait_heatmap.pdf", width = 10, height = 10)

# Adjust margins: increase left margin for module names
par(mar = c(6, 10, 4, 2))  

# Plot the heatmap
labeledHeatmap(
  Matrix = moduleTraitCor_brain,
  xLabels = names(datTraits_brain),
  yLabels = names(MEs_brain),
  ySymbols = gsub("^ME", "", names(MEs_brain)), # clean module names
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_brain,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Brain module-trait relationships"
)

dev.off()


# High-resolution PNG
png("Brain_module_trait_heatmap.png", width = 1800, height = 2400, res = 300)
par(mar = c(6, 10, 4, 2))  

labeledHeatmap(
  Matrix = moduleTraitCor_brain,
  xLabels = names(datTraits_brain),
  yLabels = names(MEs_brain),
  ySymbols = gsub("^ME", "", names(MEs_brain)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_brain,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Brain module-trait relationships"
)

dev.off()


# -------------------------------
# Fat body module-trait heatmap
# -------------------------------

# Open PDF device
pdf("FatBody_module_trait_heatmap.pdf", width = 10, height = 10)

# Adjust margins: increase left margin for module names
par(mar = c(6, 10, 4, 2))  

# Plot the heatmap
labeledHeatmap(
  Matrix = moduleTraitCor_fatbody,
  xLabels = names(datTraits_fatbody),
  yLabels = names(MEs_fatbody),
  ySymbols = gsub("^ME", "", names(MEs_fatbody)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_fatbody,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Fat body module-trait relationships"
)

dev.off()


# High-resolution PNG
png("FatBody_module_trait_heatmap.png", width = 1800, height = 2400, res = 300)
par(mar = c(6, 10, 4, 2))  

labeledHeatmap(
  Matrix = moduleTraitCor_fatbody,
  xLabels = names(datTraits_fatbody),
  yLabels = names(MEs_fatbody),
  ySymbols = gsub("^ME", "", names(MEs_fatbody)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_fatbody,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Fat body module-trait relationships"
)

dev.off()



library(png)
library(grid)
library(gridExtra)

# Read the two PNG files
brain_img <- readPNG("Brain_module_trait_heatmap.png")
fatbody_img <- readPNG("FatBody_module_trait_heatmap.png")

# Convert to raster grobs
brain_grob <- rasterGrob(brain_img, interpolate=TRUE)
fatbody_grob <- rasterGrob(fatbody_img, interpolate=TRUE)

# Arrange horizontally
pdf("Combined_Module_Trait_Heatmaps_Horizontal.pdf", width=24, height=12)  # wider for side-by-side
grid.arrange(brain_grob, fatbody_grob, ncol=2, widths=c(1,1))
dev.off()





###############
library(png)
library(grid)
library(gridExtra)

# Read the two PNG files
brain_img <- readPNG("Brain_module_trait_heatmap.png")
fatbody_img <- readPNG("FatBody_module_trait_heatmap.png")

# Convert to raster grobs
brain_grob <- rasterGrob(brain_img, interpolate=TRUE)
fatbody_grob <- rasterGrob(fatbody_img, interpolate=TRUE)

# Arrange horizontally with reduced spacing
combined <- arrangeGrob(
  brain_grob, fatbody_grob,
  ncol = 2,
  widths = c(1, 1),
  padding = unit(0.01, "cm")  # smaller gap between images
)

combined <- arrangeGrob(
  brain_grob, fatbody_grob,
  ncol = 2,
  widths = unit.c(unit(1, "npc")*0.5, unit(1, "npc")*0.5)  # exactly half-half
)
grid.arrange(
  brain_grob, fatbody_grob,
  ncol = 2,
  widths = c(1, 1),
  respect = TRUE  # keeps proportions tight
)

# Save as PDF
pdf("Combined_Module_Trait_Heatmaps_Horizontal.pdf", width=24, height=12)
grid.draw(combined)
dev.off()

# Save as PNG
png("Combined_Module_Trait_Heatmaps_Horizontal.png", width=2400, height=1200, res=150)
grid.draw(combined)
dev.off()

library(png)
library(grid)
library(gridExtra)

# Read the two PNG files
brain_img <- readPNG("Brain_module_trait_heatmap.png")
fatbody_img <- readPNG("FatBody_module_trait_heatmap.png")

# Convert to raster grobs
brain_grob <- rasterGrob(brain_img, interpolate=FALSE)  # set FALSE for crispness
fatbody_grob <- rasterGrob(fatbody_img, interpolate=FALSE)

# Arrange horizontally with almost no gap
combined <- arrangeGrob(
  brain_grob, fatbody_grob,
  ncol = 2,
  widths = unit.c(unit(1, "npc")*0.5, unit(1, "npc")*0.5)
)

# Save as PDF
pdf("Combined_Module_Trait_Heatmaps_Horizontal.pdf", width=24, height=12)
grid.draw(combined)
dev.off()

# Save as high-quality PNG
png("Combined_Module_Trait_Heatmaps_Horizontal.png",
    width=4800, height=2400, res=300)  # doubled pixels + high DPI
grid.draw(combined)
dev.off()






























png("Brain_module_trait_heatmap.png", width = 900, height = 1000, res = 150)

par(mar = c(6, 10, 4, 2))  

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = moduleNames,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Brain module-trait relationships"
)

dev.off()

############################

# Count how many genes are in each module
table(moduleColors)
# Replace with your gene name of interest
gene_of_interest <- "LMAL_12084"

# Find its module
module_of_gene <- moduleColors[which(names(moduleColors) == gene_of_interest)]

cat("Gene", gene_of_interest, "belongs to module:", module_of_gene, "\n")

head(names(moduleColors))
length(moduleColors)
names(moduleColors) <- colnames(datExpr)
head(names(moduleColors))
###########################################

datExpr <- t(datExpr)
# Select genes in the steelyellow module
module <- "yellow"
genes_in_module <- rownames(datExpr)[moduleColors == module]

# Check how many genes
length(genes_in_module)

# Save gene names to a CSV file
write.csv(genes_in_module, file = paste0("genes_", module, ".csv"),
          row.names = FALSE, quote = FALSE)


###########
library(pheatmap)
#datExpr <- t(datExpr)
# datExpr: samples x genes (as in your original), we transpose to genes x samples:
datExpr_mat <- as.matrix(datExpr)
datExpr_mat <- t(datExpr_mat)
mode(datExpr_mat) <- "numeric"

module <- "yellow"
module_genes_idx <- which(moduleColors == module)
genes_in_module <- rownames(datExpr_mat)[module_genes_idx]

datExpr_module <- datExpr_mat[module_genes_idx, , drop = FALSE]
datExpr_module_scaled <- t(scale(t(datExpr_module)))

# Prepare annotation dataframe: ensure rownames are sample names
ann_col <- datTraits
# If datTraits has SampleID column, set rownames accordingly:
if("SampleID" %in% colnames(ann_col)) {
  rownames(ann_col) <- ann_col$SampleID
}

# Keep only samples present in expression matrix (and in same order initially)
ann_col <- ann_col[intersect(rownames(ann_col), colnames(datExpr_module_scaled)), , drop = FALSE]

# --- Create a single Group column (one label per sample) ---
# Assumes datTraits has binary columns named exactly: worker, foundress, queen
# Priority: worker > foundress > queen (change if you prefer!)
ann_col$Group <- NA_character_
ann_col$Group[ann_col$worker == 1] <- "worker"
ann_col$Group[is.na(ann_col$Group) & ann_col$foundress == 1] <- "foundress"
ann_col$Group[is.na(ann_col$Group) & ann_col$queen == 1] <- "queen"

# If any samples remain NA (no flag), label them "other"
ann_col$Group[is.na(ann_col$Group)] <- "other"

# Convert to factor with desired ordering for plotting
ann_col$Group <- factor(ann_col$Group, levels = c("worker", "foundress", "queen", "other"))

# --- Build ordered sample list: all workers, then foundresses, then queens, then others ---
worker_samples <- rownames(ann_col)[ann_col$Group == "worker"]
foundress_samples <- rownames(ann_col)[ann_col$Group == "foundress"]
queen_samples <- rownames(ann_col)[ann_col$Group == "queen"]
other_samples <- rownames(ann_col)[ann_col$Group == "other"]

# Use unique() & intersect to avoid duplicates and ensure they exist in the matrix columns
ordered_samples <- unique(c(worker_samples, foundress_samples, queen_samples, other_samples))
ordered_samples <- ordered_samples[ordered_samples %in% colnames(datExpr_module_scaled)]

# Reorder expression matrix and annotation to that exact order
datExpr_module_scaled <- datExpr_module_scaled[, ordered_samples, drop = FALSE]
ann_col <- ann_col[ordered_samples, , drop = FALSE]

# --- Define annotation colours exactly as requested ---
annotation_colors <- list(
  Group = c(
    worker = "#ff1493",    # your worker color
    foundress = "#969696", # your foundress color
    queen = "dodgerblue",
    other = "white"        # color for any unlabeled samples
  )
)

# Optional: save gene names
write.csv(genes_in_module, file = paste0("genes_", module, ".csv"),
          row.names = FALSE, quote = FALSE)
dev.off()

# Plot heatmap
pheatmap(datExpr_module_scaled,
         color = colorRampPalette(c("green", "black", "#FF2400"))(75),
         cluster_rows = TRUE,       # cluster genes
         cluster_cols = FALSE,      # keep samples in the provided order
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = ann_col["Group", drop = FALSE],  # only show Group in the annotation bar
         annotation_colors = annotation_colors,
         main = paste0("Module: ", module, " fat body"))

# Save heatmap as PNG

png("module_yellow.png", width = 3000, height = 2000, res = 300)
pheatmap(datExpr_module_scaled,
         color = colorRampPalette(c("green", "black", "red"))(75),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = ann_col["Group", drop = FALSE],
         annotation_colors = annotation_colors,
         main = paste0("Module: ", module, " fat body"))


dev.off()





library(pheatmap)

# Create the pheatmap object
heatmap_obj <- pheatmap(datExpr_module_scaled,
                        color = colorRampPalette(c("green", "black", "red"))(75),
                        cluster_rows = TRUE,
                        cluster_cols = FALSE,
                        show_rownames = FALSE,
                        show_colnames = TRUE,
                        annotation_col = ann_col["Group", drop = FALSE],
                        annotation_colors = annotation_colors,
                        main = paste0("Module: ", module, " fat body"))

# Save the plot as PNG
png("module_yellow.png", width = 3000, height = 2000, res = 300)
grid::grid.newpage()          # start a new page
grid::grid.draw(heatmap_obj$gtable)  # draw the heatmap
dev.off()

while(dev.cur() > 1) dev.off()

######################high KME genes
# Load WGCNA if not already
library(WGCNA)
datExpr = t(datExpr)
# Module of interest
module <- "yellow"
ME_module <- MEs[, paste0("ME", module)]

# Indices of genes in the module
module_genes_idx <- which(moduleColors == module)
genes_in_module <- rownames(datExpr)[module_genes_idx]

# Expression matrix for module genes
gene_expr_module <- datExpr[module_genes_idx, , drop = FALSE]

# Calculate KME using bicor, handle missing values
KME_values <- apply(gene_expr_module, 1, function(gene) {
  bicor(gene, ME_module, use = "pairwise.complete.obs")
})

# Select genes with KME > 0.7
hub_genes <- genes_in_module[KME_values > 0.7]

# Create data frame with gene names + KME
hub_genes_df <- data.frame(
  Gene = hub_genes,
  KME = KME_values[KME_values > 0.7]
)

# Save to CSV
write.csv(hub_genes_df, file = paste0("hub_genes_", module, ".csv"),
          row.names = FALSE, quote = FALSE)

cat("Number of hub genes in", module, "with bicor KME > 0.7:", length(hub_genes), "\n")
# Get number of genes per module
table(moduleColors)
dev.off()



######hub gene
# Select genes in the 'yellow' module
datExpr <- t(datExpr)  
chooseTopHubInEachModule_fixed <- function(datExpr, moduleColors, power, type = "signed", colorh = NULL) {
  if (is.null(colorh)) colorh <- unique(moduleColors)
  topHubs <- sapply(colorh, function(m) {
    inModule <- (moduleColors == m)
    datModule <- datExpr[, inModule, drop = FALSE]
    adj <- adjacency(datModule, power = power, type = type)
    kWithin <- rowSums(adj, na.rm = TRUE)
    topHub <- names(which.max(kWithin))
    return(topHub)
  })
  names(topHubs) <- colorh
  return(topHubs)
}
top_hubs <- chooseTopHubInEachModule_fixed(
  datExpr = datExpr,
  moduleColors = moduleColors,
  power = 17,
  type = "signed",
  colorh = "yellow"
)

print(top_hubs)

dim(datExpr)
length(moduleColors)
############################GS vs MM
library(WGCNA)

# -----------------------------
# Define module and trait
# -----------------------------
module <- "yellow"  # change to your module
trait <- "worker"         # change to "worker" or "foundress"

# -----------------------------
# Module Membership (MM)
# -----------------------------
ME_col <- paste0("ME", module)
MM_module <- bicor(datExpr, MEs[, ME_col], use = "pairwise.complete.obs")

# -----------------------------
# Gene Significance (GS)
# -----------------------------
traitVec <- as.numeric(datTraits[[trait]])
GS_trait <- bicor(datExpr, traitVec, use = "pairwise.complete.obs")

# -----------------------------
# Combine and subset module genes
# -----------------------------
moduleGenes <- data.frame(MM = MM_module[,1],
                          GS = GS_trait[,1],
                          Module = moduleColors)
moduleGenes <- moduleGenes[moduleGenes$Module == module, ]
moduleGenes <- moduleGenes[complete.cases(moduleGenes), ]

# -----------------------------
# Plot GS vs MM
# -----------------------------
plot(moduleGenes$MM, moduleGenes$GS,
     xlab = paste("Module Membership in", module),
     ylab = paste("Gene Significance for", trait),
     main = paste("MM vs GS\nModule:", module, "| Trait:", trait),
     col = module,
     pch = 19)
abline(lm(moduleGenes$GS ~ moduleGenes$MM), col = "red", lwd = 2)









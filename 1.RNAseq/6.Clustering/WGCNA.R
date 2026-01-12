args <- commandArgs(trailingOnly = TRUE)
count <- args[1]
pheno <- args[2]
trait <- args[3]
output_dir <- args[4]

suppressMessages(library(dplyr))
library(data.table) %>% suppressMessages()
library(WGCNA) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(reshape) %>% suppressMessages()
library(gplots) %>% suppressMessages()
library(RColorBrewer) %>% suppressMessages()

cat("Trait value received:", trait, "\n")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. Preprocess data
pheno = fread(pheno)
pheno = as.data.frame(pheno)
colnames(pheno)[1] <- "sample_id"
pheno[, -1] <- lapply(pheno[, -1], function(x) as.numeric(as.factor(x)))

count = fread(count)
count = as.data.frame(count)
count <- count[!duplicated(count[, 1]), ]  # Remove duplicated genes
rownames(count) <- count[, 1]              
count <- count[, -1]                
if (any(!pheno$id %in% colnames(count))){
  stop('The sample list in the trait file contains samples that do not exist in the expression matrix.')
}
count = count[,pheno$sample_id]
datExpr0 <- as.data.frame(t(count))         # Transpose

# 2. Check the quality of genes and samples
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    #printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

raw_counts <- as.data.frame(t(datExpr0))    # Transpose back

# 3. Process trait data
# Rename pheno column
num_conditions <- nlevels(as.factor(pheno[[trait]]))
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(as.factor(pheno[[trait]]))]

# 4. Plot the original correlation heatmap
pdf(file.path(output_dir, "adjusted_heatmap.pdf"), width = 12, height = 12)
heatmap.2(cor(raw_counts), RowSideColors = cond_colors, trace = 'none',
          main = 'Sample correlations (raw)', margins = c(12, 12))
while (!is.null(dev.list()))  dev.off()

# 5. Filter low expression genes and log2 transform
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)
raw_counts <- raw_counts[!low_count_mask, ]
log_counts <- log2(raw_counts + 1)

# 6. Density plot
x <- reshape::melt(as.matrix(log_counts))
colnames(x) <- c('gene_id', 'sample', 'value')
p <- ggplot(x, aes(x = value, color = sample)) + 
  geom_density() + labs(title = "Density Plot", x = "Value", y = "Density")
ggsave(file.path(output_dir, "density_plot.pdf"), plot = p, width = 8, height = 6)

# 7. Clustering and detect abnormal samples
sampleTree <- hclust(dist(t(log_counts)), method = "average")
pdf(file.path(output_dir, "sample_clustering.pdf"), width = 6, height = 6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", 
     xlab = "", cex.lab = 1, cex.axis = 1, cex.main = 1)
while (!is.null(dev.list()))  dev.off()

# 8. Extract the intersection of expression and trait information
rownames(pheno) <- pheno$sample_id
index <- intersect(pheno$sample_id, colnames(log_counts))
pheno <- pheno[index, ]
datExpr <- as.data.frame(t(log_counts[, index]))

# 9. Choose the best Soft-Threshold Power
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
candidates <- sft$fitIndices[sft$fitIndices$SFT.R.sq >= 0.8, ]
best_power <- if (nrow(candidates) > 0) {
  candidates$Power[which.max(candidates$SFT.R.sq)]
} else {
  sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
}

# 10. Construct the network and detect modules
cor <- WGCNA::cor
net <- blockwiseModules(datExpr, power = best_power, TOMType = "unsigned", 
                        minModuleSize = 30, networkType = "signed", mergeCutHeight = 0.25,
                        saveTOMs = FALSE, saveTOMFileBase = "yourname", verbose = 3)
cor<-stats::cor
moduleColors <- labels2colors(net$colors)
pdf(file.path(output_dir, "Dendrogram_with_Module_Colors.pdf"), width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
while (!is.null(dev.list()))  dev.off()

# 11. Compute the correlation between modules and traits
#atTraits <- pheno
#conditionLevels <- unique(datTraits$condition)
#datTraits$condition <- factor(datTraits$condition,levels = conditionLevels)
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
design = pheno[-1]
moduleTraitCor <- cor(MEs, design, use = "p")
nSamples <- nrow(datExpr)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# 12. Plot the module-trait correlation heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
pdf(file.path(output_dir, "module_trait_relationships.pdf"), width = 10, height = 6)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(design),
               yLabels = names(MEs),ySymbols=names(MEs),colorLabels=FALSE,
	       textMatrix = textMatrix,cex.lab.y = 0.5, 
               colors = blueWhiteRed(50), cex.text = 0.5, zlim = c(-1, 1),
               main = "Module-trait relationships", margins = c(6, 6))
while (!is.null(dev.list()))  dev.off()

# 13. Extract the highest correlation module genes
modNames = substring(names(MEs), 3)
max_cor_index <- which(abs(moduleTraitCor) == max(abs(moduleTraitCor)), arr.ind = TRUE)[1, ]
best_module <- rownames(moduleTraitCor)[max_cor_index["row"]]
module_name <- gsub("^ME", "", best_module)
genes_in_module <- names(moduleColors)[moduleColors == module_name]
cat("Design matrix column names:", colnames(design), "\n")

if (!(trait %in% colnames(design))) {
  stop(paste("Error: Trait", trait, "not found in the design matrix columns"))
}

mylove <- as.data.frame(design[, trait])
names(mylove) <- "mylove"

geneModuleMembership = as.data.frame(cor(datExpr,   use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

### Gene Significance GS
geneTraitSignificance = as.data.frame(cor(datExpr, mylove$mylove, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(mylove), sep="")
names(GSPvalue) = paste("p.GS.", names(mylove), sep="")
for (module in unique(moduleColors)) {
  
  # Recognize the module column number 
  column = match(module, modNames)
  modNames = substring(names(MEs), 3)
  
  # Extract the module genes
  moduleGenes = moduleColors == module
  
  # Extract the module gene names
  geneNames = colnames(datExpr)[moduleGenes]
  
  # Save the module gene names
  fileName <- paste0(output_dir, "/genes_", module, ".txt")
  write.table(geneNames, file = fileName, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Save
  cat("Saved", length(geneNames), "genes for module", module, "to", fileName, "\n")
  #library(dplyr)
  #library(tibble)
  #library(tidyr)
  #datExpr <- as.data.frame(datExpr)
  #dd <- datExpr[,moduleGenes]
  #dd <- dd %>% 
  #	rownames_to_column("sample") %>% 
  #	mutate(group = datTraits$condition) %>% 
  #	select(group, everything()) %>% 
  #	pivot_longer(cols = 3:ncol(.), names_to = "gene", values_to = "expression")
  
  # Obtain the gene module membership (MM) and gene significance (GS) of this module.
  MM <- geneModuleMembership[moduleGenes, column]
  GS <- geneTraitSignificance[moduleGenes, 1]
  
  # Scatter plot
  pdf(paste0(output_dir, "/scatterplot_", module, ".pdf"), width = 7, height = 7)
  verboseScatterplot(MM, GS,
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for body weight",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  while (!is.null(dev.list()))  dev.off()
  
  cat("Scatterplot for the module with the highest correlation to the trait: ", module, "saved as PDF.\n")
  
}
cat("Done, result in：", output_dir, "\n")

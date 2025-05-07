args <- commandArgs(trailingOnly = TRUE)
count <- args[1]
pheno <- args[2]
trait <- args[3]
output_dir <- args[4]
R_Lib <- args[5]

.libPaths(c(R_Lib, .libPaths()))

# 加载必要的R包
library(data.table)
library(WGCNA)
library(ggplot2)
library(reshape)
library(gplots)
library(RColorBrewer)
# 创建输出目录
cat("Trait value received:", trait, "\n")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. 数据预处理
count = fread(count)
count = as.data.frame(count)
count <- count[!duplicated(count[, 1]), ]  # 去除重复基因
rownames(count) <- count[, 1]               # 第一列作为行名
count <- count[, -1]                        # 去除第一列
datExpr0 <- as.data.frame(t(count))         # 转置数据

pheno = fread(pheno)
pheno = as.data.frame(pheno)

# 2. 检查基因和样本质量
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

raw_counts <- as.data.frame(t(datExpr0))    # 转置回原格式

# 3. 性状数据处理
# 重命名 pheno 的列
colnames(pheno)[1] <- "sample_id"
colnames(pheno)[2] <- "condition"
num_conditions <- nlevels(as.factor(pheno$condition))
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(as.factor(pheno$condition))]

# 4. 绘制原始相关性热图
pdf(file.path(output_dir, "adjusted_heatmap.pdf"), width = 12, height = 12)
heatmap.2(cor(raw_counts), RowSideColors = cond_colors, trace = 'none',
          main = 'Sample correlations (raw)', margins = c(12, 12))
dev.off()

# 5. 低表达基因过滤和 log2 转换
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)
raw_counts <- raw_counts[!low_count_mask, ]
log_counts <- log2(raw_counts + 1)

# 6. 绘制密度图
x <- reshape::melt(as.matrix(log_counts))
colnames(x) <- c('gene_id', 'sample', 'value')
p <- ggplot(x, aes(x = value, color = sample)) + 
  geom_density() + labs(title = "Density Plot", x = "Value", y = "Density")
ggsave(file.path(output_dir, "density_plot.pdf"), plot = p, width = 8, height = 6)

# 7. 聚类并检测异常样本
sampleTree <- hclust(dist(t(log_counts)), method = "average")
pdf(file.path(output_dir, "sample_clustering.pdf"), width = 6, height = 6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", 
     xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# 8. 提取表达量和性状信息的交集
rownames(pheno) <- pheno$sample_id
index <- intersect(pheno$sample_id, colnames(log_counts))
pheno <- pheno[index, ]
datExpr <- as.data.frame(t(log_counts[, index]))

# 9. 选择最佳 Soft-Threshold Power
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
candidates <- sft$fitIndices[sft$fitIndices$SFT.R.sq >= 0.8, ]
best_power <- if (nrow(candidates) > 0) {
  candidates$Power[which.max(candidates$SFT.R.sq)]
} else {
  sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
}

# 10. 构建网络并检测模块
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
dev.off()

# 11. 计算模块与性状的相关性
#atTraits <- pheno
#conditionLevels <- unique(datTraits$condition)
#datTraits$condition <- factor(datTraits$condition,levels = conditionLevels)


MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
design <- model.matrix(~ 0 + pheno$condition)
colnames(design) <- levels(as.factor(pheno$condition))
moduleTraitCor <- cor(MEs, design, use = "p")
nSamples <- nrow(datExpr)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# 12. 绘制模块-性状相关性热图
textMatrix <- paste(signif(moduleTraitCor, 2), "(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
pdf(file.path(output_dir, "module_trait_relationships.pdf"), width = 10, height = 6)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(design),
               yLabels = names(MEs), textMatrix = textMatrix, 
               colors = blueWhiteRed(50), cex.text = 0.5, zlim = c(-1, 1),
               main = "Module-trait relationships", margins = c(6, 6))
dev.off()
cat("部分分析完成！结果保存在：", output_dir, "\n")
# 13. 提取最大相关性的模块和基因

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
  
  # 识别模块的列号
  column = match(module, modNames)
  modNames = substring(names(MEs), 3)
  
  # 提取属于该模块的基因
  moduleGenes = moduleColors == module
  
  # 提取模块的基因名称
  geneNames = colnames(datExpr)[moduleGenes]
  
  # 保存模块的基因名称为 txt 文件
  fileName <- paste0(output_dir, "/genes_", module, ".txt")
  write.table(geneNames, file = fileName, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # 打印保存信息
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
  
  # 获取该模块的基因模块成员关系 (MM) 和基因显著性 (GS)
  MM <- geneModuleMembership[moduleGenes, column]
  GS <- geneTraitSignificance[moduleGenes, 1]
  
  # 保存散点图为 PDF 文件
  pdf(paste0(output_dir, "/scatterplot_", module, ".pdf"), width = 7, height = 7)
  verboseScatterplot(MM, GS,
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for body weight",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  cat("Scatterplot for module", module, "saved as PDF.\n")
  
}
cat("分析完成！结果保存在：", output_dir, "\n")

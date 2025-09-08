rm(list = ls())
suppressMessages(library(gplots))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

LoadInputData <- function(input_seurat) {
  # 检查文件类型并加载数据
  if (grepl("\\.h5$", input_seurat)) {
    scdata <- Read10X_h5(input_seurat)
  } else if (grepl("\\.(txt.gz|csv.gz|tsv.gz|txt|csv|tsv)$", input_seurat)) {
    scdata <- fread(input_seurat, data.table = FALSE)
    rownames(scdata) <- scdata[,1]
    scdata <- scdata[, -1]
  } else if (dir.exists(input_seurat)) {
    scdata <- Read10X(data.dir = input_seurat)
  } else if (grepl("\\.rds$", input_seurat)) {
    scdata <- readRDS(input_seurat)
    if (!inherits(scdata, "Seurat")) {
      stop("The .rds file does not contain a valid Seurat object.")
    }
  } else if (grepl("\\.h5ad$", input_seurat)) {
    stop(paste("Input is an h5ad file, only run in scanpy"))
  } else {
    stop("Input ", input_seurat, " does not meet our requirements")
  }
  # 返回加载的数据
  return(scdata)
}

sys_argv <- commandArgs(trailingOnly = TRUE)
sample1 = sys_argv[1]
sample2 = sys_argv[2]
meta1 = sys_argv[3]
meta2 = sys_argv[4]
clustering1 = trimws(scan(sys_argv[5], what = 'character'))
clustering2 = trimws(scan(sys_argv[6], what = 'character'))
output_directory = sys_argv[7]
nfeatures_used = as.numeric(sys_argv[8]) #2000
integrate_features = as.numeric(sys_argv[9]) #2000
plot_width = as.numeric(sys_argv[10])
plot_height = as.numeric(sys_argv[11])

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
  message("Output directory created: ", output_directory)
}

nfeatures_used = as.numeric(sys_argv[8]) #2000
integrate_features = as.numeric(sys_argv[9]) #2000

#######################################################################################
if (!grepl("\\.rds$", sample1)){
    # 读取数据
    counts <- LoadInputData(sample1)
    # 创建 Seurat 对象
    seurat_obj1 <- CreateSeuratObject(counts = counts)
}else {
    seurat_obj1 <- LoadInputData(sample1)
}

if (!grepl("\\.rds$", sample2)){
    # 读取数据
    counts <- LoadInputData(sample2)
    # 创建 Seurat 对象
    seurat_obj2 <- CreateSeuratObject(counts = counts)
}else {
    seurat_obj2 <- LoadInputData(sample2)
}

if (! meta1 %in% colnames(seurat_obj1@meta.data)){
    stop(paste0(meta1,'not in meta data column of ',sample1))
}
if (! meta2 %in% colnames(seurat_obj2@meta.data)){
    stop(paste0(meta2,'not in meta data column of ',sample2))
}

if (any(! clustering1 %in% unique(seurat_obj1@meta.data[[meta1]]))){
    temp = clustering1[!clustering1 %in% unique(seurat_obj1@meta.data[[meta1]])]
    stop(paste0(temp,'not in ',sample1,' ',meta1))
}
if (any(! clustering2 %in% unique(seurat_obj2@meta.data[[meta2]]))){
    temp = clustering2[!clustering2 %in% unique(seurat_obj2@meta.data[[meta2]])]
    stop(paste0(temp,'not in ',sample2,' ',meta2))
}

DefaultAssay(seurat_obj1) = "RNA"
DefaultAssay(seurat_obj2) = "RNA"

common_gene <- intersect(rownames(seurat_obj1),rownames(seurat_obj2))

seurat_obj1 <- seurat_obj1[common_gene,]
seurat_obj2 <- seurat_obj2[common_gene,]

seurat_obj1 <- NormalizeData(object = seurat_obj1, normalization.method = "LogNormalize")
seurat_obj2 <- NormalizeData(object = seurat_obj2, normalization.method = "LogNormalize")

seurat_obj1 <- FindVariableFeatures(object = seurat_obj1, selection.method = "vst", nfeatures = nfeatures_used)
seurat_obj2 <- FindVariableFeatures(object = seurat_obj2, selection.method = "vst", nfeatures = nfeatures_used)

integr_list <- list(seurat_obj1,seurat_obj2)
integr_features <- SelectIntegrationFeatures(object.list = integr_list, nfeatures = integrate_features)
integr_anchors <- FindIntegrationAnchors(object.list = integr_list, reference = c(1,2), anchor.features = integr_features)
data_integrated <- IntegrateData(anchorset = integr_anchors, normalization.method = "LogNormalize", features.to.integrate = common_gene)

# R
sample_cor <- matrix(NA, nrow = length(clustering1), ncol = length(clustering2), dimnames = list(clustering1, clustering2))

for (cluster1 in clustering1) {
  for (cluster2 in clustering2) {
    cluster_cells1 <- data_integrated@meta.data[data_integrated@meta.data[[meta1]] == cluster1, ]
    cluster_cells2 <- data_integrated@meta.data[data_integrated@meta.data[[meta2]] == cluster2, ]
    
    cluster_avg_expr_1 <- rowMeans(data_integrated@assays$RNA@data[, rownames(cluster_cells1)[!grepl('NA',rownames(cluster_cells1))]])
    cluster_avg_expr_2 <- rowMeans(data_integrated@assays$RNA@data[, rownames(cluster_cells2)[!grepl('NA',rownames(cluster_cells2))]])
    
    correlation_result <- cor(cluster_avg_expr_1, cluster_avg_expr_2, method = "pearson")
    
    sample_cor[as.character(cluster1), as.character(cluster2)] <- correlation_result
  }
}

pdf(paste0(output_directory,"/cor.pdf"),width=plot_width,height=plot_height)
cor_cell <- matrix(as.character(round(sample_cor, 2)), ncol=dim(sample_cor)[2])
ColorRamp <- colorRampPalette(c("white", "#F2FAB0", "red"), bias=1)(100)
heatmap.2(sample_cor, main="Correlation Heatmap", col=ColorRamp, key=TRUE, trace="none", margins=c(10,10), cexRow=1, cexCol=1, cellnote=cor_cell, revC=TRUE)
while (!is.null(dev.list()))  dev.off()
cat('Finished, you can see the results now.\n')

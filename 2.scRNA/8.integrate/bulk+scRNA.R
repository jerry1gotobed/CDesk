suppressMessages(library(ggrepel))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(dplyr))

cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
allcol <- c("#f7acbc","#deab8a","#817936","#444693","#ef5b9c","#fedcbd","#7f7522","#2b4490","#f47920","#80752c","#2a5caa","#f05b72","#905a3d","#87843b","#224b8f","#f15b6c","#8f4b2e","#726930","#003a6c","#f8aba6","#87481f","#454926","#102b6a","#f69c9f","#5f3c23","#2e3a1f","#426ab3","#f58f98","#6b473c","#4d4f36","#46485f","#ca8687","#faa755","#b7ba6b","#4e72b8","#f391a9","#fab27b","#b2d235","#181d4b","#bd6758","#f58220","#5c7a29","#1a2933","#d71345","#843900","#bed742","#121a2a","#d64f44","#905d1d","#7fb80e","#0c212b","#d93a49","#8a5d19","#a3cf62","#6a6da9","#b3424a","#8c531b","#769149","#585eaa","#c76968","#826858","#6d8346","#494e8f","#bb505d","#64492b","#78a355","#afb4db","#987165","#ae6642","#abc88b","#9b95c9","#ac6767","#56452d","#74905d","#6950a1","#973c3f","#96582a","#cde6c7","#6f60aa","#b22c46","#705628","#1d953f","#867892","#a7324a","#4a3113","#77ac98","#918597","#aa363d","#412f1f","#007d65","#6f6d85","#ed1941","#845538","#84bf96","#594c6d","#f26522","#8e7437","#45b97c","#694d9f","#d2553d","#69541b","#225a1f","#6f599c","#b4534b","#d5c59f","#367459","#8552a1","#ef4136","#cd9a5b","#007947","#543044","#c63c26","#cd9a5b","#40835e","#63434f","#f3715c","#b36d41","#2b6447","#7d5886","#a7573b","#df9464","#005831","#401c44","#aa2116","#b76f40","#006c54","#472d56","#b64533","#ad8b3d","#375830","#45224a","#b54334","#dea32c","#274d3d","#411445","#853f04","#d1923f","#375830","#4b2f3d","#840228","#c88400","#27342b","#402e4c","#7a1723","#c37e00","#65c294","#c77eb5","#a03939","#c37e00","#73b9a2","#ea66a6","#8a2e3b","#e0861a","#72baa7","#f173ac","#8e453f","#ffce7b","#005344","#8f4b4a","#fcaf17","#122e29","#892f1b","#ba8448","#293047","#f6f5ec","#6b2c25","#896a45","#00ae9d","#d9d6c3","#733a31","#76624c","#508a88","#d1c7b7","#54211d","#6d5826","#70a19f","#f2eada","#78331e","#ffc20e","#50b7c1","#d3d7d4","#53261f","#fdb933","#00a6ac","#999d9c","#f15a22","#d3c6a6","#78cdd1","#a1a3a6","#b4533c","#c7a252","#008792","#9d9087","#84331f","#dec674","#94d6da","#8a8c8e","#f47a55","#b69968","#afdfe4","#74787c","#f15a22","#c1a173","#5e7c85","#7c8577","#f3704b","#dbce8f","#76becc","#72777b","#da765b","#ffd400","#90d7ec","#77787b","#c85d44","#ffd400","#009ad6","#4f5555","#ae5039","#ffe600","#145b7d","#6c4c49","#6a3427","#f0dc70","#11264f","#563624","#8f4b38","#fcf16e","#7bbfea","#3e4145","#8e3e1f","#decb00","#33a3dc","#3c3645","#f36c21","#cbc547","#228fbd","#464547","#b4532a","#6e6b41","#2468a2","#130c0e","#b7704f","#596032","#2570a1","#281f1d","#de773f","#525f42","#2585a6","#2f271d","#c99979","#5f5d46","#1b315e","#1d1626")
cccol01 <- paste(cccol,"01",sep="")
cccol05 <- paste(cccol,"05",sep="")
cccol30 <- paste(cccol,"30",sep="")
cccol50 <- paste(cccol,"50",sep="")
cccol80 <- paste(cccol,"80",sep="")
options(scipen = 5)
red_col <- colorRampPalette(c(cccol[1],cccol[8],cccol[3],cccol[6]), bias=1)(4)
blue_col <- c(colorRampPalette(c(cccol[4],cccol[7]), bias=1)(5),colorRampPalette(c(cccol[2],cccol[11],cccol[9],cccol[14]), bias=1)(5),colorRampPalette(c(cccol[16],cccol[15]), bias=1)(5))
grey_col <- colorRampPalette(c(cccol[8],"black"), bias=1)(10)
cccol = c(cccol,allcol)

sys_argv <- commandArgs(trailingOnly = TRUE)
bulk_data_files = trimws(scan(sys_argv[1], what = 'character'))
scRNA_data_files = trimws(scan(sys_argv[2], what = 'character'))
meta_file = sys_argv[3]
meta_data <- fread(meta_file)
output_directory = sys_argv[4]
if (!dir.exists(output_directory)) {
  # 如果目录不存在，则创建
  dir.create(output_directory, recursive = TRUE)
  cat("Directory created:", output_directory, "\n")
}
nfeatures_used = as.numeric(sys_argv[5]) #2000
integrate_features = as.numeric(sys_argv[6]) #2000
pc.num = as.numeric(sys_argv[7]) #50
pc_reduce = as.numeric(sys_argv[8]) #30
plot_width = as.numeric(sys_argv[9]) #8
plot_height = as.numeric(sys_argv[10]) #5

########################################################################################
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

# Load scRNA data
scRNA_list <- list()
for(i in seq_along(scRNA_data_files)) {
  if (!grepl("\\.rds$", scRNA_data_files[i])){
    # 读取数据
    counts <- LoadInputData(scRNA_data_files[i])
    # 创建 Seurat 对象
    seurat_obj <- CreateSeuratObject(counts = counts)
  }
  else {
    seurat_obj <- LoadInputData(scRNA_data_files[i])
  }
  # 标准化，取高变基因
  seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize")
  seurat_obj <- FindVariableFeatures(object = seurat_obj, selection.method = "vst", nfeatures = nfeatures_used)
  # 存储到列表中
  scRNA_list[[i]] <- seurat_obj
}

gene_lists <- lapply(scRNA_list, rownames)
scRNA_samples <- lapply(scRNA_list, colnames)
scRNA_samples <- unlist(scRNA_samples)

if (length(intersect(scRNA_samples,meta_data$sample))==0){
  stop('No scRNA data sample found in meta file')
}

common_gene <- Reduce(intersect, gene_lists)
if (length(common_gene) == 0){
  stop('No common gene in scRNA data list')
}

# Load Bulk data
bulk_list = list()
for(i in seq_along(bulk_data_files)) {
  temp <- fread(bulk_data_files[i], header = TRUE, data.table = FALSE) # 将 data.table 转换为 data.frame
  rownames(temp) <- temp[[1]]  # 设置第一列为行名
  temp <- temp[, -1] 
  bulk_list[[i]] <- temp
}
gene_lists <- lapply(bulk_list, rownames)
common_bulk_gene <- Reduce(intersect, gene_lists)
if (length(common_bulk_gene) == 0){
  stop('No common gene in bulkRNA data list')
}

bulk_list = list()
for(i in seq_along(bulk_data_files)) {
  temp <- fread(bulk_data_files[i], header = TRUE, data.table = FALSE) # 将 data.table 转换为 data.frame
  rownames(temp) <- temp[[1]]  # 设置第一列为行名
  temp <- temp[, -1] 
  bulk_list[[i]] <- temp[common_bulk_gene,]
}
bulk_data <- do.call(cbind,bulk_list)

if (any(colnames(bulk_data) %in% meta_data$sample)){
  bulk_data = bulk_data[,colnames(bulk_data) %in% meta_data$sample]
}else{
  stop('No bulk data sample found in meta file')
}

all_genes <- intersect(rownames(bulk_data), common_gene)
if (length(all_genes) == 0){
  stop('No common gene in scRNA and bulkRNA data list')
}

# scRNA对象
scRNA_data <- merge(scRNA_list[[1]],y=scRNA_list[-1])
scRNA_data = NormalizeData(scRNA_data) %>% FindVariableFeatures(nfeatures=nfeatures_used)%>% ScaleData()
scRNA_data = RunPCA(scRNA_data,npcs=pc.num,verbose=FALSE)

scRNA_data  <- IntegrateLayers( 
  object = scRNA_data, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE,features=all_genes
)

# 合并bulk
bulk_data = CreateSeuratObject(counts=bulk_data)

merge_list = c(scRNA_list,list(bulk_data))
merge_data <- merge(merge_list[[1]],y=merge_list[-1])
merge_data = NormalizeData(merge_data) %>% FindVariableFeatures(nfeatures=nfeatures_used)%>% ScaleData()
merge_data = RunPCA(merge_data,npcs=pc.num,verbose=FALSE)

merge_data  <- IntegrateLayers( 
  object = merge_data, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated",
  verbose = FALSE,features=all_genes
)

################# 可视化
meta_data = as.data.frame(meta_data)
rownames(meta_data) = meta_data$sample
mark_cells <- intersect(colnames(scRNA_data),rownames(meta_data))
dataType <- c(rep("scRNA", length(mark_cells)),rep("bulk", ncol(bulk_data)))

merge_data <- FindNeighbors(merge_data, reduction = "integrated", dims = 1:30)
merge_data <- FindClusters(merge_data, resolution = 1, cluster.name = "integrated.harmony")
merge_data <- RunUMAP(merge_data, reduction = "integrated", dims = 1:30, reduction.name = "integrated.umap")
merge_data <- RunTSNE(merge_data, reduction = "integrated", dims = 1:30, reduction.name = "integrated.tsne")

mark_cells <- c(mark_cells,colnames(bulk_data))
select_cells <- subset(merge_data, cells = mark_cells)

celltypes <- rep("NotSure", ncol(select_cells)); names(celltypes) <- colnames(select_cells)
celltypes[intersect(colnames(select_cells),rownames(meta_data))] <- meta_data[intersect(colnames(select_cells),rownames(meta_data)),3]
select_cells@meta.data$celltypes <- celltypes

names(dataType) <- colnames(select_cells)
select_cells@meta.data$dataType <- dataType

pdf(paste0(output_directory,'/scRNA+bulk_integrate.pdf'),width=plot_width,height=plot_height)
highlight_cells <- colnames(bulk_data)

# 创建 DimPlot 图
plot <- DimPlot(select_cells, reduction = "pca", group.by = "celltypes", shape.by = "dataType", 
                cols = cccol[1:length(table(celltypes))], 
                shuffle = TRUE,raster = FALSE)+ scale_shape_manual(values = c(11,20))

# 标记bulk data
pca_coords <- Embeddings(select_cells, "pca")
pca_df <- as.data.frame(pca_coords)
pca_df$cell_id <- colnames(select_cells)
highlight_df <- pca_df[pca_df$cell_id %in% highlight_cells, ]

# 添加标签和虚线连接
plot  +
  geom_text_repel(data = highlight_df, aes(x = PC_1, y = PC_2, label = cell_id),
                  size = 2, color = "black", 
                  segment.color = "black",         # 虚线颜色
                  segment.linetype = "dashed",    # 虚线样式
                  segment.alpha = 0.5,            # 虚线透明度
                  segment.size = 0.2,             # 虚线粗细
                  max.overlaps = Inf,             # 确保所有标签都显示
                  force = 5,                      # 增加标签的分散力度
                  force_pull = 0.5)               # 减少标签的吸附力 
# 创建交互式图形,保存为 HTML 文件
pca_df$datatype = dataType
# 根据 datatype 映射形状
interactive_plot <- plot_ly(pca_df, x = ~PC_1, y = ~PC_2, color = ~celltypes, symbol = ~datatype,
                            text = ~cell_id, hoverinfo = "text") %>%  add_markers()
saveWidget(interactive_plot, paste0(output_directory,"/PCA_integrate.html"), selfcontained = FALSE)

umap_coords <- Embeddings(select_cells, "integrated.umap")
umap_df <- as.data.frame(umap_coords)
umap_df$cell_id <- colnames(select_cells)
highlight_df <- umap_df[umap_df$cell_id %in% highlight_cells, ]
DimPlot(select_cells, reduction = "integrated.umap", group.by = "celltypes", shape.by = "dataType",cols = cccol[1:length(table(celltypes))],shuffle=T,raster = FALSE) + scale_shape_manual(values = c(11,20))+ 
  geom_text_repel(data = highlight_df, aes(x = 	integratedumap_1, y = 	integratedumap_2, label = cell_id),
                  size = 2, color = "black", 
                  segment.color = "black",         # 虚线颜色
                  segment.linetype = "dashed",    # 虚线样式
                  segment.alpha = 0.5,            # 虚线透明度
                  segment.size = 0.2,             # 虚线粗细
                  max.overlaps = Inf,             # 确保所有标签都显示
                  force = 5,                      # 增加标签的分散力度
                  force_pull = 0.5)               # 减少标签的吸附力 
# 创建交互式图形,保存为 HTML 文件
umap_df$datatype = dataType
# 根据 datatype 映射形状
interactive_plot <- plot_ly(umap_df, x = ~	integratedumap_1, y = ~	integratedumap_2, color = ~celltypes, symbol = ~datatype,
                            text = ~cell_id, hoverinfo = "text") %>%  add_markers()
saveWidget(interactive_plot, paste0(output_directory,"/umap_integrate.html"), selfcontained = FALSE)

tsne_coords <- Embeddings(select_cells, "integrated.tsne")
tsne_df <- as.data.frame(tsne_coords)
tsne_df$cell_id <- colnames(select_cells)
highlight_df <- tsne_df[tsne_df$cell_id %in% highlight_cells, ]
DimPlot(select_cells, reduction = "integrated.tsne", group.by = "celltypes",shape.by = "dataType",cols = cccol[1:length(table(celltypes))],shuffle=T,raster = FALSE) + scale_shape_manual(values = c(11,20))+
  geom_text_repel(data = highlight_df, aes(x = tSNE_1, y = tSNE_2, label = cell_id),
                  size = 2, color = "black", 
                  segment.color = "black",         # 虚线颜色
                  segment.linetype = "dashed",    # 虚线样式
                  segment.alpha = 0.5,            # 虚线透明度
                  segment.size = 0.2,             # 虚线粗细
                  max.overlaps = Inf,             # 确保所有标签都显示
                  force = 5,                      # 增加标签的分散力度
                  force_pull = 0.5)               # 减少标签的吸附力 
# 创建交互式图形,保存为 HTML 文件
tsne_df$datatype = dataType
# 根据 datatype 映射形状
interactive_plot <- plot_ly(tsne_df, x = ~tSNE_1, y = ~tSNE_2, color = ~celltypes, symbol = ~datatype,
                            text = ~cell_id, hoverinfo = "text") %>%  add_markers()
saveWidget(interactive_plot, paste0(output_directory,"/tSNE_integrate.html"), selfcontained = FALSE)
while (!is.null(dev.list()))  dev.off()

cat('Finished, you can see the results now.\n')


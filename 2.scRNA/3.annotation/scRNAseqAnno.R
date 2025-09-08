suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))


###################################################################################################
################################            GET PARAMETER           ###############################
###################################################################################################
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
# mode <- 1
print(paste0("mode is ", mode))


###################################################################################################
####################################            FUNCTION           ################################
###################################################################################################
#ref_path <- "/mnt/shiwei/Work/OTHER/CDesk/Fig3/Fig3a/ref_data.rds"
#qry_path <- "/mnt/shiwei/Work/OTHER/CDesk/Fig3/Fig3a/qry_data.rds"
#ref_meta <- "celltype_annotation_20210517"
#save_path <- "/mnt/shiwei/Work/OTHER/CDesk/Fig3/Fig3a"
# TransferLabelAnnotation(ref_data, qry_data, "celltype_annotation_20210517", "/mnt/shiwei/Work/OTHER/CDesk/Fig3")
TransferLabelAnnotation <- function(ref_path, qry_path, ref_meta, save_path,nfeatures_used,dims_used){
  # read data
  ref_data <- readRDS(ref_path)
  qry_data <- readRDS(qry_path)
  if (! ref_meta %in% colnames(ref_data@meta.data)){
    stop(paste0(ref_meta,'not in meta.data column'))
  } 
  # FindTransferAnchors & TransferData
  common_features <- intersect(rownames(ref_data), rownames(qry_data))
  ref_data <- subset(ref_data, features = common_features)
  qry_data <- subset(qry_data, features = common_features)
  
  ref_data <- NormalizeData(ref_data, normalization.method = "LogNormalize")
  ref_data <- FindVariableFeatures(ref_data, selection.method = "vst", nfeatures = nfeatures_used)
  ref_data <- ScaleData(ref_data)
  ref_data <- RunPCA(ref_data, npcs = 50)
  
  qry_data <- NormalizeData(qry_data, normalization.method = "LogNormalize")
  qry_data <- FindVariableFeatures(qry_data, selection.method = "vst", nfeatures = nfeatures_used)
  qry_data <- ScaleData(qry_data)
  qry_data <- RunPCA(qry_data, npcs = 50)
  
  anchors <- FindTransferAnchors(
    reference = ref_data,
    query = qry_data,
    dims = 1:dims_used,
    reduction = "pcaproject",
    reference.reduction = "pca",
    features = common_features
  )
  
  # Transfer cell type annotations
  TransferCelltype <- TransferData(
    anchorset = anchors,
    refdata = ref_data@meta.data[[ref_meta]],
    dims = 1:dims_used
  )
  
  # Add predicted cell type and embryonic period to metadata
  qry_data@meta.data$PredCelltype <- TransferCelltype$predicted.id
  
  result_df <- data.frame(
    Cell=rownames(qry_data@meta.data),
    PredCelltype=qry_data@meta.data$PredCelltype
  )
  
  write.table(result_df, paste0(save_path, "/PredCelltype.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

}

#data_path <- "/mnt/shiwei/Work/OTHER/CDesk/Fig3/Fig3a/qry_data.rds"
#meta_info_path <- "/mnt/shiwei/Work/OTHER/CDesk/Fig3/Fig3a/marker_info.xlsx"
#meta_col <- "umap_clusters"
#save_path <- "/mnt/shiwei/Work/OTHER/CDesk/Fig3/Fig3a"
#expression_thre <- 1.5
#percentage_thre <- 0.7
#marker_percentage <- 0.7
# MarkerAnnotation(data, meta_info, "umap_clusters", save_path, 1.5, 0.8, 0.8)
MarkerAnnotation <- function(data_path, meta_info_path, meta_col, save_path, expression_thre, percentage_thre, marker_percentage){
  # read data
  meta_info <- fread(meta_info_path)
  data <- readRDS(data_path)
  if (! meta_col %in% colnames(data@meta.data)){
    stop(paste0(meta_col,'not in meta.data column'))
  } 
  Idents(data) <- data@meta.data[[meta_col]]
  
  common_genes <- intersect(rownames(data),meta_info$Marker)
  data <- data[common_genes,]
  meta_info <- meta_info[meta_info$Marker %in% common_genes,]
  
  a <- table(meta_info$Celltype,meta_info$Marker)
  type_to_remove <- rownames(a)[rowSums(a > 0) == 1]
  meta_info <- meta_info[!meta_info$Celltype %in% type_to_remove, ]
  
  result_df <- as.data.frame(matrix(0, nrow = length(unique(meta_info$Celltype)), ncol = length(unique(data@meta.data[[meta_col]]))))
  rownames(result_df) <- unique(meta_info$Celltype)
  colnames(result_df) <- unique(data@meta.data[[meta_col]])
  
  # 得到每个cluster的celltype
  for (type in unique(meta_info$Celltype)) {
    marker <- meta_info[meta_info$Celltype==type,]$Marker
    avg_exp <- AverageExpression(
      data,
      features = marker,
      assays = "RNA",  # 默认是"RNA"，如果是其他assay请修改
      slot = "data"    # 使用标准化数据（log1p），如果用原始counts则用"counts"
    )
    
    # 计算平均表达量
    avg_exp_df <- as.data.frame(avg_exp$RNA)
    
    # 计算表达大于0的细胞比例
    pct_exp_df <- as.data.frame(matrix(0, nrow = nrow(avg_exp_df), ncol = ncol(avg_exp_df)))
    rownames(pct_exp_df) <- rownames(avg_exp_df)
    colnames(pct_exp_df) <- colnames(avg_exp_df)
    
    for (cluster in colnames(pct_exp_df)) {
      cluster_cells <- WhichCells(data, idents = cluster)
      expr_matrix <- GetAssayData(data, assay = "RNA", slot = "counts")[marker, cluster_cells]
      moreZeroCount <- rowSums(expr_matrix > 0)
      pct_exp_df[,cluster] <- moreZeroCount / length(cluster_cells)
    }
    
    avg_exp_above <- colSums(avg_exp_df > expression_thre) / length(marker)
    pct_exp_above <- colSums(pct_exp_df > percentage_thre) / length(marker)
    

    for (cluster in colnames(pct_exp_df)) {
      if (avg_exp_above[cluster] >= marker_percentage & pct_exp_above[cluster] >= marker_percentage) {
        result_df[type,cluster] = (avg_exp_above[cluster]+pct_exp_above[cluster])/2
      }
    }

  }
  
  # 给每个细胞分配type并且保存为txt
  cell_file <- data.frame(Cell = character(), PredCelltype = character())
  
  result_df <- result_df[, colSums(result_df != 0) > 0]
  for (i in colnames(result_df)) {
    temp_cell <- colnames(data[,data@meta.data[[meta_col]]==i])
    temp_file <- data.frame(
      Cell=temp_cell,
      PredCelltype=rownames(result_df)[which.max(result_df[, i])]
    )
    cell_file <- rbind(cell_file,temp_file)
  }
  
  write.table(cell_file, paste0(save_path, "/PredCelltype.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  
}

#data_path <- "/mnt/shiwei/Work/OTHER/CDesk/Fig3/Fig3a/qry_data.rds"
#save_path <- "/mnt/shiwei/Work/OTHER/CDesk/Fig3/Fig3c"
#cluster_type <- "tsne"
PlotData <- function(data_path, save_path, cluster_type,mode){
  data <- readRDS(data_path)
  
  PredCelltype_df <- read.table(paste0(save_path,"/PredCelltype.txt"), sep = "\t", header = T)
  data@meta.data$PredType <- ifelse(rownames(data@meta.data) %in% PredCelltype_df$Cell,PredCelltype_df$PredCelltype[match(rownames(data@meta.data), PredCelltype_df$Cell)],"unassigned")
  
  # samples
  p <- DimPlot(data,reduction=cluster_type,label=T,group.by = "PredType")
  ggsave(paste0(save_path,"/PredCelltype.",cluster_type,".pdf"),p,width = 7,height = 6)	
}

###################################################################################################
#####################################            RUNNING           ################################
###################################################################################################
if (mode=='marker') {
  data_path <- args[2]
  print(paste0("data_path is ", data_path))
  meta_info_path <- args[3]
  print(paste0("meta_info_path is ", meta_info_path))
  meta_col <- args[4]
  print(paste0("meta_col is ", meta_col))
  save_path <- args[5]
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
    message("Output directory created: ", save_path)
  }
  print(paste0("save_path is ", save_path))
  expression_thre <- as.numeric(args[6])
  print(paste0("expression_thre is ", expression_thre))
  percentage_thre <- as.numeric(args[7])
  print(paste0("percentage_thre is ", percentage_thre))
  marker_percentage <- as.numeric(args[8])
  print(paste0("marker_percentage is ", marker_percentage))
  cluster_type <- args[9]
  print(paste0("cluster_type is ", cluster_type))
  
  MarkerAnnotation(data_path, meta_info_path, meta_col, save_path, expression_thre, percentage_thre, marker_percentage)
  PlotData(data_path, save_path, cluster_type,mode)

} else if (mode=='transfer') {
  ref_path <- args[2]
  print(paste0("ref_path is ", ref_path))
  qry_path <- args[3]
  print(paste0("qry_path is ", qry_path))
  ref_meta <- args[4]
  print(paste0("ref_meta is ", ref_meta))
  save_path <- args[5]
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
    message("Output directory created: ", save_path)
  }
  print(paste0("save_path is ", save_path))
  cluster_type <- args[6]
  print(paste0("cluster_type is ", cluster_type))
  nfeatures_used <- as.numeric(args[7])
  print(paste0("nfeatures : ", nfeatures_used))
  dims_used <- as.numeric(args[8])
  print(paste0("dims : ", dims_used))
  TransferLabelAnnotation(ref_path, qry_path, ref_meta, save_path,nfeatures_used,dims_used)
  PlotData(qry_path, save_path, cluster_type,mode)
}

cat('Finished! You can check the result now\n')

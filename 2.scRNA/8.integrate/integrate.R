rm(list = ls())
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(harmony))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

sys_argv <- commandArgs(trailingOnly = TRUE)
###读取并处理数据
###读取文件名称
input_file = sys_argv[1]
mode = sys_argv[2]
output_directory = sys_argv[3]
if (!dir.exists(output_directory)) {
  # 如果目录不存在，则创建
  dir.create(output_directory, recursive = TRUE)
  cat("Directory created:", output_directory, "\n")
}

min_cells = as.numeric(sys_argv[4]) #3
min_features = as.numeric(sys_argv[5]) #200
nfeatures_used = as.numeric(sys_argv[6]) #2000
integrate_features = as.numeric(sys_argv[7]) #2000
pc.num = as.numeric(sys_argv[8]) #50
umap_pc = as.numeric(sys_argv[9]) #30
max_iter_harmony = as.numeric(sys_argv[10]) #20
plot_width = as.numeric(sys_argv[11])
plot_height = as.numeric(sys_argv[12])

group = fread(input_file, data.table = FALSE)
input_files = unlist(group['file'])
samples = unlist(group['group'])

# 显示输入文件列表
cat("Input files:", input_files, "\n")
# 显示样本名称
cat("Samples:", samples, "\n")
# 显示模式
cat("Mode:", mode, "\n")
# 显示输出目录
cat("Output Directory:", output_directory, "\n")
# 显示过滤参数
cat("Minimum cells:", min_cells, "\n")
cat("Minimum features:", min_features, "\n")
# 显示用于整合和分析的特征数
cat("Number of features used:", nfeatures_used, "\n")
cat("Number of features for integration:", integrate_features, "\n")
# 显示 PCA 和 UMAP 参数
cat("Number of PCs for PCA:", pc.num, "\n")
cat("Number of PCs for UMAP:", umap_pc, "\n")
# 显示 Harmony 的最大迭代次数
cat("Maximum iterations for Harmony:", max_iter_harmony, "\n")

# anchor,merge,harmony
# input_files = c('/mnt/linzejie/temp/jerry/data/GSM5050521_G1counts.csv.gz','/mnt/linzejie/temp/jerry/data/GSM5050523_G2counts.csv.gz')
# samples = c('a','b')
# mode = 'harmony'
# output_directory = '/mnt/linzejie/test/integrate'

# 定义一个函数来加载不同格式的输入文件
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


Seu_obj_list <- list()
if (mode=='merge'){
  for(i in seq_along(input_files)) {
    # 提取样本名称
    sample_name <- samples[i]
    # 读取数据
    #counts <- LoadInputData(input_files[i])
    # 创建 Seurat 对象，添加样本名称作为前缀
    #seurat_obj <- CreateSeuratObject(counts = counts,  project = sample_name,min.cells = min_cells, min.features = min_features)
    if (!grepl("\\.rds$", input_files[i])){
    	# 读取数据
    	counts <- LoadInputData(input_files[i])
    	# 创建 Seurat 对象
    	seurat_obj <- CreateSeuratObject(counts = counts,  project = sample_name,min.cells = min_cells, min.features = min_features)
    }
    else {
	seurat_obj <- LoadInputData(input_files[i])
    }
    # 用样本名称作为前缀，确保细胞名称唯一
    seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_name)
    seurat_obj@meta.data$merge_sample = sample_name
    # 存储到列表中
    Seu_obj_list[[sample_name]] <- seurat_obj
  }
  sce.all<-merge(Seu_obj_list[[1]],y=Seu_obj_list[-1])
  saveRDS(sce.all,file=paste0(output_directory,'/','merged.rds'))
  sce.all = NormalizeData(sce.all) %>% FindVariableFeatures(nfeatures=nfeatures_used)%>% ScaleData()
  sce.all = RunPCA(sce.all,npcs=pc.num,verbose=FALSE)
  sce.all = RunUMAP(sce.all,dims=1:umap_pc)
  p = DimPlot(sce.all,group.by = 'merge_sample')
  ggsave(paste0(output_directory,'/','batch_check_merge.pdf'),width=plot_width,height=plot_height,p)
} else if(mode=='CCA'){
  for(i in seq_along(input_files)) {
    sample_name <- samples[i]
    if (!grepl("\\.rds$", input_files[i])){
        # 读取数据
        counts <- LoadInputData(input_files[i])
        # 创建 Seurat 对象
        seurat_obj <- CreateSeuratObject(counts = counts,  project = sample_name,min.cells = min_cells, min.features = min_features)
    }
    else {
        seurat_obj <- LoadInputData(input_files[i])
    }
    # 标准化
    seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize")
    seurat_obj <- FindVariableFeatures(object = seurat_obj, selection.method = "vst", nfeatures = nfeatures_used)
    # 用样本名称作为前缀，确保细胞名称唯一
    seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_name)
    seurat_obj@meta.data$merge_sample = sample_name
    # 存储到列表中
    Seu_obj_list[[sample_name]] <- seurat_obj
  }
  gene_lists <- lapply(Seu_obj_list, rownames)
  common_gene <- Reduce(intersect, gene_lists)
  Integr_features <- SelectIntegrationFeatures(object.list = Seu_obj_list, nfeatures = integrate_features)
  Integr_anchors <- FindIntegrationAnchors(object.list = Seu_obj_list, anchor.features = Integr_features)
  sce.all <- IntegrateData(anchorset = Integr_anchors, normalization.method = "LogNormalize", features.to.integrate = common_gene)
  saveRDS(sce.all,file=paste0(output_directory,'/','anchor_integrate.rds'))
  DefaultAssay(sce.all) <- "integrated"
  sce.all <- ScaleData(sce.all, verbose = FALSE)
  sce.all = RunPCA(sce.all,npcs=pc.num,verbose=FALSE)
  sce.all = RunUMAP(sce.all,dims=1:umap_pc)
  p = DimPlot(sce.all,group.by = 'merge_sample')
  ggsave(paste0(output_directory,'/','batch_check_anchor.pdf'),width=plot_width,height=plot_height,p)
} else if(mode=='harmony'){
  for(i in seq_along(input_files)) {
    sample_name <- samples[i]
    if (!grepl("\\.rds$", input_files[i])){
        # 读取数据
        counts <- LoadInputData(input_files[i])
        # 创建 Seurat 对象
        seurat_obj <- CreateSeuratObject(counts = counts,  project = sample_name,min.cells = min_cells, min.features = min_features)
    }
    else {
        seurat_obj <- LoadInputData(input_files[i])
    }
    # 用样本名称作为前缀，确保细胞名称唯一
    seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_name)
    seurat_obj@meta.data$merge_sample = sample_name
    # 存储到列表中
    Seu_obj_list[[sample_name]] <- seurat_obj
  }
  sce.all<-merge(Seu_obj_list[[1]],y=Seu_obj_list[-1])
  sce.all = NormalizeData(sce.all) %>% FindVariableFeatures(nfeatures=nfeatures_used)%>% ScaleData()
  sce.all = RunPCA(sce.all,npcs=pc.num,verbose=FALSE)
  sce.all <- RunHarmony(sce.all, group.by.vars="merge_sample", max.iter.harmony=max_iter_harmony)
  saveRDS(sce.all,file=paste0(output_directory,'/','merged_harmony.rds'))
  sce.all = RunUMAP(sce.all,reduction = 'harmony',dims=1:umap_pc)
  p = DimPlot(sce.all,group.by = 'merge_sample')
  ggsave(paste0(output_directory,'/','batch_check_harmony.pdf'),width=plot_width,height=plot_height,p)
}

cat('Finished, you can see the results now.\n')

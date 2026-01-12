rm(list = ls())
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(harmony))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

sys_argv <- commandArgs(trailingOnly = TRUE)

input_file = sys_argv[1]
mode = sys_argv[2]
output_directory = sys_argv[3]
if (!dir.exists(output_directory)) {
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
k_anchor = as.numeric(sys_argv[11])
plot_width = as.numeric(sys_argv[12])
plot_height = as.numeric(sys_argv[13])

group = fread(input_file, data.table = FALSE)
input_files = unlist(group['file'])
samples = unlist(group['group'])

cat("Input files:", input_files, "\n")
cat("Samples:", samples, "\n")
cat("Mode:", mode, "\n")
cat("Output Directory:", output_directory, "\n")
cat("Minimum cells:", min_cells, "\n")
cat("Minimum features:", min_features, "\n")
cat("Number of features used:", nfeatures_used, "\n")
cat("Number of features for integration:", integrate_features, "\n")
cat("Number of PCs for PCA:", pc.num, "\n")
cat("Number of PCs for UMAP:", umap_pc, "\n")
cat("Maximum iterations for Harmony:", max_iter_harmony, "\n")

# Read input data
LoadInputData <- function(input_seurat) {
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
  return(scdata)
}


Seu_obj_list <- list()
if (mode=='merge'){
  for(i in seq_along(input_files)) {
    # Read the sample name
    sample_name <- samples[i]
    if (!grepl("\\.rds$", input_files[i])){
    	# Read data
    	counts <- LoadInputData(input_files[i])
    	# Create Seurat Object
    	seurat_obj <- CreateSeuratObject(counts = counts,  project = sample_name,min.cells = min_cells, min.features = min_features)
    }
    else {
	seurat_obj <- LoadInputData(input_files[i])
    }
    # Sample name as prefix
    seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_name)
    seurat_obj@meta.data$merge_sample = sample_name
    # Save in the list
    Seu_obj_list[[sample_name]] <- seurat_obj
  }
  sce.all<-merge(Seu_obj_list[[1]],y=Seu_obj_list[-1])
  saveRDS(sce.all,file=paste0(output_directory,'/','merged.rds'))
  sce.all = NormalizeData(sce.all,verbose=F) %>% FindVariableFeatures(nfeatures=nfeatures_used,verbose=F)%>% ScaleData(verbose=F)
  sce.all = RunPCA(sce.all,npcs=pc.num,verbose=FALSE)
  sce.all = RunUMAP(sce.all,dims=1:umap_pc,verbose=F)
  p = DimPlot(sce.all,group.by = 'merge_sample')
  ggsave(paste0(output_directory,'/','batch_check_merge.pdf'),width=plot_width,height=plot_height,p)
} else if(mode=='CCA'){
  for(i in seq_along(input_files)) {
    sample_name <- samples[i]
    if (!grepl("\\.rds$", input_files[i])){
        counts <- LoadInputData(input_files[i])
        seurat_obj <- CreateSeuratObject(counts = counts,  project = sample_name,min.cells = min_cells, min.features = min_features)
    }
    else {
        seurat_obj <- LoadInputData(input_files[i])
    }
    seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_name)
    seurat_obj@meta.data$merge_sample = sample_name
    Seu_obj_list[[sample_name]] <- seurat_obj
  }
  sce.all<-merge(Seu_obj_list[[1]],y=Seu_obj_list[-1])
  sce.all = NormalizeData(sce.all,verbose=F) %>% FindVariableFeatures(nfeatures=nfeatures_used,verbose=F)%>% ScaleData(verbose=F)
  sce.all = RunPCA(sce.all,npcs=pc.num,verbose=FALSE)
  sce.all  <- IntegrateLayers(
    object = sce.all, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "CCA",
    #group.by.vars    = "merge_sample",
    verbose = FALSE, k.anchor = k_anchor
  )
  saveRDS(sce.all,file=paste0(output_directory,'/','merged_CCA.rds'))
  sce.all = RunUMAP(sce.all,reduction = 'CCA',dims=1:umap_pc,verbose=F)
  p = DimPlot(sce.all,group.by = 'merge_sample')
  ggsave(paste0(output_directory,'/','batch_check_CCA.pdf'),width=plot_width,height=plot_height,p)
} else if(mode=='harmony'){
  for(i in seq_along(input_files)) {
    sample_name <- samples[i]
    if (!grepl("\\.rds$", input_files[i])){
        counts <- LoadInputData(input_files[i])
        seurat_obj <- CreateSeuratObject(counts = counts,  project = sample_name,min.cells = min_cells, min.features = min_features)
    }
    else {
        seurat_obj <- LoadInputData(input_files[i])
    }
    seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_name)
    seurat_obj@meta.data$merge_sample = sample_name
    Seu_obj_list[[sample_name]] <- seurat_obj
  }
  sce.all<-merge(Seu_obj_list[[1]],y=Seu_obj_list[-1])
  sce.all = NormalizeData(sce.all,verbose=F) %>% FindVariableFeatures(nfeatures=nfeatures_used,verbose=F)%>% ScaleData(verbose=F)
  sce.all = RunPCA(sce.all,npcs=pc.num,verbose=FALSE)
  sce.all <- IntegrateLayers(
    object = sce.all, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    group.by.vars    = "merge_sample",
    verbose = FALSE,max.iter.harmony = max_iter_harmony
  )
  saveRDS(sce.all,file=paste0(output_directory,'/','merged_harmony.rds'))
  sce.all = RunUMAP(sce.all,reduction = 'harmony',dims=1:umap_pc,verbose=F)
  p = DimPlot(sce.all,group.by = 'merge_sample')
  ggsave(paste0(output_directory,'/','batch_check_harmony.pdf'),width=plot_width,height=plot_height,p)
}

cat('Finished, you can see the results now.\n')


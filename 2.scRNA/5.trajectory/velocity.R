library(Seurat)

sys_argv <- commandArgs(T)
input_seurat = sys_argv[1]
output_directory =sys_argv[2]
meta = sys_argv[3]

if (length(sys_argv) != 3) {
  stop("Input, output and meta arguments should be provided.")
}

data = readRDS(input_seurat)

if (! meta %in% colnames(data@meta.data)){
  stop(paste0(meta,'not in meta data column of ',input_seurat))
}

## cell ID
write.csv(row.names(data@meta.data),file=paste0(output_directory,"/cellID_obs.csv"),row.names=FALSE)
## UMAP
write.csv(Embeddings(data,reduction="umap"),file=paste0(output_directory,"/cell_embeddings.csv"))
##seurat clustering
write.csv(data@meta.data[meta],paste0(output_directory,"/clusters.csv"))

suppressMessages(library(Seurat))
suppressMessages(library(monocle3))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
rm(list=ls())

# 输入聚类并且注释好的rds对象，轨迹分析的metadata（如celltype或其他），轨迹起点类型
sys_argv <- commandArgs(T)
input_seurat = sys_argv[1] # 
output_directory =sys_argv[2] # 
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
  message("Output directory created: ", output_directory)
}
meta = sys_argv[3] # celltype
start_p = sys_argv[4] # NK
plot_width = as.numeric(sys_argv[5])
plot_height = as.numeric(sys_argv[6])

if (length(sys_argv) != 6) {
  stop("6 command line arguments should be provided.")
}

# Import seurat object 
seu_obj <- readRDS(input_seurat) # "Seurat_single_sample_scobj.rds"

# Check meta and start_p exists or not
if (! meta %in% colnames(seu_obj@meta.data)){
  stop(paste0(meta,'not in meta data column of ',input_seurat))
}

if (! start_p %in% unique(seu_obj@meta.data[[meta]])){
  stop(paste0(start_p,'not in ',meta))
}

# Get seurat expression matrix
data <- GetAssayData(seu_obj, assay = 'RNA', layer = 'counts')

# Get meta information
cell_metadata <- seu_obj@meta.data

# Get gene annotation
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Make CDS object
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)     # Like seurat:NormalizeData+ScaleData+RunPCA

# umap
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method default: PCA

# tSNE
cds <- reduce_dimension(cds, reduction_method="tSNE")

# clustering
cds <- cluster_cells(cds) 

# Replace umap embedding
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seu_obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  

# Learn trajectory
cds <- learn_graph(cds) 

# Plot trajectory
p = plot_cells(cds,
           color_cells_by = meta,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size=4,
           cell_size=1.5)
ggsave(paste0(output_directory,'/Pseudotime_in',meta,'.pdf'), plot = p, width = plot_width, height = plot_height)

# a helper function to identify the root principal points:
get_earliest_principal_node  <- function(cds, time_bin=start_p,meta){
  cell_ids <- which(colData(cds)[, meta] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds,start_p,meta))

p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
ggsave(paste0(output_directory,'/Pseudotime.pdf'), plot = p, width = plot_width, height = plot_height)

# DEGs
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3) # ***

# Top10 DEGs
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Top DEGs expression trend
p = plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by=meta, 
                         min_expr=0.5, ncol = 2)
ggsave(paste0(output_directory,"/Top10_DEGs_expression_trend.pdf"), plot = p, width = plot_width, height = plot_height)
# FeaturePlot
p = plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
ggsave(paste0(output_directory,"/Top10_DEGs_FeaturePlot.pdf"), plot = p, width = plot_width, height = plot_height)

## Coexpressin module
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10) # ***
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)[[meta]])
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2", 
                   filename = paste0(output_directory, "/Coexpression_DEGs_Module.pdf"),width = plot_width, height = plot_height)

# Output pseudotime seurat object
pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(seu_obj@meta.data)]
pseudotime[!is.finite(pseudotime)] <- 0
seu_obj$pseudotime <- pseudotime
p = FeaturePlot(seu_obj, reduction = "umap", features = "pseudotime")
ggsave(paste0(output_directory,"/Pseudotime_Seurat.pdf"), plot = p, width = plot_width, height = plot_height)
saveRDS(seu_obj, file = paste0(output_directory,"/sco_pseudotime.rds"))
write.csv(gene_module,file = paste0(output_directory,'/gene_module.csv'))
write.csv(gene_module,file = paste0(output_directory,'/gene_module.csv'))
#while (!is.null(dev.list())) dev.off()
cat('Finished, you can check the results now.\n')

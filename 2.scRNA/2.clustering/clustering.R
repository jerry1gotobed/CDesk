rm(list=ls())
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(clustree))
suppressMessages(library(rhdf5))
suppressMessages(library(data.table))

# Command line arguments
sys_argv <- commandArgs(trailingOnly = TRUE)  # 使用 trailingOnly = TRUE 来跳过脚本路径

# Check that at least two arguments are provided (input and output directories)
if (length(sys_argv) < 2) {
  stop("Not enough command line arguments provided. At least input and output directory are required.")
}

input_seurat <- sys_argv[1]  # 第一个参数：输入文件
output_directory <- sys_argv[2]  # 第二个参数：输出目录
sample_name <- sys_argv[3]

# Set default values for optional parameters
min_cells <- ifelse(length(sys_argv) >= 4, as.numeric(sys_argv[4]), 3)  # Default: 3
min_features <- ifelse(length(sys_argv) >= 5, as.numeric(sys_argv[5]), 200)  # Default: 200
nFeature_RNA_min <- ifelse(length(sys_argv) >= 6, as.numeric(sys_argv[6]), 200)  # Default: 200
nFeature_RNA_max <- ifelse(length(sys_argv) >= 7, as.numeric(sys_argv[7]), 2500)  # Default: 2500
mt_percent <- ifelse(length(sys_argv) >= 8, as.numeric(sys_argv[8]), 5)  # Default: 5
variable_features <- ifelse(length(sys_argv) >= 9, as.numeric(sys_argv[9]), 2000)  # Default: 2000
dim_prefer <- ifelse(length(sys_argv) >= 10, as.numeric(sys_argv[10]), 30)  # Default: 30
res <- ifelse(length(sys_argv) >= 11, as.numeric(sys_argv[11]), 1)  # Default: 1
width = ifelse(length(sys_argv) >= 12, as.numeric(sys_argv[12]), 10)  # Default: 10
height = ifelse(length(sys_argv) >= 13, as.numeric(sys_argv[13]), 8)  # Default: 8

# Ensure required directories exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
  message("Output directory created: ", output_directory)
}

# Read input data
cat("Reading input data...\n")
if (grepl("\\.h5$", input_seurat)) {
  scdata <- Read10X_h5(input_seurat)
} else if (grepl("\\.(txt.gz|csv.gz|tsv.gz)$", input_seurat)) {
  scdata <- fread(input_seurat, data.table = FALSE)
  rownames(scdata) <- scdata$V1
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

cat("Input Seurat file: ", input_seurat, "\n")
cat("Output directory: ", output_directory, "\n")
pdf(paste0(output_directory,"/",sample_name,"_Seurat_clustering.pdf"),width=width,height=height)

# Create Seurat object
cat("Creating Seurat object...\n")
cat("Min cells: ", min_cells, "\n")
cat("Min features: ", min_features, "\n")
if (!grepl("\\.rds$", input_seurat)){
	seu_obj <- CreateSeuratObject(counts = scdata, 
                              project = "seurat", 
                              min.cells = min_cells, 
                              min.features = min_features)
} else {
	seu_obj = scdata
}

# Add percentage of mitochondrial genes
#cat("Calculating mitochondrial gene percentage...\n")
seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")

# Filter cells based on thresholds
cat("nFeature_RNA min: ", nFeature_RNA_min, "\n")
cat("nFeature_RNA max: ", nFeature_RNA_max, "\n")
cat("MT percentage threshold: ", mt_percent, "\n")
cat("Subsetting Seurat object based on thresholds...\n")
seu_obj <- subset(seu_obj, 
                  subset = nFeature_RNA > nFeature_RNA_min & 
                    nFeature_RNA < nFeature_RNA_max & 
                    percent.mt < mt_percent)

# Normalize, scale data and find variable features
cat("Normalizing and scaling data...\n")
seu_obj <- NormalizeData(seu_obj,verbose=F)
cat("Variable features: ", variable_features, "\n")
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = variable_features,verbose=F)

# Scale the data
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes,verbose=F)

# Perform PCA
cat("Running PCA...\n")
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj),verbose=F)

# Elbow plot to determine the number of dimensions to use
cat("Generating Elbow plot...\n")
ElbowPlot(seu_obj,ndims = 50)

# Find optimal number of dimensions for clustering
cat("Use", dim_prefer, "dimensions for clustering.\n")

# Find neighbors and clusters
cat("Finding neighbors...\n")
seu_obj <- FindNeighbors(seu_obj, dims = 1:dim_prefer,verbose=F)
seu_obj <- FindClusters(seu_obj, resolution = res,verbose=F)

# Run UMAP
cat("Running UMAP...\n")
seu_obj <- RunUMAP(seu_obj, dims = 1:dim_prefer,verbose=F)

# Visualize UMAP
cat("Generating UMAP plot...\n")
umap_plot <- DimPlot(seu_obj, reduction = "umap", label = TRUE)
print(umap_plot)

# Run t-SNE
cat("Running t-SNE...\n")
seu_obj <- RunTSNE(seu_obj, dims = 1:dim_prefer,verbose=F)

# Visualize t-SNE
cat("Generating t-SNE plot...\n")
tsne_plot <- DimPlot(seu_obj, reduction = "tsne", label = TRUE)
print(tsne_plot)

# Visualize cluster tree
seu_obj <- FindClusters(seu_obj, resolution = seq(0.2, 1.2, 0.1),verbose=F)
cat("Generating cluster tree...\n")
clustree(seu_obj) + ggtitle(paste('Use', dim_prefer, 'dimensions for dimensionality reduction'))

# Save plots to output PDF
while (!is.null(dev.list()))  dev.off()
saveRDS(seu_obj,file=paste0(output_directory,'/',sample_name,'_seurat_clustering.rds'))
cat('Finished! You can check the result now.\n')

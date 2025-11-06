suppressMessages(library(destiny))
suppressMessages(library(Seurat))
suppressMessages(library(ggthemes))
suppressMessages(library(ggbeeswarm))
rm(list=ls())

sys_argv <- commandArgs(T)
input_seurat = sys_argv[1] # 
output_directory =sys_argv[2] # 
meta = sys_argv[3] # celltype
start_p = sys_argv[4] # NK
plot_width = as.numeric(sys_argv[5])
plot_height = as.numeric(sys_argv[6])

if (length(sys_argv) != 6) {
  stop("6 command line arguments should be provided.")
}

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
  message("Output directory created: ", output_directory)
}

seurat_obj <- readRDS(input_seurat)

if (! meta %in% colnames(seurat_obj@meta.data)){
  stop(paste0(meta,'not in meta data column of ',input_seurat))
}

if (! start_p %in% unique(seurat_obj@meta.data[[meta]])){
  stop(paste0(start_p,'not in ',meta))
}

# 获取 Seurat 对象中的表达矩阵（基因表达数据）
data.matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
data.matrix = t(as.matrix(data.matrix))

pdf(paste0(output_directory,'/DiffusionMap.pdf'),width=plot_width,height=plot_height)
sigs <- find_sigmas(data.matrix, verbose = TRUE)
dm = DiffusionMap(data.matrix, sigs)
plot(dm)
plot(dm,c(1,2))
plot(dm,c(1,3))
plot(dm,c(2,3))

# 指定起始细胞
tips = rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[meta]] == start_p][1]
tips_index <- which(rownames(dm@eigenvectors)==tips)

# 计算伪时间
dpt <- DPT(dm, tips = tips_index)
plot(dpt)

# 
tmp <- data.frame(DC1 = dm$DC1,DC2 = dm$DC2,Timepoint = seurat_obj@meta.data[[meta]],dpt = dpt$DPT1)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +  
  geom_point() + scale_color_tableau() +   
  xlab("Diffusion component 1") +   
  ylab("Diffusion component 2") +  
  theme_classic()

saveRDS(dm, paste0(output_directory,"/dm.rds"))
while (!is.null(dev.list()))  dev.off()
cat('Finished, you can see the results now!\n')

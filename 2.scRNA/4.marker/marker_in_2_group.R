#!/usr/bin/env Rscript
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))

# 定义解析命令行参数的函数
option_list <- list( 
		    make_option("--width", type = "character", help = "Plot width", metavar = "character"),
		    make_option("--height", type = "character", help = "Plot height", metavar = "character"),
		    make_option(c("-i", "--input"), type = "character", help = "Input Seurat object file (.rds or .Rdata)", metavar = "character"),
		      make_option(c("-g", "--group"), type = "character", help = "Specify the identity for comparison", metavar = "character"),
		        make_option(c("-t", "--cell_type_1"), type = "character", help = "First cell type for comparison", metavar = "character"),
		        make_option(c("-e", "--cell_type_2"), type = "character", help = "Second cell type for comparison", metavar = "character"),
			  make_option(c("-f", "--logFCfilter"), type = "numeric", default = 0.5, help = "Log fold change filter (default: 0.5)", metavar = "numeric"),
			  make_option(c("-p", "--adjPvalFilter"), type = "numeric", default = 0.05, help = "Adjusted p-value filter (default: 0.05)", metavar = "numeric"),
			    make_option(c("-m", "--min_pct"), type = "numeric", default = 0.25, help = "Min percentage of cells expressing the gene (default: 0.25)", metavar = "numeric"),
			    make_option(c("-o", "--output"), type = "character", default = "results", help = "Output directory", metavar = "character")
			    )

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 加载 Seurat 对象
seurat_obj_path <- opt$input
scRNA <- readRDS(seurat_obj_path)
Group = opt$group
width = as.numeric(opt$width)
height = as.numeric(opt$height)
output_dir = opt$output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
  message("Output directory created: ", output_dir)
}

# 确保群体标识符存在
if (! Group %in% colnames(scRNA@meta.data)) {
	stop(paste0(Group,'not in meta data column'))
  Idents(scRNA) <- Group
}

# 获取两种细胞类型的细胞
ident.1 <- WhichCells(scRNA, idents = as.character(opt$cell_type_1))
ident.2 <- WhichCells(scRNA, idents = as.character(opt$cell_type_2))

# 查找 marker genes
mfb.markers <- FindMarkers(object = scRNA, 
			                              min.pct = opt$min_pct,
						                                 logfc.threshold = opt$logFCfilter,
						                                 ident.1 = ident.1, 
										                            ident.2 = ident.2)

# 筛选显著的 marker genes
sig.markers <- mfb.markers[(abs(as.numeric(as.vector(mfb.markers$avg_log2FC))) > opt$logFCfilter & 
			                                as.numeric(as.vector(mfb.markers$p_val_adj)) < opt$adjPvalFilter),]

# 提取前10个上调和下调基因
top_up <- sig.markers %>% 
  filter(avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(10)

top_down <- sig.markers %>% 
  filter(avg_log2FC < 0) %>% 
  arrange(avg_log2FC) %>% 
  head(10)

top_genes <- rbind(top_up, top_down)

# 绘制热图
if(nrow(top_genes) > 0) {
	#cells_to_plot <- c(ident.1, ident.2)
	#subset_data <- subset(scRNA, cells = cells_to_plot)
	#subset_data <- subset_data[, Idents(subset_data) %in% c(opt$cell_type_1, opt$cell_type_2)]
	# 设置细胞身份为比较组，并筛选只保留两个组
	#Idents(subset_data) <- Group
	subset_data <- subset(scRNA, idents = c(opt$cell_type_1, opt$cell_type_2))  # 添加这行来筛选

	# 绘制热图
	heatmap_file <- paste0(output_dir, '/', opt$cell_type_1, '_vs_', opt$cell_type_2, '_heatmap.pdf')
	pdf(file = heatmap_file, width = width, height = height)
	print(DoHeatmap(object = subset_data,
	                features = rownames(top_genes),
	                group.by = "ident",  # 使用ident而不是Group
	                label = TRUE) +
	        scale_fill_gradientn(colors = c("#BEBEBE", "#F5F5F5", "#CD2626")))
	dev.off()
  
  cat("Heatmap has been saved to", heatmap_file, "\n")
}

# 保存显著的 marker genes
write.csv(sig.markers, file = paste0(output_dir,'/',opt$cell_type_1,'_vs_',opt$cell_type_2,'.csv'),row.names = TRUE, quote = FALSE)

cat("Significant markers have been saved to", opt$output, "\n")

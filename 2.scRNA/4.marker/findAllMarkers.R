suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

find_significant_markers <- function(scRNA_filt,cluster_group, logFCfilter = 0.5, adjPvalFilter = 0.05, min_pct = 0.25, output_dir = "results",plot_width=16,plot_height=14) {
	  # 确保输出目录存在
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
      message("Output directory created: ", output_dir)
    }

    # 计算 marker genes
    pbmc.markers <- FindAllMarkers(object = scRNA_filt,
                                   group.by = cluster_group,
				                           only.pos = FALSE,
								                   min.pct = min_pct,
								                   logfc.threshold = logFCfilter)

    sig.markers <- pbmc.markers[(as.numeric(as.vector(pbmc.markers$avg_log2FC)) > logFCfilter & 
				                                      as.numeric(as.vector(pbmc.markers$p_val_adj)) < adjPvalFilter),]

    # 选择每个 cluster 中 top10 基因
    top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

	  # 生成输出文件路径
	  heatmap_file <- file.path(output_dir, "Heatmap.pdf")
	  output_file_path <- file.path(output_dir, "markers.csv")

	  # 生成热图
	  pdf(file = heatmap_file, width = plot_width, height = plot_height)
    print(DoHeatmap(object = scRNA_filt, features = top10$gene, label = FALSE) +
              scale_fill_gradientn(colors = c("#BEBEBE", "#F5F5F5", "#CD2626")) +
	          theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12)))
	  while (!is.null(dev.list()))  dev.off()

    # 保存显著性 marker 基因
    write.table(sig.markers, file = output_file_path, sep = ",", row.names = FALSE, quote = FALSE)

    cat("The significant markers have been saved to", output_file_path, "\n")
    cat("Heatmap has been saved to", heatmap_file, "\n")

    return(sig.markers)
}

# 主程序入口
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Please provide the Seurat object file path as a command line argument.")
}

seurat_obj_path <- args[1]
scRNA_filt <- readRDS(seurat_obj_path)
cluster_group <- args[2]

if (! cluster_group %in% colnames(scRNA_filt@meta.data)){
    stop(paste0(cluster_group,'not in meta.data column'))
}

logFCfilter <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.5)
adjPvalFilter <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)
min_pct <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.25)
output_dir <- ifelse(length(args) >= 6, args[6], "results")  # 设置默认输出目录
plot_width <- ifelse(length(args) >= 7, as.numeric(args[7]), 16)
plot_height <- ifelse(length(args) >= 8, as.numeric(args[8]), 14)


# 调用函数
result <- find_significant_markers(scRNA_filt,cluster_group, logFCfilter, adjPvalFilter, min_pct, output_dir,plot_width,plot_height)
cat('Finished, you can see the results now.\n')

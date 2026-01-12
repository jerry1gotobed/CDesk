suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

find_significant_markers <- function(scRNA_filt,cluster_group, logFCfilter = 0.5, adjPvalFilter = 0.05, min_pct = 0.25, output_dir = "results",plot_width=16,plot_height=14) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE) 
      message("Output directory created: ", output_dir)
    }
    #Idents(scRNA_filt) = scRNA_filt@meta.data$tsne_clusters
    # Calculate marker genes
    pbmc.markers <- FindAllMarkers(object = scRNA_filt,
                                   group.by = cluster_group,
				                           only.pos = FALSE,
								                   min.pct = min_pct,
								                   logfc.threshold = logFCfilter)

    sig.markers <- pbmc.markers[(as.numeric(as.vector(pbmc.markers$avg_log2FC)) > logFCfilter & 
				                                      as.numeric(as.vector(pbmc.markers$p_val_adj)) < adjPvalFilter),]

    # Top10 genes in each cluster
    top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

	  heatmap_file <- file.path(output_dir, "Heatmap.pdf")
	  output_file_path <- file.path(output_dir, "markers.csv")

	  # Generate heatmap
	  pdf(file = heatmap_file, width = plot_width, height = plot_height)
    print(DoHeatmap(object = scRNA_filt, features = top10$gene, label = FALSE,group.by = cluster_group) +
              scale_fill_gradientn(colors = c("#BEBEBE", "#F5F5F5", "#CD2626")) +
	          theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12)))
	  while (!is.null(dev.list()))  dev.off()

    # Save significant marker genes
    write.table(sig.markers, file = output_file_path, sep = ",", row.names = FALSE, quote = FALSE)

    cat("The significant markers have been saved to", output_file_path, "\n")
    cat("Heatmap has been saved to", heatmap_file, "\n")

    return(sig.markers)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Please provide the Seurat object file path as a command line argument.")
}

seurat_obj_path <- args[1]
scRNA_filt <- readRDS(seurat_obj_path)
cluster_group <- args[2]

# Check the meta column
if (! cluster_group %in% colnames(scRNA_filt@meta.data)){
    stop(paste0(cluster_group,'not in meta.data column'))
}

logFCfilter <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.5)
adjPvalFilter <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)
min_pct <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.25)
output_dir <- ifelse(length(args) >= 6, args[6], "results") 
plot_width <- ifelse(length(args) >= 7, as.numeric(args[7]), 16)
plot_height <- ifelse(length(args) >= 8, as.numeric(args[8]), 14)


# Call the function
result <- find_significant_markers(scRNA_filt,cluster_group, logFCfilter, adjPvalFilter, min_pct, output_dir,plot_width,plot_height)
cat('Finished, you can see the results now.\n')

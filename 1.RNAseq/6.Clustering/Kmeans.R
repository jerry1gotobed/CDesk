args <- commandArgs(trailingOnly = TRUE) 
count_file <- args[1]
meta_file <- args[2] 
output_dir <- args[3] 
clusters <- as.numeric(args[4]) 
method <- args[5] # correlation 
figure_height <- as.numeric(args[6]) # 10
figure_width <- as.numeric(args[7]) # 8
gene_list <- args[8] # no

k = clusters
# Load required libraries
suppressMessages(library(dplyr))
library(stats) %>% suppressMessages()
library(amap) %>% suppressMessages()
library(data.table) %>% suppressMessages()
library(randomcoloR) %>% suppressMessages()

# Check for required arguments
if (is.null(count_file) || is.null(meta_file)) {
  stop("Both count_file and meta_file are required arguments.")
}

# Define gene_expression_clustering function
gene_expression_clustering <- function(count_file, meta_file, output_dir, k,
                                       sample_col = "sample", group_col = "group",
                                       figure_height = 10, figure_width = 8,gene_list = 'no',
                                       method = "correlation") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Set working directory to output directory
  original_dir <- getwd()
  setwd(output_dir)

  # Read input files
  meta_data <- fread(meta_file, data.table = FALSE)
  Group = factor(meta_data[[group_col]],levels = unique(meta_data[[group_col]]))
  count_data <- read.csv(count_file,check.names = FALSE)
  row.names(count_data) = count_data[,1]
  count_data = count_data[,-1]
  
  if (!sample_col %in% colnames(meta_data)) {
      stop(paste("Error: Column", sample_col, "is missing in meta_data."))
  }
  if (!group_col %in% colnames(meta_data)) {
      stop(paste("Error: Column", group_col, "is missing in meta_data."))
  }

  # Extract sample IDs and group information
  sample_ids <- meta_data[[sample_col]]
  sample_group <- meta_data[[group_col]]

  # Select columns from count data
  selected_columns <- count_data[, c(sample_ids)]
  if (gene_list != 'no'){
    gene_list <- readLines(gene_list)
    selected_columns <- selected_columns[gene_list[gene_list %in% row.names(selected_columns)], ]
  }

  # Log2 transformation of expression data
  logcount <- log2(selected_columns + 1)
  logcount = as.data.frame(t(apply(logcount, 1, function(x) tapply(x, Group, mean, na.rm = TRUE))))
  pdf(file = "Elbow_plot.pdf")
  wss <- (nrow(logcount)-1)*sum(apply(logcount,2,var))
  for (i in 2:70) wss[i] <- sum(kmeans(logcount,centers=i)$withinss)
  plot(1:70, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
  while (!is.null(dev.list()))  dev.off()
  
  # K-means clustering
  # Perform K-means clustering
  km <- Kmeans(logcount, k, method = method, iter.max = 100, nstart = 50)

  clusterCol <- distinctColorPalette(k)

  # Write cluster gene lists to files
  for (id in seq(k)) {
    modGenes <- names(which(km$cluster == id))  # Get genes in current cluster
    fileName <- paste("cluster_", id, "_gene.Data.txt", sep = "")
    write.table(cbind(modGenes), file = fileName, row.names = FALSE,
                col.names = FALSE, quote = FALSE, sep = "\t")
  }

  # Generate cluster line plots (Figure 1)
  selected_cluster <- seq(k)
  pdf(file = "Clusters_line_plot.pdf", height = figure_height,width = figure_width)
  all_path <- levels(Group)
  for (each_i in seq(length(selected_cluster))) {
    each <- selected_cluster[each_i]
    modGenes <- names(which(km$cluster == each))

    # Skip if no genes in this cluster
    if (length(modGenes) == 0) {
      warning(paste("Cluster", each, "contains no genes. Skipping this cluster in plot."))
      next
    }

    v1 <- apply(logcount[modGenes, ], 2, mean)
    n <- length(modGenes)
    sd <- apply(logcount[modGenes, ], 2, sd)
    alpha <- 0.05
    v2 <- v1 - sd / sqrt(n) * qt(1 - alpha / 2, n - 1)
    v3 <- v1 + sd / sqrt(n) * qt(1 - alpha / 2, n - 1)
    par(mar = c(20, 4, 2, 2))
    plot(v1, lwd = 3, type = "l", col = clusterCol[each_i], ylim = c(min(v2), max(v3)),
         xlab = NA, ylab = "Log2(expression+1)", xaxt = "n",
         main = paste("C", each, "(n=", length(modGenes), ")", sep = ""))
    box(lwd = 2)
    axis(side = 1, 1:length(all_path), all_path, las = 2)
    polygon(c(1, 1:length(all_path), length(all_path):2),
            c(v2[1], v3, v2[length(all_path):2]),
            col = adjustcolor("grey", alpha.f = 0.4), border = NA)
  }
  while (!is.null(dev.list()))  dev.off()

  # Generate heatmap plots (Figure 2)
  pdf(file = "Clusters_heatmap.pdf", height = figure_height,width = figure_width)
  plot_data <- logcount
  for (each_i in seq(length(selected_cluster))) {
    each <- selected_cluster[each_i]
    modGenes <- names(which(km$cluster == each))
    plotMatrix <- plot_data[names(which(km$cluster == each)), ]

    # Skip if no genes in this cluster
    if (nrow(plotMatrix) == 0) {
      warning(paste("Cluster", each, "contains no genes. Skipping this cluster in plot."))
      next
    }

    all_exp <- c(as.matrix(plotMatrix))
    zmax <- quantile(all_exp, 0.99, na.rm = TRUE)
    zmin <- quantile(all_exp, 0.01, na.rm = TRUE)
    ColorRamp <- colorRampPalette(c("lightblue", "red"), bias = 1)(10000)
    ColorLevels <- seq(to = zmax, from = zmin, length = 10000)
    plotMatrix[plotMatrix < zmin] <- zmin
    plotMatrix[plotMatrix > zmax] <- zmax
    par(oma = c(0.5, 0.5, 0.5, 0.5), mar = c(20, 2, 2, 2))
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2), ncol = 9, nrow = 1, byrow = TRUE))
    image(1:ncol(plotMatrix), 1:nrow(plotMatrix), t(plotMatrix),
          xaxt = "n", yaxt = "n", col = ColorRamp, xlab = "", ylab = "")
    title(main = paste("C", each, "(n=", length(modGenes), ")", sep = ""), cex.main = 1)
    box(lwd = 2)
    axis(side = 1, 1:ncol(plotMatrix), labels = all_path, cex.axis = 1.2, las = 2)
    image(1, ColorLevels, t(matrix(data = ColorLevels, nrow = length(ColorLevels), ncol = 1)),
          col = t(ColorRamp), xlab = "", ylab = "", cex.axis = 2,
          xaxt = "n", yaxt = "n", useRaster = TRUE)
    #title(main = paste("C", each, "(n=", length(modGenes), ")", sep = ""), cex.main = 1)
    box(lwd = 2)
    axis(side = 2, c(zmin, round((zmax + zmin) / 2, 1), zmax),
         labels = c(round(zmin, 2), round((zmax + zmin) / 2, 1), round(zmax, 1)))
  }
  while (!is.null(dev.list()))  dev.off()

  # Reset working directory
  setwd(original_dir)

  # Return clustering results
  return(list(
    cluster = km$cluster,
    centers = km$centers,
    size = km$size,
    cluster_files = paste(output_dir, "/cluster_", seq(k), "_gene.Data.txt", sep = ""),
    figure1 = paste(output_dir, "/Clusters_line_plot.pdf", sep = ""),
    figure2 = paste(output_dir, "/Clusters_heatmap.pdf", sep = "")
  ))
}

valid_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary",
                   "pearson", "abspearson", "abscorrelation", "correlation",
                   "spearman", "kendall")
if (!any(sapply(valid_methods, function(x) grepl(x, method, ignore.case = TRUE)))) {
  stop("Invalid method: '", method, "'. Must be one of: ", paste(valid_methods, collapse=", "))
}

# Print execution information
cat("\n============= Gene Expression Clustering =============\n")
cat("Running with the following parameters:\n")
cat("Count file:", count_file, "\n")
cat("Metadata file:", meta_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Number of clusters (k):", k, "\n")
cat("Figure height:", figure_height, "\n")
cat("Figure width:", figure_width, "\n")
cat("Clustering method:", method, "\n")
cat("Gene list:", gene_list, "\n")
cat("=====================================================\n\n")

# Run the clustering function
results <- gene_expression_clustering(
  count_file = count_file,
  meta_file = meta_file,
  output_dir = output_dir,
  k = k,
  figure_height = figure_height,
  figure_width = figure_width,
  method = method,
  gene_list = gene_list
)

# Print results summary
cat("\nClustering completed successfully!\n")
cat("Generated", k, "clusters with the following sizes:\n")
print(results$size)
cat("\nOutput files:\n")
cat("- Cluster gene lists have been saved to:", output_dir, "as cluster_*_gene.Data.txt\n")
cat("- Line plots figure has been saved to:", results$figure1, "\n")
cat("- Heatmap figure has been saved to:", results$figure2, "\n")

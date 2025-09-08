suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))

rm(list = ls())

draw_heatmap = function(exp,gene,meta,width,height,output_dir,top_num){
  exp = exp[,meta$sample]
  gene = gene[gene %in% row.names(exp)]
  if (any(gene != '')){
    exp = exp[gene,]
    data = log2(exp+1)
    annotation_col = data.frame(group=meta$group)
    row.names(annotation_col) = meta$sample
    pheatmap::pheatmap(data, 
                       show_colnames = T,
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       annotation_col =annotation_col, 
                       annotation_legend=TRUE, 
                       show_rownames = T,
                       #scale = "row", 
                       color = colorRampPalette(c('#318fc4',"white" ,'#ca2b2b'))(200),
                       cellwidth = 10, 
                       cellheight = 10,
                       fontsize = 10, 
                       filename = paste0(output_dir,'/','heatmap.pdf'),
                       width = width,
                       height = height
    )
  }else{
    row_variances <- apply(exp, 1, var)
    top_genes <- exp[order(row_variances, decreasing = TRUE)[1:top_num], ]
    gene = row.names(top_genes)
    exp = exp[gene,]
    data = log2(exp+1)
    column_annotation <- HeatmapAnnotation(
      Group = meta$group)
    pdf(paste0(output_dir,'/','heatmap.pdf'),width = width,height = height)
    ht = Heatmap(as.matrix(data),
                 cluster_columns = TRUE,  
                 cluster_rows = TRUE,     
                 show_row_names = FALSE,  
                 show_column_names = TRUE, 
                 top_annotation = column_annotation
    )
    draw(ht)
    while (!is.null(dev.list()))  dev.off()
  }
}

sys_argv <- commandArgs(T)
method <- sys_argv[1]
if (method=='heatmap'){
  count <- sys_argv[2]
  gene <- sys_argv[3]
  meta <- sys_argv[4]
  top_num <- as.numeric(sys_argv[5])
  output_dir <- sys_argv[6]
  width <- as.numeric(sys_argv[7])
  height <- as.numeric(sys_argv[8])
  
  exp = read.csv(count,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
  row.names(exp) = exp[,1]
  exp = exp[,-1]
  exp = exp[rowSums(exp)>0,]  
  
  meta = read.csv(meta,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
  if (any(!(meta$sample %in% colnames(exp)))){
    stop(paste0('Sample not found:',meta$sample[!(meta$sample %in% colnames(exp))]))
  }
  
  if (any(! c('sample','group') %in% colnames(meta))){
    stop('sample group columns must be in the sample information file')
  }

  if (gene != 'NO'){
    gene = trimws(readLines(gene))
  }
  
  draw_heatmap(exp,gene,meta,width,height,output_dir,top_num)
}else if(method=='MA'){
  gfold <- sys_argv[2]
  gene <- sys_argv[3]
  gfold_thre <- as.numeric(sys_argv[4])
  name <- sys_argv[5]
  output_dir <- sys_argv[6]
  width <- as.numeric(sys_argv[7])
  height <- as.numeric(sys_argv[8])
  
  diff <- read.csv(gfold, sep=',', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Significant
  k1 <- diff$GFOLD > gfold_thre
  k2 <- diff$GFOLD < -gfold_thre
  diff$sig <- ifelse(k1, "Up", ifelse(k2, "Down", "NotSig"))
  
  # Count
  up_count <- sum(diff$sig == "Up")
  down_count <- sum(diff$sig == "Down")
  notsig_count <- sum(diff$sig == "NotSig")
  
  # Custom labels
  custom_labels <- c(
    paste0("Up (n=", up_count, ")"),
    "NotSig",
    paste0("Down (n=", down_count, ")")
  )
  
  # Calculate log2 mean RPKM
  diff$log2mean <- (log2(diff$RPKM1+1) + log2(diff$RPKM2+1)) / 2
  diff$sig <- factor(diff$sig, levels = c("Down", "NotSig", "Up"))
  
  if (gene != 'NO'){
    gene = trimws(readLines(gene))
  }
  missing = setdiff(gene,diff$Gene)
  cat('Gene not in result:',missing,'\n')
  mark = gene[gene%in%diff$Gene]
  diff$label <- ifelse(diff$Gene %in% mark, diff$Gene, NA)
  
  # Plot
  p <- ggplot(diff, aes(x = log2mean, y = GFOLD)) +
    geom_point(aes(color = sig)) +
    scale_color_manual(
      values = c("Up" = "red", "NotSig" = "gray", "Down" = "green"),
      labels = custom_labels,
      limits = c("Up", "NotSig", "Down")
    ) +
    labs(
      x = "Log2 Mean RPKM Expression",
      y = "GFOLD",
      title = paste0("MA plot of ",name)
    ) +
    geom_point(data = filter(diff, diff$Gene %in% mark),
               aes(color = sig),
               #size = 2, 
               shape = 21,  
               stroke = 0.2,  
               color = gray(0.5)) +
    geom_text_repel(data = filter(diff, !is.na(label)),
                    aes(x = log2mean,
                        y = GFOLD,
                        label = label),                       
                    size=3,                            
                    box.padding=unit(0.5,'lines'),          
                    point.padding=unit(0.1, 'lines'), 
                    segment.color='black',                 
                    show.legend=FALSE,
                    max.overlaps = Inf,na.rm = TRUE) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = 1, color = "blue", size = 0.8) +
    geom_hline(yintercept = c(-gfold_thre, gfold_thre), linetype = 2, color = "black", size = 0.6)+ theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.background = element_blank(),   
      panel.border = element_rect(fill = NA, color = "black")  
    )
  ggsave(filename = paste0(output_dir,'/',name,'_MA.pdf'),p,width = width, height = height)
}

library(ggplot2)
library(ggrepel)
library(dplyr)

sys_argv <- commandArgs(T)  
input_file <- sys_argv[1] 
pval_threshold <- as.numeric(sys_argv[2])
fc_threshold <- as.numeric(sys_argv[3])
gene_list <- unlist(strsplit(sys_argv[4], split = ","))
name <- sys_argv[5]
output_dir <- sys_argv[6]
y_type <- sys_argv[7]

data = read.csv(input_file,head=TRUE)
missing_genes <- setdiff(gene_list, data$gene_name)
# 如果有缺失的基因，输出到屏幕
if (length(missing_genes) > 0) {
  cat("The following genes are not found:\n")
  print(missing_genes)
}

down_count <- sum(data$sig == "down")
up_count <- sum(data$sig == "up")

# 自定义图例标签
custom_labels <- c(
  paste0("Down (n=", down_count, ")"),
  paste0("Normal"),
  paste0("Up (n=", up_count, ")")
)

pdf(paste0(output_dir,'/',name,'.pdf'))
if (y_type == 'padj'){ 
  ggplot(data,aes(log2FC, -log10(qvalue)))+
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "#999999")+
    geom_vline(xintercept = c(-fc_threshold,fc_threshold), linetype = "dashed", color = "#999999")+
    geom_point(aes(color = sig),
               size =0.6,
               alpha = 2) +
    theme_bw(base_size = 12)+
    #ggsci::scale_color_jama() +
    scale_color_manual(values = c("blue", "grey", "red"), labels = custom_labels)+
    theme(panel.grid = element_blank(),
          legend.position = 'right') +
    # 添加标签：
    geom_point(data = filter(data, gene_name %in% gene_list),
               aes(color = sig),
               size = 0.6,  # 调整点的大小
               shape = 21,  # 设置点为带边框的圆形
               stroke = 0.2,  # 边框宽度
               color = gray(0.5)) +  # 边框颜色
    geom_text_repel(data = filter(data, gene_name %in% gene_list),force_pull=100,
                    max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                    aes(label = gene_name, 
                        color = sig),
                    size = 2.5,show.legend = FALSE) +
    xlab("Log2FC")+
    ylab("-Log10(padj)")
}


if (y_type == 'GFOLD'){
  ggplot(data,aes(GFOLD, -log10(qvalue)))+
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "#999999")+
    geom_vline(xintercept = c(-fc_threshold,fc_threshold), linetype = "dashed", color = "#999999")+
    geom_point(aes(color = sig),
               size =0.6,
               alpha = 2) +
    theme_bw(base_size = 12)+
    #ggsci::scale_color_jama() +
    scale_color_manual(values = c("blue", "grey", "red"), labels = custom_labels)+
    theme(panel.grid = element_blank(),
          legend.position = 'right') +
    # 添加标签：
    geom_point(data = filter(data, gene_name %in% gene_list),
               aes(color = sig),
               size = 0.6,  # 调整点的大小
               shape = 21,  # 设置点为带边框的圆形
               stroke = 0.2,  # 边框宽度
               color = gray(0.5)) +  # 边框颜色
    geom_text_repel(data = filter(data, gene_name %in% gene_list),force_pull=100,
                    max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                    aes(label = gene_name, 
                        color = sig),
                    size = 2.5,show.legend = FALSE) +
    xlab("GFOLD")+
    ylab("-Log10(qvalue)")
}

dev.off()

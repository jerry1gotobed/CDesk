#! /usr/bin/Rscript


# Description : Perform the enrichment analysis based on the customization gene sets
# Author      : Zheting Zhang
# Time        : 2024-11-17 


args <- commandArgs(trailingOnly = TRUE)
gene_file <- args[1]
customer_term <- args[2]
output <- args[3]
high_color <- args[4]
low_color <- args[5]
R_Lib <- args[6]

.libPaths(c(R_Lib, .libPaths()))

# 检测并安装必要的包
suppressMessages(library(dplyr))
packages <- c( "ggplot2","clusterProfiler")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "BiocManager") {
      install.packages("BiocManager")
    } else {
      BiocManager::install(pkg)
    }
  }
}


if (!dir.exists(output)) dir.create(output, recursive = TRUE)
#Load packages
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()
# 定义物种数据库映射



  # 定义函数进行富集分析

  perform_enrichment <- function(gene_file,customer_term) {
  
  # 读取基因列表文件

  gene_list <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
  
  # 定义自定义term

  customer_term <- read.table(customer_term,header = FALSE,stringsAsFactors = FALSE,col.names = c("Name","Gene"),sep = "\t")
  
  if (ncol(customer_term)<2) {
    stop("customer file needs two colums: functional term and gene(symbol)")
  }
  
  customer_term_id <- unique(customer_term$Name)
  
  customer_term_id <- data.frame(Term=paste0("term_",c(1:length(customer_term_id))),Name=customer_term_id)
  
  customer_term <- merge(customer_term,customer_term_id,by="Name") %>% select(Term,Gene)
  

  # 进行富集分析

    enrich_result <- enricher(gene = gene_list$V1,
        		                  pAdjustMethod = "BH",
			                        pvalueCutoff = 0.05,
			                        qvalueCutoff = 0.5,
			                        TERM2GENE = customer_term,
			                        TERM2NAME = customer_term_id)
  
  

  
  # 输出结果表格

#  output_file <- paste0(gene_file,analysis_type, "_enrichment_results.txt")
#  write.table(as.data.frame(enrich_result), file = output_file, row.names = FALSE, sep = "\t")
  
  # 返回结果对象

  return(enrich_result)
}

file_name <- basename(gene_file)
file_prefix <- sub("\\.[^.]*$", "", file_name)



# 执行富集分析
enrichment_result <- perform_enrichment(gene_file, customer_term)

# 打印结果
output_file <- paste0(output,'/',file_prefix,"_custom_enrichment_results.txt")
write.table(as.data.frame(enrichment_result), file = output_file, row.names = FALSE, sep = "\t",quote=FALSE)

#################################################################################

plot_bubble <- function(enrich_result,high_color,low_color) {
  
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
    
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
      theme(text = element_text(size = 25))
    return(p)
  }
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      theme(text = element_text(size = 25)) 
    return(p)
  }}
# 定义函数绘制柱状图

plot_bar <- function(enrich_result) {
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description,fill=ONTOLOGY),stat = "identity")+
      theme_bw() +
      facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
      theme(text = element_text(size = 25))
    return(p)
  }
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description),stat = "identity")+
      theme_bw() +
      theme(text = element_text(size = 25))
    return(p)  
  }}

enrich_result = read.table(output_file,header = TRUE,sep = "\t")
enrich_result$Description <- factor(enrich_result$Description,levels = as.character(enrich_result$Description[rev(order(enrich_result$p.adjust))]))
enrich_result$GeneRatio <-  sapply(enrich_result$GeneRatio, cal <- function(x){out=eval(parse(text = x)); return(out)})

file_name = basename(output_file)
file_prefix <- sub("\\.[^.]*$", "", file_name)
# 绘制气泡图
pdf(paste0(output,'/',file_prefix,"_bubble.pdf"),width = 20,height = 10)
plot_bubble(enrich_result,high_color,low_color)
dev.off()

# 绘制柱状图
pdf(paste0(output,'/',file_prefix,"_bar.pdf"),width = 20,height = 10)
plot_bar(enrich_result)
dev.off()


cat("Thanks for using ! ^_^ , any question please contact with Bioinformatic Team of Pei\n")

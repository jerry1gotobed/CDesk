#! /usr/bin/Rscript

# Description : Perform the GO/KEGG enrichment analysis based on the pubsihed GO/KEGG gene sets
# Author      : Zheting Zhang
# Time        : 2024-11-17
args <- commandArgs(trailingOnly = TRUE)

gene_file <- args[1]
analysis_type <- args[2]
if (analysis_type == "GO") {
  species <- args[3]
  cutoff <- as.numeric(args[4])
  output = args[5]
  if (cutoff <= 0 || cutoff > 1){
          stop("cutoff value must be [0:1]")
  }
  high_color <- args[6]
  low_color <- args[7]
  R_Lib <- args[8]
}

if (analysis_type == "KEGG") {
  species <- args[3]
  output = args[4]
  high_color <- args[5]
  low_color <- args[6]
  R_Lib <- args[7]
}
 
.libPaths(c(R_Lib, .libPaths()))
# 检测并安装必要的包
packages <- c("BiocManager", "clusterProfiler", "org.Hs.eg.db", "org.Rn.eg.db", "org.Mm.eg.db", "org.Gg.eg.db", "org.Ss.eg.db", "enrichplot", "ggplot2","clusterProfiler")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "BiocManager") {
      install.packages("BiocManager")
    } else {
      BiocManager::install(pkg)
    }
  }
}

#Load packages
suppressMessages(library(dplyr))
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(org.Mm.eg.db) %>% suppressMessages() 
library(org.Gg.eg.db)
library(org.Ss.eg.db)
library(enrichplot) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()

if (!dir.exists(output)) dir.create(output, recursive = TRUE)

# 定义物种数据库映射
species_dbs <- list(
  "human" = org.Hs.eg.db,
  "rat" = org.Rn.eg.db,
  "mouse" = org.Mm.eg.db,
  "chicken" = org.Gg.eg.db,
  "pig" = org.Ss.eg.db)


# 定义函数进行GO或KEGG富集分析
organisms = c('hsa','mmu','rno','gga','ssc')
names(organisms) = c('human','mouse','rat','chicken','pig')

perform_enrichment <- function(gene_file, analysis_type, species) {
  # 读取基因列表文件
  gene_list <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
  # 选择物种数据库
  org_db <- species_dbs[[species]]

  # transform geneID
  gene.df <- bitr(geneID = gene_list$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)

  # 进行GO富集分析
  if (analysis_type == "GO") {
    enrich_result <- enrichGO(gene = gene.df$ENTREZID,
                              OrgDb = org_db,
                              ont = "ALL",
        		                  pAdjustMethod = "BH",
			                       pvalueCutoff = 0.05,
			                       qvalueCutoff = 0.5,
			                       keyType = "ENTREZID",
                              readable = TRUE)
  }
  
  # 进行KEGG富集分析
  if (analysis_type == "KEGG") {
    enrich_result <- enrichKEGG(gene = gene.df$ENTREZID,
                                organism = organisms[[species]],
			                         	keyType = "kegg",
			                        	pAdjustMethod = "BH",
			                         	minGSSize = 1, 
			                         	maxGSSize = 5000,
			                         	pvalueCutoff = 0.05,
		                       		qvalueCutoff = 0.5)
  }
  
  # 输出结果表格
#  output_file <- paste0(gene_file,analysis_type, "_enrichment_results.txt")
#  write.table(as.data.frame(enrich_result), file = output_file, row.names = FALSE, sep = "\t")
  
  # 返回结果对象
  return(enrich_result)
}

# 去冗余函数
simplify_enrichment <- function(enrich_result, cutoff) {
  if (nrow(enrich_result) > 0) {
    simplified_result <- simplify(enrich_result, cutoff = cutoff)
  } else {
    simplified_result <- enrich_result
    
  }
  return(simplified_result)
}

# 执行富集分析
enrichment_result <- perform_enrichment(gene_file, analysis_type, species)

# 合并冗余项
if (analysis_type == "GO") {
  enrichment_result <- simplify_enrichment(enrichment_result,cutoff)
}

# 层级筛选
#enrichment_result <- gofilter(enrichment_result,level = level_low:level_high)

file_name <- basename(output)
file_prefix <- sub("\\.[^.]*$", "", file_name)

# 打印结果
if (analysis_type == "GO") {
  output_file <- paste0(output,'/',analysis_type,'_',cutoff,"cutoff","_enrichment_results.txt")
}
if (analysis_type == "KEGG") {
  output_file <- paste0(output,'/',analysis_type,"_enrichment_results.txt")
}

write.table(as.data.frame(enrichment_result), file = output_file, row.names = FALSE, sep = ";",quote=FALSE)

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

enrich_result = read.csv(output_file,header = TRUE,sep = ";", stringsAsFactors = FALSE)
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

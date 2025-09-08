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
  plot_width <- as.numeric(args[6])
  plot_height <- as.numeric(args[7])
}

if (analysis_type == "KEGG") {
  species <- args[3]
  output = args[4]
  plot_width <- as.numeric(args[5])
  plot_height <- as.numeric(args[6])
}

high_color = 'red'
low_color = 'blue'

#Load packages
suppressMessages(library(dplyr))
library(org.Hs.eg.db) %>% suppressMessages()
library(org.Rn.eg.db) %>% suppressMessages()
library(org.Mm.eg.db) %>% suppressMessages() 
library(org.Gg.eg.db) %>% suppressMessages()
library(org.Ss.eg.db) %>% suppressMessages()
library(enrichplot) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()

if (!dir.exists(output)) dir.create(output, recursive = TRUE)

species_dbs <- list(
  "human" = org.Hs.eg.db,
  "rat" = org.Rn.eg.db,
  "mouse" = org.Mm.eg.db,
  "chicken" = org.Gg.eg.db,
  "pig" = org.Ss.eg.db)


organisms = c('hsa','mmu','rno','gga','ssc')
names(organisms) = c('human','mouse','rat','chicken','pig')

perform_enrichment <- function(gene_file, analysis_type, species) {
  # Read gene list
  gene_list <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
  # Species database
  org_db <- species_dbs[[species]]

  # transform geneID
  gene.df <- bitr(geneID = gene_list$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)

  # GO
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
  
  # KEGG
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
  
  return(enrich_result)
}

# Simplify
simplify_enrichment <- function(enrich_result, cutoff) {
  if (nrow(enrich_result) > 0) {
    simplified_result <- simplify(enrich_result, cutoff = cutoff)
  } else {
    simplified_result <- enrich_result
    
  }
  return(simplified_result)
}

# Perform enrichment
enrichment_result <- perform_enrichment(gene_file, analysis_type, species)

# Merge
if (analysis_type == "GO") {
  enrichment_result <- simplify_enrichment(enrichment_result,cutoff)
}

file_name <- basename(output)
file_prefix <- sub("\\.[^.]*$", "", file_name)

# Output
if (analysis_type == "GO") {
  output_file <- paste0(output,'/',analysis_type,'_',cutoff,"cutoff","_enrichment_results.csv")
}
if (analysis_type == "KEGG") {
  output_file <- paste0(output,'/',analysis_type,"_enrichment_results.csv")
}

enrichment_result <- as.data.frame(enrichment_result)  
enrichment_result <- data.frame(lapply(enrichment_result, function(col) {
  if (is.character(col)) {  
    return(gsub(",", ";", col))  # , -> ;
  } else {
    return(col)  
  }
}), stringsAsFactors = FALSE)  
write.table(as.data.frame(enrichment_result), file = output_file, row.names = FALSE, sep = ",",quote=FALSE)

#################################################################################

plot_bubble <- function(enrich_result,high_color,low_color) {
  
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
    
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
      theme(text = element_text(size = 20))
    return(p)
  }
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      theme(text = element_text(size = 20)) 
    return(p)
  }}
# 定义函数绘制柱状图

plot_bar <- function(enrich_result) {
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description,fill=ONTOLOGY),stat = "identity")+
      theme_bw() +
      facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
      theme(text = element_text(size = 20))
    return(p)
  }
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description,fill='orange'),stat = "identity",show.legend = FALSE)+
      theme_bw() +
      theme(text = element_text(size = 20))
    return(p)  
  }}

enrich_result = read.csv(output_file,header = TRUE,sep = ",", stringsAsFactors = FALSE)
if (length(enrich_result$Description)==0) {
  stop("No pathways meeting the criteria were found.")
}

enrich_result <- enrich_result[order(enrich_result$p.adjust), ]
enrich_result <- enrich_result[!duplicated(enrich_result$Description), ]
enrich_result$Description <- factor(enrich_result$Description,levels = as.character(enrich_result$Description[rev(order(enrich_result$p.adjust))]))
enrich_result$GeneRatio <-  sapply(enrich_result$GeneRatio, cal <- function(x){out=eval(parse(text = x)); return(out)})

file_name = basename(output_file)
file_prefix <- sub("\\.[^.]*$", "", file_name)
# Bubble plot
pdf(paste0(output,'/',file_prefix,"_bubble.pdf"),width = plot_width,height = plot_height)
plot_bubble(enrich_result,high_color,low_color)
while (!is.null(dev.list()))  dev.off()

# Bar plot
pdf(paste0(output,'/',file_prefix,"_bar.pdf"),width = plot_width,height = plot_height)
plot_bar(enrich_result)
while (!is.null(dev.list()))  dev.off()
cat("Thanks for using ! ^_^ , any question please contact with Bioinformatic Team of Pei\n")

args <- commandArgs(trailingOnly = TRUE)
gene_file <- args[1]
customer_term <- args[2]
output <- args[3]
plot_width <- as.numeric(args[4])
plot_height <- as.numeric(args[5])

high_color='red'
low_color='blue'

suppressMessages(library(dplyr))

if (!dir.exists(output)) dir.create(output, recursive = TRUE)
#Load packages
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()

perform_enrichment <- function(gene_file,customer_term) {
  gene_list <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)  
 
  # Custom term
  customer_term <- read.table(customer_term,header = FALSE,stringsAsFactors = FALSE,col.names = c("Name","Gene"),sep = "\t")
  
  if (ncol(customer_term)<2) {
    stop("customer file needs two colums: functional term and gene(symbol)")
  }
  
  customer_term_id <- unique(customer_term$Name)
  customer_term_id <- data.frame(Term=paste0("term_",c(1:length(customer_term_id))),Name=customer_term_id)
  customer_term <- merge(customer_term,customer_term_id,by="Name") %>% select(Term,Gene)
  
  # Enrichment analysis
  enrich_result <- enricher(gene = gene_list$V1,
        		                  pAdjustMethod = "BH",
			                        pvalueCutoff = 0.05,
			                        qvalueCutoff = 0.5,
			                        TERM2GENE = customer_term,
			                        TERM2NAME = customer_term_id)
  return(enrich_result)
}

file_name <- basename(gene_file)
file_prefix <- sub("\\.[^.]*$", "", file_name)

# Perform enrichment
enrichment_result <- perform_enrichment(gene_file, customer_term)

output_file <- paste0(output,'/',file_prefix,"_custom_enrichment_results.csv")
enrichment_result <- as.data.frame(enrichment_result)  
enrichment_result <- data.frame(lapply(enrichment_result, function(col) {
  if (is.character(col)) { 
    return(gsub(",", ";", col)) # , -> ;
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

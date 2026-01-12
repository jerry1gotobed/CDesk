args <- commandArgs(trailingOnly = TRUE)

gene_file <- args[1]
custom <- args[2]
if (custom == "no") {
  species <- args[3]
  output = args[4]
  plot_width <- as.numeric(args[5])
  plot_height <- as.numeric(args[6])
} else {
  customer_term = custom
  species = args[3]
  output = args[4]
  plot_width <- as.numeric(args[5])
  plot_height <- as.numeric(args[6])
}

high_color = 'red'
low_color = 'blue'

#Load packages
suppressMessages(library(dplyr))

library(enrichplot) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()
library(KEGG.db) %>% suppressMessages()
if (!dir.exists(output)) dir.create(output, recursive = TRUE)

if (species != 'no'){
  library(org.Hs.eg.db) %>% suppressMessages()
  library(org.Mm.eg.db) %>% suppressMessages()
  library(org.Rn.eg.db) %>% suppressMessages()
  library(org.Gg.eg.db) %>% suppressMessages()
  library(org.Ss.eg.db) %>% suppressMessages()
  if (!(species %in% c('human','mouse','rat','chicken','pig'))) {
    stop(paste('Species allowed:', toString(c('human','mouse','rat','chicken','pig'))))
  }
  species_dbs <- list(
    "human" = org.Hs.eg.db,
    "rat" = org.Rn.eg.db,
    "mouse" = org.Mm.eg.db,
    "chicken" = org.Gg.eg.db,
    "pig" = org.Ss.eg.db)
  
  organisms = c('hsa','mmu','rno','gga','ssc')
  names(organisms) = c('human','mouse','rat','chicken','pig')
}

perform_enrichment <- function(gene_file, species) {
  # Read gene list
  gene_list <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
  # Species database
  org_db <- species_dbs[[species]]
  
  # transform geneID
  gene.df <- bitr(geneID = gene_list$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
  
  # GO
  enrich_GO <- enrichGO(gene = gene.df$ENTREZID,
                        OrgDb = org_db,
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.5,
                        keyType = "ENTREZID",
                        readable = TRUE)
  
  enrich_GO_df <- as.data.frame(enrich_GO)
  go_empty <- nrow(enrich_GO_df) == 0
  
  # KEGG
  enrich_KEGG <- enrichKEGG(gene = gene.df$ENTREZID,
                            organism = organisms[[species]],
                            keyType = "kegg",
                            pAdjustMethod = "BH",
                            minGSSize = 1,
                            maxGSSize = 5000,
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.5,
			    use_internal_data = T)
  
  enrich_KEGG_df <- setReadable(enrich_KEGG, OrgDb = org_db, keyType="ENTREZID") 
  enrich_KEGG_df = enrich_KEGG_df@result

  kegg_empty <- nrow(enrich_KEGG_df) == 0
  
  # Check
  if (go_empty && kegg_empty) {
    stop("Error: Both GO and KEGG enrichment results are empty.")
  }
  
  if (!go_empty) {
    enrich_GO_final <- enrich_GO_df
  }
  
  if (!kegg_empty) {
    enrich_KEGG_df$ONTOLOGY <- 'KEGG'
    enrich_KEGG_final <- enrich_KEGG_df[, colnames(enrich_GO_df)]
  }

  if (!go_empty && !kegg_empty) {
    enrich_result <- rbind(enrich_GO_final, enrich_KEGG_final)
  } else if (!go_empty) {
    enrich_result <- enrich_GO_final
  } else {
    enrich_result <- enrich_KEGG_final
  }
  
  enrich_result$ONTOLOGY = factor(enrich_result$ONTOLOGY,levels = c('BP','CC','MF','KEGG'))
  enrich_result = enrich_result[enrich_result[['p.adjust']]<0.05,]
  return(enrich_result)
}

perform_custom_enrichment <- function(gene_file,customer_term) {
  gene_list <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)  
  
  # Custom term
  customer_term <- read.table(customer_term,header = FALSE,stringsAsFactors = FALSE,col.names = c("Name","Gene"),sep = "\t",quote = "")
  
  if (ncol(customer_term)<2) {
    stop("customer file needs two colums: functional term and gene(symbol)")
  }
  
  customer_term_id <- unique(customer_term$Name)
  customer_term_id <- data.frame(Term=paste0("term_",c(1:length(customer_term_id))),Name=customer_term_id)
  customer_term <- merge(customer_term,customer_term_id,by="Name") %>% dplyr::select(Term,Gene)
  
  # Enrichment analysis
  enrich_result <- enricher(gene = gene_list$V1,
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.5,
                            TERM2GENE = customer_term,
                            TERM2NAME = customer_term_id)
  
  enrich_result <- as.data.frame(enrich_result)
  enrich_result$ONTOLOGY = 'Custom' 
  enrich_result = enrich_result[enrich_result[['p.adjust']]<0.05,]
  return(enrich_result)
}

# Perform enrichment
if (custom == "no"){
  enrichment_result <- perform_enrichment(gene_file, species)
} else {
  enrichment_result <- perform_custom_enrichment(gene_file, customer_term)
}

# Output
output_file <- file.path(output,'enrichment_results.csv')

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
    enrich_result$ONTOLOGY = factor(enrich_result$ONTOLOGY,levels = c('BP','CC','MF','KEGG','Custom'))
    enrich_result <- enrich_result %>%
      group_by(ONTOLOGY) %>%
      arrange(pvalue, .by_group = TRUE) %>%
      slice_head(n = 10) %>%
      ungroup()
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
      theme(text = element_text(size = 10))
    return(p)
  }
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    enrich_result <- enrich_result[order(enrich_result$pvalue),]
    enrich_result <- head(enrich_result, 20)
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      theme(text = element_text(size = 10)) 
    return(p)
  }}

plot_bar <- function(enrich_result) {
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
    enrich_result$ONTOLOGY = factor(enrich_result$ONTOLOGY,levels = c('BP','CC','MF','KEGG','Custom'))
    enrich_result <- enrich_result %>%
      group_by(ONTOLOGY) %>%
      arrange(pvalue, .by_group = TRUE) %>%
      slice_head(n = 10) %>%
      ungroup()
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description,fill=ONTOLOGY),stat = "identity")+
      theme_bw() +
      facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
      theme(text = element_text(size = 10))
    return(p)
  }
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    enrich_result <- enrich_result[order(enrich_result$pvalue),]
    enrich_result <- head(enrich_result, 20)
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description,fill='orange'),stat = "identity",show.legend = FALSE)+
      theme_bw() +
      theme(text = element_text(size = 10))
    return(p)  
  }}

enrich_result = enrichment_result
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
pdf(file.path(output,'enrichment_bubble.pdf'),width = plot_width,height = plot_height)
plot_bubble(enrich_result,high_color,low_color)
while (!is.null(dev.list()))  dev.off()

# Bar plot
pdf(file.path(output,'enrichment_bar.pdf'),width = plot_width,height = plot_height)
plot_bar(enrich_result)
while (!is.null(dev.list()))  dev.off()
cat("Thanks for using ! ^_^ , any question please contact with Bioinformatic Team of Pei\n")

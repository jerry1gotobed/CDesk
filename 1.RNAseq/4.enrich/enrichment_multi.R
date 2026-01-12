args <- commandArgs(trailingOnly = TRUE)

gene_files <- args[1]
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
suppressMessages(library(data.table))
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

gene_files = fread(gene_files)
required_cols <- c("file", "tag")
missing_cols <- setdiff(required_cols, names(gene_files))

if (length(missing_cols) > 0) {
  stop("Miss column: ", paste(missing_cols, collapse = ", "))
}

gene_files$file <- as.character(gene_files$file)
gene_files$tag <- as.character(gene_files$tag)

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
                            qvalueCutoff = 0.5,use_internal_data = T)

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
  if (nrow(enrich_result) == 0){
    stop("Error: Empty enrichment results")
  }
  enrich_result$ONTOLOGY = 'Custom'
  enrich_result = enrich_result[enrich_result[['p.adjust']]<0.05,] 
  return(enrich_result)
}


result = data.frame()
# Perform enrichment
for (i in 1:nrow(gene_files)){
  gene_file = gene_files$file[i]
  if (custom == "no"){
    enrichment_result <- perform_enrichment(gene_file, species)
  } else {
    enrichment_result <- perform_custom_enrichment(gene_file, customer_term)
  }
  enrichment_result <- data.frame(lapply(enrichment_result, function(col) {
    if (is.character(col)) {  
      return(gsub(",", ";", col))  # , -> ;
    } else {
      return(col)  
    }
  }), stringsAsFactors = FALSE)  
  enrichment_result$Sample = gene_files$tag[i]
  result = rbind(result,enrichment_result)
}

# Output
output_file <- file.path(output,"enrichment_results_combine.csv")
write.table(as.data.frame(result), file = output_file, row.names = FALSE, sep = ",",quote=FALSE)

result1 <- result %>%
  group_by(Description, ONTOLOGY) %>%
  summarise(
    count = n(),
    mean_pvalue = mean(pvalue, na.rm = TRUE)
  ) %>%
  group_by(ONTOLOGY) %>%  # 按ONTOLOGY分组
  arrange(desc(count), mean_pvalue, .by_group = TRUE) %>%
  slice_head(n = 10) %>%  # 取每个组的前10个
  ungroup()  # 取消分组

result$Sample = factor(result$Sample,levels=gene_files$tag)
result = result[result$Description %in% result1$Description,]
result$Description = factor(result$Description,levels=rev(result1$Description))

pdf(file.path(output,"enrich_combine.pdf"),width = plot_width,height = plot_height)
if ('ONTOLOGY' %in% colnames(result)){
  result$ONTOLOGY = factor(result$ONTOLOGY,levels = c('BP','CC','MF','KEGG','Custom'))
  ggplot(data = result, mapping = aes(x = Description, y = Sample)) +
    geom_point(aes(size = Count, color = p.adjust)) +  # 图中点保持实心
    scale_size(range = c(2, 6),
               guide = guide_legend(override.aes = list(shape = 21,  # 图例显示为空心圆
                                                        fill = "white",
                                                        color = "black",
                                                        stroke = 0.8),
                                    order = 2)) +
    scale_colour_gradientn(
      colours = c("#E90F44", "#63ADEE"),
      limits = c(0, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free", switch = "y") +
    coord_flip() +
    labs(x='Pathway') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.1),
      strip.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
} else {
  ggplot(data = result, mapping = aes(x = Description, y = Sample)) +
    geom_point(aes(size = Count, color = p.adjust)) +  # 图中点保持实心
    scale_size(range = c(2, 6),
               guide = guide_legend(override.aes = list(shape = 21,  # 图例显示为空心圆
                                                        fill = "white",
                                                        color = "black",
                                                        stroke = 0.8),
                                    order = 2)) +
    scale_colour_gradientn(
      colours = c("#E90F44", "#63ADEE"),
      limits = c(0, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    coord_flip() +
    labs(x='Pathway') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.1),
      strip.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
}
while (!is.null(dev.list()))  dev.off()

cat("Thanks for using ! ^_^ , any question please contact with Bioinformatic Team of Pei\n")

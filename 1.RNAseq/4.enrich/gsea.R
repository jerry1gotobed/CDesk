rm(list = ls())

suppressMessages(library(dplyr))
library(org.Hs.eg.db) %>% suppressMessages()
library(org.Rn.eg.db) %>% suppressMessages()
library(org.Mm.eg.db) %>% suppressMessages()
library(org.Gg.eg.db) %>% suppressMessages()
library(org.Ss.eg.db) %>% suppressMessages()
library(enrichplot) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()
library(msigdbr) %>% suppressMessages()
library(patchwork) %>% suppressMessages()

sys_argv <- commandArgs(T)
deg <- sys_argv[1]
path <- sys_argv[2]
species <- sys_argv[3]
cols <- sys_argv[4]
output_dir <- sys_argv[5]
width <- as.numeric(sys_argv[6])
height <- as.numeric(sys_argv[7])

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

deg = read.csv(deg)
gene_col = strsplit(cols,',')[[1]][1]
FC_col = strsplit(cols,',')[[1]][2]
deg[[gene_col]] = toupper(deg[[gene_col]])

Species <- list(
  human = 'org.Hs.eg.db',
  rat = 'org.Rn.eg.db',
  mouse =  'org.Mm.eg.db',
  chicken = 'org.Gg.eg.db',
  pig = 'org.Ss.eg.db'
)

gene_entrezid <- bitr(geneID = deg[[gene_col]], 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = Species[[species]]
)

gene_entrezid$logFC <- deg[[FC_col]][match(gene_entrezid$SYMBOL, deg[[gene_col]])]
genelist = gene_entrezid$logFC
names(genelist) = gene_entrezid$ENTREZID 

Species <- list(
  human = 'Homo sapiens',
  rat = 'Rattus norvegicus',
  mouse =  'Mus musculus',
  chicken = 'Gallus gallus',
  pig = 'Sus scrofa'
)

m_t2g <- msigdbr(species = Species[[species]], category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

genelist = sort(genelist,decreasing = TRUE)

gsea_res <- GSEA(genelist, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH"
)

if (path=='NO'){
  pdf(paste0(output_dir,'/gsea.pdf'),width = width,height=height)
  p <- gseaplot2(gsea_res, geneSetID = 1:5)
  p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
  print(p[[1]] / p[[2]] / p[[3]])
  while (!is.null(dev.list()))  dev.off()
}else{
  path = trimws(readLines(path))
  missing = setdiff(path,gsea_res@result$Description)
  cat('Path not in result:',missing,'\n')
  mark = intersect(path,gsea_res@result$Description)
  if (length(mark)==0){
    stop('0 corresponding path')
  }
  pdf(paste0(output_dir,'/gsea.pdf'),width = width,height=height)
  p <- gseaplot2(gsea_res, geneSetID = mark)
  p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
  print(p[[1]] / p[[2]] / p[[3]])
  while (!is.null(dev.list()))  dev.off()
}

write.csv(gsea_res@result,file=paste0(output_dir,'/gsea_res.csv'))
cat('Finished, you can see the results now \n')

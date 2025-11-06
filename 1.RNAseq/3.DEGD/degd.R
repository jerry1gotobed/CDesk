rm(list = ls())

suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(rstatix))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggrepel))

sys_argv <- commandArgs(T)
method <- sys_argv[1]
count <- sys_argv[2]
output_dir <- sys_argv[3]
meta <- sys_argv[4]
Plot <- sys_argv[5]
gene <- sys_argv[6]
fc_thre <- as.numeric(sys_argv[7])
pval_thre <- as.numeric(sys_argv[8])
top_num <- as.numeric(sys_argv[9])
width <- as.numeric(sys_argv[10])
height <- as.numeric(sys_argv[11])

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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

deseq2 =function(exp,Group,pvalue_t,logFC_t){
  colData <- data.frame(row.names =colnames(exp),
                        condition=Group)
  dds <- suppressMessages(DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition))
  dds <- suppressMessages(DESeq(dds,quiet = TRUE))
  res <- suppressMessages(results(dds, contrast = c("condition",rev(levels(Group)))))
  DEG <- as.data.frame(res)
  DEG <- arrange(DEG,padj)
  DEG = na.omit(DEG)
  k1 = (DEG$padj < pvalue_t)&(DEG$log2FoldChange < -logFC_t)
  k2 = (DEG$padj < pvalue_t)&(DEG$log2FoldChange > logFC_t)
  DEG$sig = ifelse(k1,"Down",ifelse(k2,"Up","NotSig"))
  return(DEG)
}

adjust_t = function(exp,Group,pvalue_t,logFC_t){
  Group_temp <- factor(Group, levels = c(-1, 1), labels = c("Control", "Treat"))
  dfData = exp
  dfData$gene_name = row.names(exp)
  dfData <- select(dfData, gene_name, everything())
  dfClass = data.frame(Sample=names(Group),Class=as.character(Group_temp))
  df = dfData %>%
    as_tibble() %>%
    pivot_longer(-1,names_to = "Sample",values_to = "value")%>%
    left_join(dfClass,by=c("Sample" = "Sample"))

  dfFC = df %>%
    group_by(gene_name,Class) %>%
    summarise(mean = mean(value,na.rm=T)) %>%
    pivot_wider(names_from = Class,values_from = mean) %>%
    summarise(FC = Treat/Control)
  dfFC = dfFC[dfFC$FC!=0,]
  dfFC$log2FoldChange = log2(dfFC$FC)

  dfP <- df %>%
    group_by(gene_name) %>%
    summarise(
      p.value = tryCatch({
        t.test(value ~ Class, var.equal = TRUE)$p.value
      }, error = function(e) {
        return(NA)
      })
    )
  dfP = na.omit(dfP)

  dfP_FDR = dfP %>%
    select(1,last_col()) %>%
    mutate(padj = p.adjust(.$p.value,method = "BH"))  # p.adjust FDR,BH algorithm

  deg = merge(dfFC,dfP_FDR,by='gene_name')

  k1 = (deg$padj < pvalue_t)&(deg$log2FoldChange < -logFC_t)
  k2 = (deg$padj < pvalue_t)&(deg$log2FoldChange > logFC_t)
  deg$sig = ifelse(k1,"Down",ifelse(k2,"Up","NotSig"))
  return(deg)
}

draw_volcano = function(deg,mark,pvalue_t,logFC_t,width,height,output_dir){
  deg <- deg %>%
    mutate(log2FoldChange = ifelse(is.infinite(log2FoldChange), max(log2FoldChange[is.finite(log2FoldChange)]) * 2, log2FoldChange))
  down_count <- sum(deg$sig == "Down")
  up_count <- sum(deg$sig == "Up")
  custom_labels <- c(
    paste0("Down (n=", down_count, ")"),
    paste0("NotSig"),
    paste0("Up (n=", up_count, ")")
  )
  deg$label <- ifelse(row.names(deg) %in% mark, row.names(deg), NA)
  deg$sig <- factor(deg$sig, levels = c("Down", "NotSig", "Up"))
  deg <- deg %>% arrange(sig != "NotSig")
  p = ggplot(deg,aes(log2FoldChange,-1*log10(padj))) + 
    geom_point(aes(color = sig)) +                          
    labs(x="log2(FC)",
         y="-log10(Padj)") +
    scale_color_manual(values = c('Down'="blue", 'NotSig'="gray",'Up'= "red"),labels = custom_labels,limits = c("Down", "NotSig", "Up")) + 
    geom_hline(yintercept=-log10(pval_thre),linetype=2)+        
    geom_vline(xintercept=c(-fc_thre,fc_thre),linetype=2)+ 
    geom_point(data = filter(deg, row.names(deg) %in% mark),
               aes(color = sig),
               #size = 2,  
               shape = 21,  
               stroke = 0.2,  
               color = gray(0.5)) +
    geom_text_repel(data = filter(deg, !is.na(label)),
                    aes(x = log2FoldChange,
                        y = -1*log10(padj),
                        label = label),
                    size=3,                                  
                    box.padding=unit(0.5,'lines'),           
                    point.padding=unit(0.1, 'lines'),
                    segment.color='black',                   
                    show.legend=FALSE,
                    max.overlaps = Inf,na.rm = TRUE) + theme_bw()+
                    theme(
                      panel.grid.major = element_blank(),   
                      panel.grid.minor = element_blank(),  
                      panel.background = element_blank(),  
                      panel.border = element_rect(fill = NA, color = "black") 
                    )
  ggsave(filename = paste0(output_dir,'/',i,'_volcano.pdf'),p,width = width, height = height)
}

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

Cols = setdiff(colnames(meta),c('sample','group'))

for (i in Cols){
  if (any(! meta[[i]] %in% c(-1, 0, 1))) {
    stop("Only -1,1,0 is allowded in the compare columns")
    cat('\n')
  }
}

for (i in Cols){
  samples = meta[[i]]
  names(samples) = meta$sample
  samples = samples[samples!=0]
  count = exp[,names(samples)]
  Group <- factor(samples, levels = c(-1, 1))
  if (method=='deseq2'){
    deg <- deseq2(count, Group,pval_thre,fc_thre)
    deg$gene_name = row.names(deg)
    deg = select(deg, gene_name, everything())
  }else if(method=='adjusted_t'){
    deg <- adjust_t(count, Group,pval_thre,fc_thre)
    row.names(deg) = deg$gene_name
  }
  write.csv(deg,file = paste0(output_dir,'/',i,'_',method,'.csv'))

  if (Plot == 'True'){
    missing = setdiff(gene,row.names(deg))
    if (length(missing) > 0 && gene != "NO"){
      cat('Gene not in result:',missing,'\n')
    }
    mark = intersect(gene,row.names(deg))
    deg$label <- ifelse(row.names(deg) %in% mark, row.names(deg), NA)
    draw_volcano(deg,mark,pvalue_t,logFC_t,width,height,output_dir)
  }
}

if (Plot == 'True'){
  draw_heatmap(exp,gene,meta,width,height,output_dir,top_num)
}

cat('\n')
cat('Finished, you can see the results now')
cat('\n')

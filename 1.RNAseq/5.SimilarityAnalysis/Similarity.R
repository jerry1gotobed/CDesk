rm(list = ls())

sys_argv <- commandArgs(T)
sample1 <- sys_argv[1]
sample2 <- sys_argv[2]
group <- sys_argv[3]
output_dir <- sys_argv[4]
gene_list <- sys_argv[5]  # All
batch <- sys_argv[6]
plot_width <- as.numeric(sys_argv[7])
plot_height <- as.numeric(sys_argv[8])

suppressMessages(library("gplots"))
suppressMessages(library(edgeR))
suppressMessages(library(beeswarm))
suppressMessages(library(data.table))
suppressMessages(library(sva))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

cat("sample1:",sample1,"\n")
cat("sample2:",sample2,"\n")
cat("Grouping file:",group,"\n")
cat("Output directory:",output_dir,"\n")
cat("Gene list:",gene_list,"\n")
cat("Batch method:",batch,"\n")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
} 

options(scipen = 200)

###################################################################################################
##########################################     FUNCTION     #######################################
###################################################################################################
GeneExpressionPlot <-function(gene_list,out_name,matrix,CV2,unique_samples,bec_data,meta,plot_width,plot_height){
  
  pdf(out_name,width=plot_width,height=plot_height)
  for (each in gene_list){
    par(mar=c(14,4,4,2))
    layout(matrix(c(rep(1,7),rep(2,3)),nrow=1,ncol=10))
    
    plot_list = list()
    plot_list <- lapply(unique_samples,function(colo) {
      sample_names = meta[meta$tag == colo,]$sample
      data_subset <- as.matrix(bec_data[each, sample_names, drop = FALSE])
    })
    
    unique_meta <- meta %>% distinct(tag, .keep_all = TRUE)
    unique_meta$tag <- factor(unique_meta$tag, levels = unique_samples)
    unique_meta <- unique_meta[order(unique_meta$tag), ]
    
    beeswarm(plot_list,xaxt="n",col = sort(as.numeric(factor(unique_meta$group))),pch = 16,method = "swarm",corral = "wrap",main=each,xlab="",ylab="log(fpkm+1)",ylim=c(0,max(unlist(plot_list))+min(unlist(plot_list))),las=2)
    box(lwd=2)
    axis(side=1,at=1:length(plot_list),labels=unique_samples,las=2)
    plot_mean <- matrix[each,]
    points(seq(length(plot_mean)),plot_mean,type="b",col="black",pch="*",lwd=2)
    
    plot_vector <- CV2[each,]
    plot(seq(length(plot_vector)),plot_vector,xaxt="n",col=sort(as.numeric(factor(unique_meta$group))),type="h",lwd=3,xlab="",ylab="CV2",main="",ylim=c(0,max(plot_vector)+min(plot_vector)))
    axis(side=1,at=1:length(plot_vector),labels=unique_samples,las=2)
    box(lwd=2)
  }
  while (!is.null(dev.list()))  dev.off()
}

AllGeneExpressionPlot <-function(gene_list,out_name,matrix,CV2,unique_samples,bec_data,meta,plot_width,plot_height){
  pdf(out_name,width=plot_width,height=plot_height)
  par(mar=c(14,4,4,2))
  layout(matrix(c(rep(1,7),rep(2,3)),nrow=1,ncol=10))
  
  plot_list <- lapply(unique_samples, function(colo) {
    sample_names = meta[meta$tag == colo,]$sample
    #sample_names <- colnames(bec_data)[grepl(paste0("^", colo,"_"), colnames(bec_data))]
    data_subset <- as.matrix(bec_data[gene_list, sample_names, drop = FALSE])
    apply(data_subset, 2, mean, na.rm = TRUE)
  })
  
  # 计算每个unique_sample开头的样本数量
  unique_meta <- meta %>% distinct(tag, .keep_all = TRUE)
  unique_meta$tag <- factor(unique_meta$tag, levels = unique_samples)
  unique_meta <- unique_meta[order(unique_meta$tag), ]
  
  beeswarm(plot_list,xaxt="n",col = sort(as.numeric(factor(unique_meta$group))),pch = 16,method = "swarm",corral = "wrap",main="gene_list",xlab="",ylab="log(fpkm+1)",ylim=c(0,12),las=2)
  box(lwd=2)
  axis(side=1,at=1:length(plot_list),labels=unique_samples,las=2)
  plot_mean1 <- matrix[gene_list,]
  plot_mean <- apply(plot_mean1, 2, mean, na.rm = TRUE)
  points(seq(length(plot_mean)),plot_mean,type="b",col="black",pch="*",lwd=2)
  plot_vector1 <- CV2[gene_list,]
  plot_vector <- apply(plot_vector1, 2, mean, na.rm = TRUE)
  plot(seq(length(plot_vector)),plot_vector,xaxt="n",col = sort(as.numeric(factor(unique_meta$group))),type="h",lwd=3,xlab="",ylab="CV2",main="",ylim=c(0,12))
  axis(side=1,at=1:length(plot_vector),labels=unique_samples,las=2)
  box(lwd=2)
  while (!is.null(dev.list()))  dev.off()
}
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
data1 = fread(sample1, data.table = FALSE)
rownames(data1) = data1[[1]];data1 = data1[,-1]
data2 = fread(sample2, data.table = FALSE)
rownames(data2) = data2[[1]];data2 = data2[,-1]
meta = fread(group, data.table = FALSE)
meta$group = as.character(meta$group)
meta$tag = as.character(meta$tag)
meta <- meta[order(meta$group), ]

###################################################################################################
##############################             Remove batch            ################################
###################################################################################################
# Merge data
common_genes <- intersect(rownames(data1),rownames(data2))
data <- cbind(data1[common_genes,],data2[common_genes,])
all_samples <- colnames(data)
all_genes <- rownames(data)

# Sort by meta sample
data = data[,colnames(data) %in% meta$sample]
data = data[,meta$sample]

logfpkm = log2(data+1)
data = logfpkm

if (any(!meta$sample %in% colnames(data))){
  error_message <- paste0(
    "Can not find the sample in the sample information file: ",
    paste(meta$sample[!meta$sample %in% colnames(data)], collapse = ", "), "\n"
  )
  stop(error_message)
}

if (batch!='no'){
  check_col_names = colnames(meta)
  if (any(!sort(c("sample","group","tag","batch" )) == sort(check_col_names))){
    stop('Should have 4 columns in the sample information file：sample,group,tag,batch')
  }
}else {
  check_col_names = colnames(meta)
  if (any(!sort(c("sample","group","tag")) == sort(check_col_names[1:3]))){
    stop('Should have 3 columns in the sample information file：sample,group,tag')
  }
}

# Check batch
if (batch!='no'){
  if (!"batch" %in% colnames(meta)) {
    stop("Error: 'batch' column is missing in the meta data!")
  }
  # NA
  if (any(is.na(meta$batch))) {
    stop("Error: 'batch' column contains missing values!")
  }
  # Empty
  if (length(meta$batch) == 0 || all(meta$batch == "")) {
    stop("Error: 'batch' column is empty!")
  }
  if (length(unique(meta$batch)) == 1) {
    stop("Error: 'More than 1 batch")
  }
  if (batch=='removeBatchEffect'){
    experiment_batch <- as.factor(meta$batch)
    data <- removeBatchEffect(data,batch=experiment_batch)
  }else if (batch=='combat'){
    experiment_batch <- as.factor(meta$batch)
    data = ComBat(dat = data, batch = experiment_batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  }
}

bec_data = data

###################################################################################################
#######################     Get the average value of each sample type       #######################
###################################################################################################
mean_matrix <- c()
for (each in unique(meta$tag)) {
  select_group = meta[meta$tag == each,]$sample
  mean_matrix <- cbind(mean_matrix, rowMeans(bec_data[, select_group, drop = FALSE], na.rm = TRUE))
}
# data and dataref，variance，all gene matrix
var_matrix <- c()
for (each in unique(meta$tag)) {
  select_group = meta[meta$tag == each,]$sample
  var_matrix <- cbind(var_matrix, rowMeans(bec_data[, select_group, drop = FALSE], na.rm = TRUE))
}

CV2 <- var_matrix/(mean_matrix)^2
CV2[is.nan(CV2)] <- 0
unique_samples = unique(meta$tag)
colnames(CV2) <- unique_samples
colnames(var_matrix) <- unique_samples
colnames(mean_matrix) <- unique_samples

unique_meta <- meta %>% distinct(tag, .keep_all = TRUE)
unique_meta$tag <- factor(unique_meta$tag, levels = unique_samples)
unique_meta <- unique_meta[order(unique_meta$tag), ]

###################################################################################################
###############################                PCA               ##################################
###################################################################################################
# PCA
if (gene_list == "All") {
  plot_data <- mean_matrix
} else {
  genes = trimws(readLines(gene_list))
  plot_data <- mean_matrix[rownames(mean_matrix) %in% genes, ]
  genes = genes[genes %in% rownames(mean_matrix)]
}

plot_data = na.omit(plot_data)
# plot
pdf(paste0(output_dir,"/pca.pdf"),width=plot_width,height=plot_height)
pca <- prcomp(t(plot_data))

pca_data <- data.frame(
  PC1 = pca$x[, "PC1"],
  PC2 = pca$x[, "PC2"],
  Samples = unique_samples,
  group = unique_meta$group
)
# Plot
ggplot(pca_data, aes(x = PC1, y = PC2, label = Samples, color = group)) +
  geom_point(size = 3) +  
  geom_text_repel(size = 3, max.overlaps = Inf,segment.linetype = "dashed") + 
  labs(title = "PCA", x = "PC1", y = "PC2") +
  theme_classic()
while (!is.null(dev.list()))  dev.off()

###################################################################################################
####################                Correlation               #####################################
###################################################################################################
# Correlation
sample_cor <- cor(plot_data)
pdf(paste0(output_dir,"/cor.pdf"),width=plot_width,height=plot_height)
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",margins=c(13,13), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
cor_cell <- matrix(as.character(round(sample_cor, 2)), ncol=dim(sample_cor)[2])
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",cellnote=cor_cell, notecol="black", notecex=0.5, margins=c(13,13), cexRow=1, cexCol=1, revC=T, symm=T, distfun=
            function(c) as.dist(1 - c))
#dev.off()
while (!is.null(dev.list()))  dev.off()

pdf(paste0(output_dir,"/cor.legend.pdf"),width=3,height=1.8)
par(mar=c(4,2,4,2))
plotMatrix <- sample_cor
all_exp <- as.matrix(plotMatrix) # using same scale bar
zmax <- max(na.omit(all_exp))
zmin <- min(na.omit(all_exp))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="Correlation cofficient",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);
box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,2),round(zmax,2)))
while (!is.null(dev.list()))  dev.off()

###################################################################################################
###########################                line chart               ###############################
###################################################################################################
cat("-------------------------------------------------------------","\n")
sample = colnames(data)
if (gene_list=="All"){
  cat("No gene list provided, np bar figure plot.","\n")
} else{
  if (length(genes)<=20){
    GeneExpressionPlot(genes, paste0(output_dir,"/Afterbec.genelist.logfpkm.bar.pdf"),mean_matrix,CV2,unique_samples,bec_data,meta,plot_width,plot_height)
    cat('Plot the expression of each gene in the gene list \n')
  } else{
    AllGeneExpressionPlot(genes, paste0(output_dir,"/Afterbec.logfpkm.bar.pdf"),mean_matrix,CV2,unique_samples,bec_data,meta,plot_width,plot_height)
    cat("More than 20 genes in the gene list, plot the overall expression of the genes")
  }
}

cat("Done, you can see the results now \n")

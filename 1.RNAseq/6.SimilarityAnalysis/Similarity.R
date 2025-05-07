rm(list = ls())
###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
sys_argv <- commandArgs(T)  
input1 <- sys_argv[1] 
input2 <- sys_argv[2] 
Group <- sys_argv[3] # NO
output_dir <- sys_argv[4] 
gene_list <- sys_argv[5]  # ALL
R_Lib <- sys_argv[6]

.libPaths(c(R_Lib, .libPaths()))

library("gplots")
library(edgeR)
library(beeswarm)
library(readxl)

cat("Input base:", input1, "\n")
cat("Group information:", Group, "\n")
cat("Input ref:", input2, "\n")
cat("Output directory:", output_dir, "\n")
cat("Gene list:", gene_list, "\n")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建文件夹：", output_dir, "\n")
} 

cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)



###################################################################################################
##########################################     FUNCTION     #######################################
###################################################################################################
GeneExpressionPlot <-function(gene_list,out_name,matrix,CV2,unique_samples,bec_data,sample){
  pdf(paste(out_name,".pdf",sep=""),width=10,height=6)
  for (each in gene_list){
    par(mar=c(6,4,4,2))
    layout(matrix(c(rep(1,7),rep(2,3)),nrow=1,ncol=10))
    
    plot_list <- lapply(unique_samples, function(colo) {
      sample_names <- colnames(bec_data)[grepl(paste0("^", colo,"_"), colnames(bec_data))]
      data_subset <- as.matrix(bec_data[each, sample_names, drop = FALSE])
      apply(data_subset, 2, mean, na.rm = TRUE)
    })
    
    # 计算每个unique_sample开头的样本数量
    result <- c()
    for (us in unique_samples) {
      count <- sum(grepl(paste0("^", us, "_"), sample))
      result <- c(result,count)
    }
    
    ccco <- c()
    for (i in 1:length(result)) {
      ccco <- c(ccco, rep(i, result[i]))
    }
    
    beeswarm(plot_list,xaxt="n",col = ccco,pch = 16,method = "swarm",corral = "wrap",main=each,xlab="",ylab="log(fpkm+1)",ylim=c(0,12),las=2)
    box(lwd=2)
    axis(side=1,at=1:length(plot_list),labels=unique_samples,las=2)
    plot_mean <- matrix[each,]
    points(seq(length(plot_mean)),plot_mean,type="b",col="black",pch="*",lwd=2)
    plot_vector <- CV2[each,]
    plot(seq(length(plot_vector)),plot_vector,xaxt="n",col=ccco,type="h",lwd=3,xlab="",ylab="CV2",main="",ylim=c(0,12))
    axis(side=1,at=1:length(plot_vector),labels=unique_samples,las=2)
    box(lwd=2)
  }
  dev.off()
}

AllGeneExpressionPlot <-function(gene_list,out_name,matrix,CV2,unique_samples,bec_data,sample){
  pdf(paste(out_name,".pdf",sep=""),width=10,height=3.4)
  par(mar=c(6,4,4,2))
  layout(matrix(c(rep(1,7),rep(2,3)),nrow=1,ncol=10))
  
  plot_list <- lapply(unique_samples, function(colo) {
    sample_names <- colnames(bec_data)[grepl(paste0("^", colo,"_"), colnames(bec_data))]
    data_subset <- as.matrix(bec_data[gene_list, sample_names, drop = FALSE])
    apply(data_subset, 2, mean, na.rm = TRUE)
  })
  
  # 计算每个unique_sample开头的样本数量
  result <- c()
  for (us in unique_samples) {
    count <- sum(grepl(paste0("^", us, "_"), sample))
    result <- c(result,count)
  }
  
  ccco <- c()
  for (i in 1:length(result)) {
    ccco <- c(ccco, rep(i, result[i]))
  }
  
  beeswarm(plot_list,xaxt="n",col = ccco,pch = 16,method = "swarm",corral = "wrap",main="gene_list",xlab="",ylab="log(fpkm+1)",ylim=c(0,12),las=2)
  box(lwd=2)
  axis(side=1,at=1:length(plot_list),labels=unique_samples,las=2)
  plot_mean1 <- matrix[gene_list,]
  plot_mean <- apply(plot_mean1, 2, mean, na.rm = TRUE)
  points(seq(length(plot_mean)),plot_mean,type="b",col="black",pch="*",lwd=2)
  plot_vector1 <- CV2[gene_list,]
  plot_vector <- apply(plot_vector1, 2, mean, na.rm = TRUE)
  plot(seq(length(plot_vector)),plot_vector,xaxt="n",col=ccco,type="h",lwd=3,xlab="",ylab="CV2",main="",ylim=c(0,12))
  axis(side=1,at=1:length(plot_vector),labels=unique_samples,las=2)
  box(lwd=2)
  dev.off()
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
# 读取数据
data <- read.csv(input1, row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
dataref <- read.csv(input2, row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)

data_sample <- colnames(data)
dataref_sample <- colnames(dataref)
sample <- c(data_sample,dataref_sample)

if (Group != "NO"){
  group = read_excel(Group)
  group$tag = as.character(group$tag)
  data_unique_samples <- unique(group$tag[group$sample %in% data_sample])
  dataref_unique_samples <- unique(group$tag[group$sample %in% dataref_sample])
} else{
  data_unique_samples = data_sample
  dataref_unique_samples = dataref_sample
}


unique_samples <- c(data_unique_samples,dataref_unique_samples)



###################################################################################################
##############################                去批次               ################################
###################################################################################################
cat("-------------------------------------------------------------","\n")
# 对数据进行合并
common_genes <- intersect(rownames(data),rownames(dataref))
fpkm <- cbind(data[common_genes,],dataref[common_genes,])
all_samples <- colnames(fpkm)
all_genes <- rownames(fpkm)
logfpkm <- log2(fpkm+1)

# 构造样本类型
class1 <- all_samples[1:length(colnames(data))]
class2 <- all_samples[(length(colnames(data))+1):length(all_samples)]

# batch effect correction
experiment_batch <- as.factor(c(rep(1,length(class1)),rep(2,length(class2))))
bec_data <- removeBatchEffect(logfpkm,batch=experiment_batch)
cat("使用removeBatchEffect进行去批次","\n")



###################################################################################################
#######################                获取样本类型的平均数值               #######################
###################################################################################################
mean_matrix <- c()
for (prefix in unique_samples) {
  if (Group == "NO") {
    data_sample_names <- prefix
    dataref_sample_names <- prefix
  } else {
    data_sample_names <- group$sample[group$tag==prefix]
    dataref_sample_names <- group$sample[group$tag==prefix]
  }
  sample_names <- c(data_sample_names,dataref_sample_names)
  mean_matrix <- cbind(mean_matrix, rowMeans(bec_data[, sample_names, drop = FALSE], na.rm = TRUE))
}

# data和dataref，方差，全基因矩阵
var_matrix <- c()
for (prefix in unique_samples) {
  if (Group == "NO") {
    data_sample_names <- prefix
    dataref_sample_names <- prefix
  } else {
    data_sample_names <- group$sample[group$tag==prefix]
    dataref_sample_names <- group$sample[group$tag==prefix]
  }
  sample_names <- c(data_sample_names,dataref_sample_names)
  var_matrix <- cbind(var_matrix, apply(as.data.frame(bec_data[,sample_names]),1,var,na.rm=T))
}

CV2 <- var_matrix/(mean_matrix)^2
colnames(CV2) <- unique_samples
colnames(var_matrix) <- unique_samples
colnames(mean_matrix) <- unique_samples



###################################################################################################
###############################                PCA               ##################################
###################################################################################################
cat("-------------------------------------------------------------","\n")
# 根据不同样本组设置颜色
COLS <- rep("black",length(unique_samples))
names(COLS) <- unique_samples
COLS[data_unique_samples] <- cccol[1]
COLS[dataref_unique_samples] <- cccol[2]

# 根据不同样本组设置样式
PCHS <- rep(19,length(unique_samples))
names(PCHS) <- unique_samples

# PCA
if (gene_list == "ALL") {
  logfpkm_selected <- mean_matrix
} else {
  gene_table <- read.table(gene_list, sep = "\t")
  genes <- c(gene_table[[1]])
  logfpkm_selected <- mean_matrix[rownames(mean_matrix) %in% genes, ]
}

pdf(paste0(output_dir,"/Afterbec.genelist.logfpkm.pca.pdf"),width=10,height=10)
pca2 <- prcomp(t(logfpkm_selected))
plot(pca2$x[,"PC1"],pca2$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
text(pca2$x[,"PC1"],pca2$x[,"PC2"]+0.02,unique_samples,col=COLS,cex=0.3)
dev.off()
cat("PCA图绘制完成","\n")

###################################################################################################
####################                Correlation               #####################################
###################################################################################################
cat("-------------------------------------------------------------","\n")
# 对象数据的处理
data_logfpkm <- mean_matrix[,data_unique_samples]
dataref_logfpkm <- mean_matrix[,dataref_unique_samples]
if (gene_list == "ALL"){
    genes = row.names(dataref_logfpkm)
}else{
    genes = genes[genes %in% row.names(dataref_logfpkm)]
}

n_cor <- c()
if (gene_list == "ALL") {
  for (each_sample in data_unique_samples){
    tmp_cor <- c()
    for (each in rev(dataref_unique_samples)){
      tmp_cor <- c(tmp_cor,cor(dataref_logfpkm[,each],data_logfpkm[,each_sample]))
    }
    n_cor <- cbind(n_cor,tmp_cor)
  }
} else {
  for (each_sample in data_unique_samples){
    tmp_cor <- c()
    for (each in rev(dataref_unique_samples)){
      tmp_cor <- c(tmp_cor,cor(dataref_logfpkm[genes,each],data_logfpkm[genes,each_sample]))
    }
    n_cor <- cbind(n_cor,tmp_cor)
  }
}

rownames(n_cor) <- rev(dataref_unique_samples)
colnames(n_cor) <- data_unique_samples
plot_matrix <- n_cor

pdf(paste0(output_dir,"/Afterbec.genelist.logfpkm.correlation.pdf"),width=25,height=18)
all_exp <- c(as.matrix(plot_matrix))
zmax <- max(all_exp)
zmin <- min(all_exp)
par(mar=c(8,5,2,2))
layout(matrix(c(1,1,1,1,1,1,1,1,1,2),nrow=10,ncol=1,byrow=F))
ColorRamp <- colorRampPalette(c("white",cccol[2]), bias=1)(100)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=100)   #number sequence
image(1:ncol(plot_matrix), 1:nrow(plot_matrix), t(plot_matrix), axes=F, col=ColorRamp, xlab="",ylab="");box(lwd=2)
axis(side=1,1:ncol(plot_matrix),labels=data_unique_samples,las=2)
axis(side=2,1:nrow(plot_matrix),labels=rev(dataref_unique_samples),las=2)
text(matrix(rep(1:ncol(plot_matrix),nrow(plot_matrix)),ncol=ncol(plot_matrix),nrow=nrow(plot_matrix),byrow = T), matrix(rep(1:nrow(plot_matrix),ncol(plot_matrix)),ncol=ncol(plot_matrix),nrow=nrow(plot_matrix)),as.matrix(round(plot_matrix,2)))
par(mar=c(4,10,2,10))
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="log2(FPKM+1)",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,2),round(zmax,1)))
dev.off()
cat("相关性热力图绘制完成","\n")



###################################################################################################
###########################                line chart               ###############################
###################################################################################################
cat("-------------------------------------------------------------","\n")
if (gene_list=="ALL"){
  cat("因为没有gene_list,所以没有绘制bar figure...","\n")
} else{
  if (length(gene_list)<=20){
    GeneExpressionPlot(genes, paste0(output_dir,"/Afterbec.genelist.logfpkm.bar.pdf"),mean_matrix,CV2,unique_samples,bec_data,sample)
    cat("gene_list的每个基因的表达量绘制完成","\n")
  } else{
    AllGeneExpressionPlot(genes, paste0(output_dir,"/Afterbec.genelist.logfpkm.bar.pdf"),mean_matrix,CV2,unique_samples,bec_data,sample)
    cat("gene_list的长度大于等于20，所以绘制了所有gene_list中的整体表达量","\n")
  }
}

cat("分析完成，您可以查看结果了 \n")

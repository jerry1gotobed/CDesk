library("gplots")
library(readxl)
library(stringr)
###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
sys_argv <- commandArgs(T)  
input_file <- sys_argv[1] 
output_directory <- sys_argv[2]
group_file <- sys_argv[3] 

# 检查目录是否存在，如果不存在则创建该目录及其所有父目录
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)  # 这将创建所有必要的父目录
}

# 样本标识符最多14个，因为只设置了14个颜色
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
# 读取数据
data <- read.csv(input_file, row.names = 1)

# 对数据进行处理
all_samples <- colnames(data)
all_genes <- row.names(data)
logfpkm <- log2(data+1)

meta = read_excel(group_file)
meta$group = as.character(meta$group)

if (!identical(sort(meta$sample),sort(colnames(data)))){
  error_message <- paste0(
    "分组信息和样本不一致。\n",
    "以下分组信息元素不在样本中: ", 
    paste(meta$sample[!meta$sample %in% colnames(data)], collapse = ", "), "\n"
  )
  # 报错并退出程序
  stop(error_message)
}

# 构造样本类型
COLS = as.numeric(factor(meta$group))
COLS = cccol[COLS]
names(COLS) = all_samples

# 根据不同样本组设置样式,这里全部设置为19
PCHS <- rep(19,length(all_samples))
names(PCHS) <- all_samples

###################################################################################################
################################                PLOT               ################################
###################################################################################################

# PCA(Before版本)
pdf(paste0(output_directory,"/logfpkm.pca.pdf"),width=10,height=10)
pca <- prcomp(t(logfpkm))
plot(pca$x[,"PC1"],pca$x[,"PC2"],pch=PCHS,col=COLS,main="PCA",xlab="PC1",ylab="PC2")
text(pca$x[,"PC1"],pca$x[,"PC2"]+1,all_samples,col=COLS,cex=0.5)
dev.off()

# MDS(Before版本)
d <- dist(t(logfpkm)) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf(paste0(output_directory,"/logfpkm.MDS.pdf"),width=6,height=6)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS",col=COLS,pch=PCHS,lwd=4)
legend("center", unique(meta$group), col=unique(COLS),pch=19,bty="n",ncol = 2)
dev.off()

# correlation(Before版本)
sample_cor <- cor(data)
pdf(paste0(output_directory,"/fpkm.cor.pdf"),width=8,height=8)
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",margins=c(10,10), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
cor_cell <- matrix(as.character(round(sample_cor, 2)), ncol=dim(sample_cor)[2])
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",cellnote=cor_cell, notecol="black", notecex=0.5, margins=c(10,10), cexRow=1, cexCol=1, revC=T, symm=T, distfun=
            function(c) as.dist(1 - c))
dev.off()
pdf(paste0(output_directory,"/fpkm.cor.legend.pdf"),width=3,height=1.8)
par(mar=c(4,2,4,2))
plotMatrix <- sample_cor
all_exp <- as.matrix(plotMatrix) # using same scale bar
zmax <- max(na.omit(all_exp))
zmin <- min(na.omit(all_exp))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="Correlation cofficient",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);
box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,1),round((zmax+zmin)/2,1),round(zmax,1)))
dev.off()

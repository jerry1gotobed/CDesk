rm(list = ls())

sys_argv <- commandArgs(T)  # input_file,output_file,replicate_index
input_folder <- sys_argv[1]
output_dir <- sys_argv[1]
meta_file = sys_argv[2]
plot_width <- as.numeric(sys_argv[3])
plot_height <- as.numeric(sys_argv[4])

# 自动化脚本：从指定文件夹中读取所有以 `.SignalAcrossGenome.txt` 结尾的文件
library(tools)
suppressMessages(library(dplyr))

library("gplots") %>% suppressMessages()
library(stringr) %>% suppressMessages()
library(data.table) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
suppressMessages(library(edgeR))
suppressMessages(library(sva))

# 获取文件夹下所有以 `.SignalAcrossGenome.txt` 结尾的文件
file_list <- list.files(input_folder, pattern = "\\.SignalOnPeaks\\.txt$", full.names = TRUE)

# 读取文件并合并为数据框
data <- NULL
all_samples <- NULL

for (file in file_list) {
  # 读取文件的第一列qq
  sample_data <- read.table(file)[, 1]
  # 提取样本名称（去掉文件路径和后缀）
  sample_name <- file_path_sans_ext(basename(file))
  sample_name = gsub('\\.SignalOnPeaks','',sample_name)
  # 将数据合并到主数据框
  if (is.null(data)) {
    data <- sample_data
  } else {
    data <- cbind(data, sample_data)
  }
  # 保存样本名称
  all_samples <- c(all_samples, sample_name)
}

# 设置列名为样本名称
colnames(data) <- all_samples

meta = fread(meta_file, data.table = FALSE)

if (any(!meta$sample %in% colnames(data))){
  error_message <- paste0(
    "分组信息和样本不一致。\n",
    "以下分组信息元素不在样本中: ",
    paste(meta$sample[!meta$sample %in% colnames(data)], collapse = ", "), "\n"
  )
  # 报错并退出程序
  stop(error_message)
}

data = data[,colnames(data) %in% meta$sample]
data = data[,meta$sample]


check_col_names = colnames(meta)
if (any(!sort(c("sample","group","tag")) == sort(check_col_names[1:3]))){
  stop('分组信息表3列：sample,group,tag')
}

###################################################################################################
################################                PLOT               ################################
###################################################################################################
data = na.omit(data)
pca <- prcomp(t(data))

pdf(paste0(output_dir,"/CombinedPeaks.PCA.pdf"),width=plot_width,height=plot_height)
#par(mar=c(6,6,3,3))
#plot(pca$x[,"PC1"],pca$x[,"PC2"],col=COLS,pch=19,lwd=1,main="PCA",xlab=paste("PC1 (",round(summary(pca)$importance[2,1]*100,1),"%)",sep=""),ylab=paste("PC2 (",round(summary(pca)$importance[2,2]*100,1),"%)",sep=""))

#pdf(paste0(output_dir,"/CombinedPeals.PCA.legend.pdf"),width=5,height=5)
#plot(1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
#legend("center", unique(meta$rep), col=unique(COLS),pch=19,bty="n",ncol = 1)
#dev.off()

# 准备数据
pca_data <- data.frame(
  PC1 = pca$x[, "PC1"],
  PC2 = pca$x[, "PC2"],
  COLs = meta$group,
  Tag = meta$tag
)
# 使用 ggplot 绘制 PCA 图
ggplot(pca_data, aes(x = PC1, y = PC2, label = Tag, color = COLs)) +
  geom_point(size = 3) +  # 绘制散点
  geom_text_repel(size = 3, max.overlaps = Inf,segment.linetype = "dashed") +  # 自动调整标签位置，防止遮挡
  labs(title = "PCA", x = "PC1", y = "PC2") +
  theme_classic()
while (!is.null(dev.list()))  dev.off()


# correlation of all the samples
colnames(data) = str_replace(colnames(data),'.SignalOnPeaks','')
sample_cor <- cor(data,use="pair")
pdf(paste0(output_dir,"/CombinedPeaks.COR.pdf"),width = plot_width, height = plot_height)
cor_cell <- matrix(as.character(round(sample_cor, 2)), ncol=dim(sample_cor)[2])
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",cellnote=cor_cell, notecol="black", notecex=0.5, margins=c(8,8), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
while (!is.null(dev.list()))  dev.off()

pdf(paste0(output_dir,"/CombinedPeaks.legend.pdf"),width=3,height=1.8)
par(mar=c(4,2,4,2))
plotMatrix <- sample_cor
all_exp <- as.matrix(plotMatrix) # using same scale bar
zmax <- max(na.omit(all_exp))
zmin <- min(na.omit(all_exp))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="Pearson correlation cofficient",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,2),round(zmax,2)))
while (!is.null(dev.list()))  dev.off()

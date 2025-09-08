#! /usr/bin/env Rscript

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) != 2) {
  stop("Usage: Rscript script.R <input_file> <R_lib>")
}

# 读取输入文件
input_file <- args[1]
R_Lib <- args[2]
.libPaths(c(.libPaths(),R_Lib))

suppressMessages(library(ggplot2),classes = c("message","warning"))
suppressMessages(library(dplyr),classes = c("message","warning"))
suppressMessages(library(tidyr),classes = c("message","warning"))

data <- read.table(input_file, header = TRUE)

# 过滤数据
data <- data[rowSums(data[, 3:4]) > 0, ]

# 获取样本名称
sample1 <- colnames(data)[4]
sample2 <- colnames(data)[5]

# 重命名列
colnames(data)[c(4, 5)] <- c("s1", "s2")

# 计算 x 和 y 轴的限制
limitx <- quantile(data$s1, probs = c(0.0001, 0.9999))
limity <- quantile(data$s2, probs = c(0.0001, 0.9999))
axis_max <- max(limitx[2],limity[2])
axis_min <- min(limitx[1],limity[1])

# 计算 Spearman 和 Pearson 相关系数
spcor <- cor(data$s1, data$s2, method = "spearman") %>% round(digits = 2)
pscor <- cor(data$s1, data$s2, method = "pearson") %>% round(digits = 2)

# 设置输出文件路径
output_file <- gsub(".txt$", "", input_file)

# 绘制图形
pdf(paste0(output_file,".pdf"), width = 10, height = 8)
ggplot(data) +
  geom_point(aes(x = s1, y = s2, col = I("DeepSkyBlue"), alpha = I(0.02), size = I(10))) +
  geom_abline(linetype = "dashed", linewidth = 0.2, alpha = 0.5, intercept = 0, slope = 1) +
  theme_bw() +
  xlab(sample1) +
  ylab(sample2) +
  xlim(axis_min,axis_max) +
  ylim(axis_min,axis_max) +
  coord_fixed(ratio = 1) +
  annotate("text", x = axis_max * 0.25, y = axis_max * 0.9, size = 10, label = paste0("Pearson Cor=", pscor)) +
  annotate("text", x = axis_max * 0.25, y = axis_max*0.8,size=10,label=paste0("Spearman Cor=",spcor))+
  theme(text = element_text(size = 25))
  #+ggtitle(output_file)
while (!is.null(dev.list()))  dev.off()

# 打印完成信息

cat("Plot saved to", paste0(output_file,".pdf"), "\n")

#scaled_data

data$s1 <- log2(data$s1+0.1) %>% scale()
data$s2 <- log2(data$s2+0.1) %>% scale()

colnames(data)[c(4, 5)] <- c("s1", "s2")

limitx <- quantile(data$s1, probs = c(0.0001, 0.9999))
limity <- quantile(data$s2, probs = c(0.0001, 0.9999))
axis_max <- max(limitx[2],limity[2])
axis_min <- min(limitx[1],limity[1])

spcor <- cor(data$s1, data$s2, method = "spearman") %>% round(digits = 2)
pscor <- cor(data$s1, data$s2, method = "pearson") %>% round(digits = 2)

pdf(paste0(output_file,"scaled.pdf"),width = 10,height = 8)

ggplot(data)+geom_point(aes(x=s1,y=s2,col=I("DeepSkyBlue"),alpha=I(0.02),size=I(10)))+
  theme_bw()+
  xlab(sample1)+ylab(sample2)+
  xlim(limitx) +
  ylim(limity) +
  annotate("text",x=limitx[2]*0.2+limitx[1],y=limity[2]*0.9,size=10,label=paste0("Pearson Cor=",pscor))+
  annotate("text",x=limitx[2]*0.2+limitx[1],y=limity[2]*0.8,size=10,label=paste0("Spearman Cor=",spcor))+
  theme(text = element_text(size = 25))
  #+ggtitle(paste0(output_file,"_scaled"))

while (!is.null(dev.list()))  dev.off()

cat("Plot saved to", paste0(output_file,"scaled.pdf"), "\n")

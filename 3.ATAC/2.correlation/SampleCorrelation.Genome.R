rm(list = ls())

sys_argv <- commandArgs(T)  # input_file,output_file,replicate_index
input_folder <- sys_argv[1]
output_dir <- sys_argv[1]
meta_file = sys_argv[2]
plot_width <- as.numeric(sys_argv[3])
plot_height <- as.numeric(sys_argv[4])

# 自动化脚本：从指定文件夹中读取所有以 `.SignalAcrossGenome.txt` 结尾的文件
suppressMessages(library(dplyr))
library(tools)
library("gplots") %>% suppressMessages()
library(stringr) %>% suppressMessages()
library(data.table) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
suppressMessages(library(pheatmap))

input_folder = file.path(input_folder,'tmp')
# 获取文件夹下所有以 `.SignalAcrossGenome.txt` 结尾的文件
file_list <- list.files(input_folder, pattern = "\\.SignalAcrossGenome\\.txt$", full.names = TRUE)

# 读取文件并合并为数据框
data <- NULL
all_samples <- NULL

for (file in file_list) {
  # 读取文件的第一列qq
  sample_data <- read.table(file)[, 1]
  # 提取样本名称（去掉文件路径和后缀）
  sample_name <- file_path_sans_ext(basename(file))
  sample_name = gsub('\\.SignalAcrossGenome','',sample_name)
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
meta$sample = sub("\\..*$", "", basename(meta$bw))

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
colnames(data) = meta$tag

required <- c("sample", "group", "tag")
if (!any(required %in% colnames(meta))) {
  stop("Error: 数据框中必须至少包含以下一列: ",
       paste(required, collapse = ", "))
}

###################################################################################################
################################                PLOT               ################################
###################################################################################################
data = na.omit(data)
pca <- prcomp(t(data))

pdf(paste0(output_dir,"/GenomeBin.PCA.pdf"),width = plot_width, height = plot_height)
var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_pct <- round(var_explained[1], 1)
pc2_pct <- round(var_explained[2], 1)

pca_data <- data.frame(
  PC1     = pca$x[, "PC1"],
  PC2     = pca$x[, "PC2"]
)

# Plot
ggplot(pca_data, aes(x = PC1, y = PC2, label = meta$tag, color = meta$group)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf, segment.linetype = "dashed") +
  labs(
    title = "PCA",
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + scale_color_discrete(name = "Group")
while (!is.null(dev.list()))  dev.off()

# correlation of all the samples
sample_cor <- cor(data)
pdf(paste0(output_dir,"/GenomeBin.COR.pdf"),width = plot_width, height = plot_height)
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"))(100)

pheatmap(
  sample_cor,
  color            = ColorRamp,
  display_numbers  = TRUE,
  number_format    = "%.2f",
  fontsize_number  = 8,
  cluster_rows     = T,
  cluster_cols     = T,
  legend            = TRUE,
  legend_breaks    = c(0,0.5,1),
  legend_labels    = c("0","0.5","1"),
  border_color     = NA,
  show_rownames    = TRUE,
  show_colnames    = TRUE
)
while (!is.null(dev.list()))  dev.off()

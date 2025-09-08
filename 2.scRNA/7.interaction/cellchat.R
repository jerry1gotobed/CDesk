#!/usr/bin/env Rscript
suppressMessages(library(argparse))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(CellChat))
suppressMessages(library(NMF))
suppressMessages(library(ggalluvial))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(svglite))
suppressMessages(library(future))

# 创建命令行参数解析器
parser <- ArgumentParser()
parser$add_argument("--seurat_file", type="character", help="Path to the Seurat object file (RDS format)")
parser$add_argument("--output_directory", type="character", help="Output directory")
parser$add_argument("--species", type="character", help="Species (human or mouse)")
parser$add_argument("--threads", type="integer", default=16, help="Number of threads for parallel processing")
parser$add_argument("--group_by", type="character", help="Grouping variable (seurat_clusters or celltype......)")
parser$add_argument("--width", type="integer", default=8,help="Plot width, default: 8")
parser$add_argument("--height", type="integer", default=8,help="Plot height, default: 8")

# 解析命令行参数
args <- parser$parse_args()
options(future.globals.maxSize = 10 * 1024 * 1024 * 1024)  # 1 GiB

if (!dir.exists(args$output_directory)) {
  dir.create(args$output_directory, recursive = TRUE)  # 使用recursive = TRUE来确保创建多级目录
  message("Output directory created: ", args$output_directory)
}

# 读取Seurat对象
data_integrated <- readRDS(args$seurat_file)

if (!args$group_by %in% colnames(data_integrated@meta.data)){
  stop(paste0(args$group_by,' not in meta data column of ',args$seurat_file))
}

# 防止因为全是数字而报错
if (grepl('cluster',args$group_by)){
  data_integrated@meta.data$seurat_clusters <- factor(paste0("c", data_integrated@meta.data$seurat_clusters))
}
#data_integrated$celltype <- data_integrated@active.ident
data_integrated <- SetIdent(data_integrated, value = data_integrated@meta.data[[args$group_by]])


# 创建CellChat对象
cellchat <- createCellChat(data_integrated@assays$RNA@data, meta = data_integrated@meta.data, group.by = args$group_by)
groupSize <- as.numeric(table(cellchat@idents))

# 导入配体受体数据库
if (args$species == "human") {
  CellChatDB <- CellChatDB.human

} else if (args$species == "mouse") {
  CellChatDB <- CellChatDB.mouse

} else {
  stop("Unsupported species. Please specify 'human' or 'mouse'.")
}

unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# CellChat分析
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = args$threads)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- smoothData(cellchat, adj=PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
saveRDS(cellchat, file = paste0(args$output_directory, "/cellchat.rds"))

cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = paste0(args$output_directory, "/net_lr.csv"))

if (args$species == "human" || args$species == "mouse") {
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, file = paste0(args$output_directory, "/net_pathway.csv"))
}

# 可视化
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf(file = paste0(args$output_directory, "/cell_net_circle.pdf"), width=args$width, height=args$height)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Interaction weights/strength")
while (!is.null(dev.list()))  dev.off()

mat <- cellchat@net$count

pdf(file = paste0(args$output_directory, "/TIL_net_number_individual.pdf"), width=args$width, height=args$height)
par(mar=c(2,2,2,2))
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
while (!is.null(dev.list()))  dev.off()

mat <- cellchat@net$weight

pdf(file = paste0(args$output_directory, "/TIL_net_strength_individual.pdf"), width=args$width, height=args$height)
par(mar=c(2,2,2,2))
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
while (!is.null(dev.list()))  dev.off()

#gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = args$group_by, measure = "count")
#gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = args$group_by, measure = "weight")
#p <- gg1 + gg2

#ggsave(paste0(args$output_prefix, "_Overview_number_strength.pdf"), p, width = 8, height = 6)

cat("Pipeline completed successfully.\n")

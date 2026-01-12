suppressMessages(library(dplyr))
library("gplots") %>% suppressMessages()
library(stringr) %>% suppressMessages()
library(data.table) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
library(pheatmap) %>% suppressMessages()

sys_argv <- commandArgs(T)
input_dir = sys_argv[1]
meta = sys_argv[2]
meta = fread(meta)
output_directory = sys_argv[3]
width = as.numeric(sys_argv[4])
height = as.numeric(sys_argv[5])

tags <- meta$tag

dt_list <- lapply(tags, function(tag) {
  dt <- fread(file.path(input_dir, paste0(tag, ".oe.txt")))
  setnames(dt, names(dt)[7], tag)
  dt[, c(names(dt)[1:6], tag), with = FALSE]
})

res <- Reduce(function(x, y) {
  merge(x, y, by = names(dt_list[[1]])[1:6], all = FALSE)
}, dt_list)

data = res[,-1:-6]

###################################################################################################
################################                PLOT               ################################
###################################################################################################
options(scipen = 200)
data = na.omit(data)
meta = as.data.frame(meta)
row.names(meta) = meta$tag
lables <- meta[colnames(data),'tag']
colors <- meta[colnames(data),'group']
# PCA
pdf(paste0(output_directory,"/pca.pdf"),width=width,height=height)
pca <- prcomp(t(data))

var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_pct <- round(var_explained[1], 1)
pc2_pct <- round(var_explained[2], 1)

pca_data <- data.frame(
  PC1     = pca$x[, "PC1"],
  PC2     = pca$x[, "PC2"],
  Label  = lables,
  Color = colors
)

# Plot
ggplot(pca_data, aes(x = PC1, y = PC2, label = Label, color = Color)) +
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

# correlation
sample_cor <- cor(data)
pdf(paste0(output_directory,"/cor_heatmap.pdf"),width=width,height=height)
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
  legend            = T,      
  legend_breaks    = c(0,0.5,1),  
  legend_labels    = c("0","0.5","1"),
  border_color     = NA,
  show_rownames    = TRUE,
  show_colnames    = TRUE
)

while (!is.null(dev.list()))  dev.off()


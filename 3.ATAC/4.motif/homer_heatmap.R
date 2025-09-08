rm(list=ls())
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(circlize))
suppressMessages(library(stringr))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
input_file <- fread(args[1])
output_dir <- args[2]
homers = args[3]
width=as.numeric(args[4])
height=as.numeric(args[5])

homer_filter_final = data.frame()
for (file in input_file$dir){
  homer = fread(file.path(file,'knownResults.txt'))
  colnames(homer) = sub("\\(.+", "", colnames(homer))
  homer$Sample= which(input_file$dir == file)
  homer_filter_final = rbind(homer,homer_filter_final)
}

homer_filter_final$Motif_temp = homer_filter_final$`Motif Name`
homer_filter_final <- homer_filter_final %>%
  separate(
    col = Motif_temp,
    into = c("TF_Family", "Source", "Lib"),
    sep = "/",
    extra = "drop",
    fill = "right"
  )
# Fetch the most significant TF_family, remove duplicate
homer_no_dup <- homer_filter_final %>%
  group_by(Source, TF_Family) %>%
  mutate(avg_pvalue = mean(`Log P-value`, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(avg_pvalue)

homer_no_dup = homer_no_dup[!duplicated(homer_no_dup$TF_Family),]
homer_no_dup$TF = sub("\\(.+", "", homer_no_dup$TF_Family)
homer_no_dup = homer_no_dup[!duplicated(homer_no_dup$TF),]
if (homers != "no"){
  homers = homers <- readLines(homers)
  homer_no_dup = homer_no_dup[homer_no_dup$TF %in% homers,]
} else{
  homer_no_dup = head(homer_no_dup, 20)
}

homer_filter_final = homer_filter_final[homer_filter_final$`Motif Name` %in% homer_no_dup$`Motif Name`,]
homer_filter_final$`Log P-value` = -homer_filter_final$`Log P-value`

# Grab Sample,TF_Family and Log_P_value
df_logp <- homer_filter_final %>%
  select(Sample, TF_Family, `Log P-value`)

# Construct matrix
matrix_logp <- df_logp %>%
  pivot_wider(names_from = Sample, values_from = `Log P-value`)

matrix_logp = as.data.frame(matrix_logp)
matrix_logp = matrix_logp[order(matrix_logp$TF_Family),]
row.names(matrix_logp) = matrix_logp$TF_Family
matrix_logp <- matrix_logp[, colnames(matrix_logp) != "TF_Family"]

mycol <- colorRamp2(seq(0, 40, length.out = 4), 
                    brewer.pal(n = 7, name = "RdBu")[4:1])

###########
homer_filter_final <- homer_filter_final %>%
  rename(
    q_value = `q-value `  
  )

df_logp <- homer_filter_final %>%
  select(Sample, TF_Family, q_value)

# 构建矩阵
matrix_tag <- df_logp %>%
  pivot_wider(names_from = Sample, values_from = q_value)

#matrix_tag$TF_Family = sapply(str_split(matrix_tag$TF_Family, "\\("), `[`, 1)
matrix_tag = as.data.frame(matrix_tag)
matrix_tag = matrix_tag[order(matrix_tag$TF_Family),]
row.names(matrix_tag) = matrix_tag$TF_Family
matrix_tag = matrix_tag[, colnames(matrix_tag) != "TF_Family"]

###########
colnames(matrix_logp) = input_file$sample[as.numeric(colnames(matrix_logp))]
row_groups = sub(".*\\((.*)\\).*", "\\1", rownames(matrix_logp))
if (length(levels(factor(input_file$group)))==1){
  column_split = NULL
}else{
  column_split = input_file$group
}

pdf(file.path(output_dir,'homer_heatmap.pdf'),width=width,height=height)
row.names(matrix_logp) = sub("\\(.+", "", rownames(matrix_logp))
Heatmap(
  matrix = matrix_logp,
  col = mycol,
  border = TRUE, 
  rect_gp = gpar(col = "white", lwd = 1), 
  width = unit(10, "cm"),
  height = unit(10, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  split = row_groups,
  row_title_rot = 0,
  column_split = column_split,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (matrix_tag[i, j] < 0.001) {
      grid.text("+++", x, y, gp = gpar(fontsize = 6, col = 'gray20'))
    } else if (matrix_tag[i, j] < 0.01) {
      grid.text("++", x, y, gp = gpar(fontsize = 6, col = 'gray20'))
    } else if (matrix_tag[i, j] < 0.05) {
      grid.text("+", x, y, gp = gpar(fontsize = 6, col = 'gray20'))
    }
  },
  heatmap_legend_param = list(title = "-Log10(P.Value)")
)
while (!is.null(dev.list()))  dev.off()
cat('Finished, you can check the results now\n')

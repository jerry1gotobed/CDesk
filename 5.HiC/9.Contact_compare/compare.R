suppressMessages(library(HiCcompare))
suppressMessages(library(BiocParallel))
suppressMessages(library(Cairo))
suppressMessages(library(ggplot2))

sys_argv <- commandArgs(T)
output_dir = sys_argv[1]
prefix1 = sys_argv[2]
prefix2 = sys_argv[3]
chr = sys_argv[4]

####################################  HiC compare  ####################################
hic1 = read.csv(file.path(output_dir,'tmp',paste0(prefix1,'.txt')),sep = '\t',header = FALSE)
hic1 = hic1[c(3,6,7)]
colnames(hic1) = c('region1','region2','IF')
hic2 = read.csv(file.path(output_dir,'tmp',paste0(prefix2,'.txt')),sep = '\t',header = FALSE)
hic2 = hic2[c(3,6,7)]
colnames(hic2) = c('region1','region2','IF')

cells_table <- create.hic.table(hic1, hic2, chr = chr,scale = FALSE)

options(bitmapType = 'cairo')
#show parameters choices for "A" cutoff
CairoPDF(file.path(output_dir,"parameter_filter.pdf"))
filter_params(cells_table,Plot = TRUE)
dev.off()
cells_norm_table <- hic_loess(cells_table, Plot = FALSE ,parallel = TRUE)

#A cutoff can be 10th percentile
CairoPDF(file.path(output_dir,"MD_plot_1.pdf"))
hic_compare_table <- hic_compare(cells_norm_table,adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)
dev.off()

CairoPDF(file.path(output_dir,"MD_plot_2.pdf"))
MD.plot2(M = hic_compare_table$adj.M, D = hic_compare_table$D,hic_compare_table$p.value, smooth = FALSE)
dev.off()

####################################  Delta Plot  ####################################
hic1 = read.csv(file.path(output_dir,'tmp',paste0(prefix1,'.oe.txt')),sep = '\t',header = FALSE)
hic1 = hic1[c(3,6,7)]
colnames(hic1) = c('region1','region2','IF')
hic2 = read.csv(file.path(output_dir,'tmp',paste0(prefix2,'.oe.txt')),sep = '\t',header = FALSE)
hic2 = hic2[c(3,6,7)]
colnames(hic2) = c('region1','region2','IF')

cells_table <- create.hic.table(hic1, hic2, chr = chr,scale = FALSE)
resolution = cells_table$end1[1] - cells_table$start1[1]
positions_index1 = (cells_table$start1 - min(cells_table$start1)+ resolution) / resolution
positions_index1 = as.integer(positions_index1)
positions_index2 = (cells_table$start2 - min(cells_table$start2)+ resolution) / resolution
positions_index2 = as.integer(positions_index2)

hic = data.frame(pos1 = positions_index1, pos2 = positions_index2, delta = cells_table$IF2 - cells_table$IF1)

# Get upper triangle data
all_positions <- unique(c(hic$pos1, hic$pos2))
complete_grid <- expand.grid(Var1 = all_positions, Var2 = all_positions)

# add 0
melted_corr <- complete_grid %>%
  left_join(hic %>% rename(Var1 = pos1, Var2 = pos2, value = delta), 
            by = c("Var1", "Var2")) %>%
  mutate(value = ifelse(is.na(value), 0, value))

melted_corr <- hic %>%
  rename(Var1 = pos1,
         Var2 = pos2,
         value = delta)

melted_corr$value[melted_corr$value>4] = 4
melted_corr$value[melted_corr$value< -4] = -4

p1 <- ggplot(data = melted_corr, aes(Var2, Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(
    low = "#BBDEEC", 
    high = "#802520",
    limits = c(-4, 4)
  ) +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),   
    axis.title = element_blank()
  ) +
  coord_fixed()

pdf(file.path(output_dir,"delta.pdf"))
p1
dev.off()


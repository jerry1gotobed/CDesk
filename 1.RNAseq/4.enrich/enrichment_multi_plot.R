args <- commandArgs(trailingOnly = TRUE)
enrich_file <- args[1]
output = args[2]
plot_width <- as.numeric(args[3])
plot_height <- as.numeric(args[4])

#Load packages
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
library(enrichplot) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()

if (!dir.exists(output)) dir.create(output, recursive = TRUE)

result = read.csv(enrich_file,header = TRUE,sep = ",", stringsAsFactors = FALSE)
result$Sample = as.character(result$Sample)

result1 <- result %>%
  group_by(Description, ONTOLOGY) %>%
  summarise(
    count = n(),
    mean_pvalue = mean(pvalue, na.rm = TRUE)
  ) %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(count), mean_pvalue, .by_group = TRUE) %>%
  ungroup() 

result$Description = factor(result$Description,levels=rev(result1$Description))

pdf(file.path(output,"enrich_combine.pdf"),width = plot_width,height = plot_height)
if ('ONTOLOGY' %in% colnames(result)){
  result$ONTOLOGY = factor(result$ONTOLOGY,levels = c('BP','CC','MF','KEGG','Custom'))
  ggplot(data = result, mapping = aes(x = Description, y = Sample)) +
    geom_point(aes(size = Count, color = p.adjust)) +  # 图中点保持实心
    scale_size(range = c(2, 6),
               guide = guide_legend(override.aes = list(shape = 21,  # 图例显示为空心圆
                                                        fill = "white",
                                                        color = "black",
                                                        stroke = 0.8),
                                    order = 2)) +
    scale_colour_gradientn(
      colours = c("#E90F44", "#63ADEE"),
      limits = c(0, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free", switch = "y") +
    coord_flip() +
    labs(x='Pathway') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.1),
      strip.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
} else {
  ggplot(data = result, mapping = aes(x = Description, y = Sample)) +
    geom_point(aes(size = Count, color = p.adjust)) +  # 图中点保持实心
    scale_size(range = c(2, 6),
               guide = guide_legend(override.aes = list(shape = 21,  # 图例显示为空心圆
                                                        fill = "white",
                                                        color = "black",
                                                        stroke = 0.8),
                                    order = 2)) +
    scale_colour_gradientn(
      colours = c("#E90F44", "#63ADEE"),
      limits = c(0, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    coord_flip() +
    labs(x='Pathway') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.1),
      strip.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
}
while (!is.null(dev.list()))  dev.off()
cat("Thanks for using ! ^_^ , any question please contact with Bioinformatic Team of Pei\n")

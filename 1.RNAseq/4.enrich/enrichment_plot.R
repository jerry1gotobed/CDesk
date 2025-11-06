suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
enrich_file <- args[1]
output = args[2]
plot_width <- as.numeric(args[3])
plot_height <- as.numeric(args[4])

if (!dir.exists(output)) dir.create(output, recursive = TRUE)
high_color = 'red'
low_color = 'blue'

plot_bubble <- function(enrich_result,high_color,low_color) {
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
      theme(text = element_text(size = 10))
    return(p)
  }
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    p <- ggplot(enrich_result)+
      geom_point(aes(x=GeneRatio,y=Description,col=-log10(p.adjust),size=Count))+
      theme_bw()+
      scale_color_gradient(high = high_color,low = low_color)+
      theme(text = element_text(size = 10)) 
    return(p)
  }}

plot_bar <- function(enrich_result) {
  if (colnames(enrich_result)[1]=="ONTOLOGY") {
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description,fill=ONTOLOGY),stat = "identity")+
      theme_bw() +
      facet_grid(ONTOLOGY~.,scales = "free_y",space = "free_y")+
      theme(text = element_text(size = 10))
    return(p)
  }
  if (colnames(enrich_result)[1]!="ONTOLOGY") {
    p <- ggplot(enrich_result) +
      geom_bar(aes(x=-log10(p.adjust),y=Description,fill='orange'),stat = "identity",show.legend = FALSE)+
      theme_bw() +
      theme(text = element_text(size = 10))
    return(p)  
  }}

enrich_result = read.csv(enrich_file,header = TRUE,sep = ",", stringsAsFactors = FALSE)
if (length(enrich_result$Description)==0) {
  stop("No pathways meeting the criteria were found.")
}

enrich_result <- enrich_result[order(enrich_result$p.adjust), ]
enrich_result <- enrich_result[!duplicated(enrich_result$Description), ]
enrich_result$Description <- factor(enrich_result$Description,levels = as.character(enrich_result$Description[rev(order(enrich_result$p.adjust))]))
enrich_result$GeneRatio <-  sapply(enrich_result$GeneRatio, cal <- function(x){out=eval(parse(text = x)); return(out)})

file_name = basename(enrich_file)
file_prefix <- sub("\\.[^.]*$", "", file_name)

# Bubble plot
pdf(paste0(output,'/',file_prefix,"_bubble.pdf"),width = plot_width,height = plot_height)
plot_bubble(enrich_result,high_color,low_color)
while (!is.null(dev.list()))  dev.off()

# Bar plot
pdf(paste0(output,'/',file_prefix,"_bar.pdf"),width = plot_width,height = plot_height)
plot_bar(enrich_result)
while (!is.null(dev.list()))  dev.off()
cat("Thanks for using ! ^_^ , any question please contact with Bioinformatic Team of Pei\n")

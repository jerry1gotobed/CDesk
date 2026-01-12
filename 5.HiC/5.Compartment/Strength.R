suppressMessages(library(readxl))
suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggsignif))
suppressMessages(library(ggplot2))
suppressMessages(library(enrichplot))
suppressMessages(library(plotly))

calculate_ab_strength <- function(x,sample){
  sum_mtx <- read.table(paste0(x,".s_matrix.csv"),header = TRUE)
  count_mtx <- read.table(paste0(x,".c_matrix.csv"),header = TRUE)
  AA <- sum(sum_mtx[42:51,42:51])/sum(count_mtx[42:51,42:51])
  BB <- sum(sum_mtx[2:11,2:11])/sum(count_mtx[2:11,2:11])
  AB <- sum(sum_mtx[42:51,2:11])/sum(count_mtx[42:51,2:11])
  new_result <- data.frame(
    type = c("AA","BB","AB"),
    cell = rep(sample,3),
    oe = c(AA,BB,AB),
    strength = (AA*BB/(AB*AB))
  )
  return(new_result)
}

args <- commandArgs(trailingOnly = TRUE)

path <- args[1]
samples <- args[2]
tag <- args[3]
width = args[4]
height = args[5]

file_path = file.path(path,'interaction')
samples <- unlist(strsplit(samples, ","))
tag <- unlist(strsplit(tag, ","))
save_path <- file.path(path,'img')

result <- data.frame(
  type = 1,
  cell = 1,
  oe = 1,
  strength = 1
)
for (i in 1:length(samples)){
  new_result <- calculate_ab_strength(paste0(file_path,"/",samples[i]),tag[i])
  result <- rbind(result,new_result)
}
result <- result[-1, ]

color <- c("#FF0000","#FF7F00","#FFFF00","#00FF00","#0000FF","#4B0082","#8B00FF")
colo <- colorRampPalette(color, bias=1)(length(samples))

pdf(paste0(save_path,"/Compartment_transition.pdf"),width = width,height = height)
ggplot(result)+
  geom_bar(aes(x=type,y=oe,fill=cell),stat = "identity",position = "dodge")+
  theme_linedraw(base_line_size = I(0))+scale_fill_manual(values = colo)+
  theme(text = element_text(size = 20),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ylab("Average O/E contacts")+ggtitle("Compartment transition")
while (!is.null(dev.list()))  dev.off()

strength <- result[,c(2,4)]
strength <- unique(strength)

pdf(paste0(save_path,"/Compartment_strength.pdf"),width = width,height = height)
ggplot(strength)+
  geom_bar(aes(x=cell,y=strength,fill=cell),stat ="identity")+
  theme_linedraw(base_line_size = I(0))+scale_fill_manual(values = colo)+
  theme(text = element_text(size = 20),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")+
  ylab("Compartment strength")+xlab('Sample')
while (!is.null(dev.list()))  dev.off()

save_path <- file.path(path,'AB')
write.csv(result,paste0(save_path,"/Compartment_strength_oe.csv"),row.names = FALSE)

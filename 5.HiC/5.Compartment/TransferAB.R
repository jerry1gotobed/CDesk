suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

file_path <- args[1]
resolution <- args[2]
samples <- args[3]
tag = args[4]
width = args[5]
height = args[6]

samples = unlist(strsplit(samples, split = ","))
tag = unlist(strsplit(tag, split = ","))
save_path = file.path(file_path,'AB')

# Read the files
expected_cis_path = file.path(file_path,'expected_cis')
file_paths <- paste0(expected_cis_path,"/",samples,'.eigs.',as.character(resolution),'.cis.vecs.tsv')
data_list <- lapply(file_paths, function(x) {
  read.table(x, header = TRUE, fill = TRUE, na.strings = "NA")
})
names(data_list) <- tag

e1 <- data.frame(
  chr = data_list[[1]][, "chrom"],
  start = data_list[[1]][, "start"],
  end = data_list[[1]][, "end"]
)

# Get A/B from E1 column
for (samp in names(data_list)) {
  e1[[samp]] <- data_list[[samp]]$E1
}

ab <- e1

ab[,4:ncol(ab)] <- ab[,4:ncol(ab)]/abs(ab[,4:ncol(ab)])

sample_vector <- names(data_list)
ab[sample_vector] <- lapply(ab[sample_vector], function(x) {
  x <- sub("-1", "B", x)
  sub("1", "A", x)
})

ab <- ab[complete.cases(ab),]

# Get A/B transfer
abtransfer <- data.frame()
n <- length(ab)

for (i in 4:(n-1)) {
  temp_df <- data.frame(
    status = paste0(ab[[i]], ">", ab[[i+1]]),
    cmp = paste0(names(data_list)[i-3], "_vs_", names(data_list)[i-2])
  )
  abtransfer <- rbind(abtransfer, temp_df)
}

write.csv(ab,paste0(save_path,"/AB.csv"),row.names = FALSE)
abtransfer_count <- tapply(paste0(abtransfer$cmp,":",abtransfer$status),paste0(abtransfer$cmp,":",abtransfer$status),length)%>%as.data.frame()
abtransfer_count <- data.frame(info=rownames(abtransfer_count),count=abtransfer_count$.)

abtransfer_count <- separate(abtransfer_count,"info",c("cmp","status"),sep = ":")
write.csv(abtransfer_count,paste0(save_path,"/A_B_transfer_count.csv"),row.names = FALSE)

save_path = file.path(file_path,'img')
pdf(paste0(save_path,"/A_B_transfer_count.pdf"),width = width,height = height)
ggplot(abtransfer_count)+
  geom_bar(aes(x=cmp,fill=status,y=count),position="stack",stat = "identity")+
  geom_text(aes(x=cmp,label=count,y=count),position = "stack")+
  coord_polar(theta = "y")+
  theme_linedraw(base_line_size = I(0))+xlab("")+ylab("")
while (!is.null(dev.list()))  dev.off()

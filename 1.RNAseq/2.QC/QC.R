###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
sys_argv <- commandArgs(T)
input_file <- sys_argv[1]
output_directory <- sys_argv[2]
group_file <- sys_argv[3]
batch <- sys_argv[4]
plot_width <- as.numeric(sys_argv[5])
plot_height <- as.numeric(sys_argv[6])

suppressMessages(library(dplyr))

library("gplots") %>% suppressMessages()
library(stringr) %>% suppressMessages()
library(data.table) %>% suppressMessages()
library(ggplot2) %>% suppressMessages()
library(ggrepel) %>% suppressMessages()
suppressMessages(library(edgeR))
suppressMessages(library(sva))

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
# Read data
data <- fread(input_file, data.table = FALSE)
rownames(data) = data[[1]];data = data[,-1]

meta = fread(group_file, data.table = FALSE)
if (batch!='no'){
  check_col_names = colnames(meta)
  if (any(!sort(c("sample","group","tag","batch" )) == sort(check_col_names))){
    stop('Should have 4 columns in the sample information file：sample,group,tag,batch')
  }
}else {
  check_col_names = colnames(meta)
  if (any(!sort(c("sample","group","tag")) == sort(check_col_names[1:3]))){
    stop('Should have 3 columns in the grouping file：sample,group,tag')
  }
}

meta$group = as.character(meta$group)
meta$tag = as.character(meta$tag)

if (any(!meta$sample %in% colnames(data))){
  error_message <- paste0(
    "Can not find the sample in the sample information file\n",
    paste(meta$sample[!meta$sample %in% colnames(data)], collapse = ", "), "\n"
  )
  stop(error_message)
}

# Sort 
data = data[,colnames(data) %in% meta$sample]
data = data[,meta$sample]

data = log2(data+1)

if (batch!='no'){
  if (!"batch" %in% colnames(meta)) {
    stop("Error: 'batch' column is missing in the meta data!")
  }
  # Check NA
  if (any(is.na(meta$batch))) {
    stop("Error: 'batch' column contains missing values!")
  }
  # Check empty
  if (length(meta$batch) == 0 || all(meta$batch == "")) {
    stop("Error: 'batch' column is empty!")
  }
  if (length(unique(meta$batch)) == 1) {
    stop("Error: 'More than 1 batch")
  }
  if (batch=='removeBatchEffect'){
    experiment_batch <- as.factor(meta$batch)
    data <- removeBatchEffect(data,batch=experiment_batch)
  }else if (batch=='combat'){
    experiment_batch <- as.factor(meta$batch)
    data = ComBat(dat = data, batch = experiment_batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  }
}

# Process
row.names(meta) = meta$sample
all_samples <- meta[colnames(data),'tag']

###################################################################################################
################################                PLOT               ################################
###################################################################################################
options(scipen = 200)
data = na.omit(data)
# PCA
pdf(paste0(output_directory,"/pca.pdf"),width=plot_width,height=plot_height)
pca <- prcomp(t(data))

pca_data <- data.frame(
  PC1 = pca$x[, "PC1"],
  PC2 = pca$x[, "PC2"],
  Samples = all_samples
)
# Plot
ggplot(pca_data, aes(x = PC1, y = PC2, label = Samples, color = Samples)) +
  geom_point(size = 3) +  
  geom_text_repel(size = 3, max.overlaps = Inf,segment.linetype = "dashed") +  
  labs(title = "PCA", x = "PC1", y = "PC2") +
  theme_classic()
while (!is.null(dev.list()))  dev.off()


# MDS
d <- dist(t(data)) 
fit <- cmdscale(d, eig = TRUE, k = 2)
x <- fit$points[,1]
y <- fit$points[,2]

pdf(paste0(output_directory,"/MDS.pdf"), width = plot_width, height = plot_height)
mds_data <- data.frame(
  c1 = x,
  c2 = y,
  Samples = all_samples
)
# Plot
ggplot(mds_data, aes(x = c1, y = c2, label = Samples, color = Samples)) +
  geom_point(size = 3) +  
  geom_text_repel(size = 3, max.overlaps = Inf,segment.linetype = "dashed") + 
  labs(title = "MDS", x = "Coordinate 1", y = "Coordinate 2") +
  theme_classic()

while (!is.null(dev.list()))  dev.off()

# correlation
colnames(data) = meta$tag
sample_cor <- cor(data)
pdf(paste0(output_directory,"/cor.pdf"),width=plot_width,height=plot_height)
ColorRamp <- colorRampPalette(c("white","#F2FAB0","red"), bias=1)(100)
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",margins=c(13,13), cexRow=1, cexCol=1, revC=T, symm=T, distfun=function(c) as.dist(1 - c))
cor_cell <- matrix(as.character(round(sample_cor, 2)), ncol=dim(sample_cor)[2])
heatmap.2(sample_cor,main="",col=ColorRamp,key=F,trace="none",cellnote=cor_cell, notecol="black", notecex=0.5, margins=c(13,13), cexRow=1, cexCol=1, revC=T, symm=T, distfun=
            function(c) as.dist(1 - c))
while (!is.null(dev.list()))  dev.off()

pdf(paste0(output_directory,"/cor.legend.pdf"),width=3,height=1.8)
par(mar=c(4,2,4,2))
plotMatrix <- sample_cor
all_exp <- as.matrix(plotMatrix) # using same scale bar
zmax <- max(na.omit(all_exp))
zmin <- min(na.omit(all_exp))
ColorLevels <- seq(to=zmax,from=zmin, length=1000)   #number sequence
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="Correlation cofficient",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T);
box(lwd=2)
axis(side=1,c(zmin,(zmax+zmin)/2,zmax),labels=c(round(zmin,2),round((zmax+zmin)/2,2),round(zmax,2)))
#dev.off()
while (!is.null(dev.list()))  dev.off()
cat('Done, you can check the results now\n')

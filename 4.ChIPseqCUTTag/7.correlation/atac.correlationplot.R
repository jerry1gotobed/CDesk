#! /usr/bin/env Rscript

# Parameters
args <- commandArgs(trailingOnly = TRUE)

# Check
if (length(args) != 1) {
  stop("Usage: Rscript script.R <input_file>")
}

# Read input
input_file <- args[1]

suppressMessages(library(ggplot2),classes = c("message","warning"))
suppressMessages(library(dplyr),classes = c("message","warning"))
suppressMessages(library(tidyr),classes = c("message","warning"))

data <- read.table(input_file, header = TRUE)

# Filter data
data <- data[rowSums(data[, 3:4]) > 0, ]

# Get sample names
sample1 <- colnames(data)[4]
sample2 <- colnames(data)[5]

# Rename columns
colnames(data)[c(4, 5)] <- c("s1", "s2")

# Calculate x and y axis limit
limitx <- quantile(data$s1, probs = c(0.0001, 0.9999))
limity <- quantile(data$s2, probs = c(0.0001, 0.9999))
axis_max <- max(limitx[2],limity[2])
axis_min <- min(limitx[1],limity[1])

# Calculate Spearman and Pearson correlation coefficient
spcor <- cor(data$s1, data$s2, method = "spearman") %>% round(digits = 2)
pscor <- cor(data$s1, data$s2, method = "pearson") %>% round(digits = 2)

# Output
output_file <- gsub(".txt$", "", input_file)

# Plot
pdf(paste0(output_file,".pdf"), width = 10, height = 8)
ggplot(data) +
  geom_point(aes(x = s1, y = s2, col = I("DeepSkyBlue"), alpha = I(0.02), size = I(10))) +
  geom_abline(linetype = "dashed", linewidth = 0.2, alpha = 0.5, intercept = 0, slope = 1) +
  theme_bw() +
  xlab(sample1) +
  ylab(sample2) +
  xlim(axis_min,axis_max) +
  ylim(axis_min,axis_max) +
  coord_fixed(ratio = 1) +
  annotate("text", x = axis_max * 0.25, y = axis_max * 0.9, size = 10, label = paste0("Pearson Cor=", pscor)) +
  annotate("text", x = axis_max * 0.25, y = axis_max*0.8,size=10,label=paste0("Spearman Cor=",spcor))+
  theme(text = element_text(size = 25))
  #+ggtitle(output_file)
while (!is.null(dev.list()))  dev.off()

# Done
cat("Plot saved to", paste0(output_file,".pdf"), "\n")

#scaled_data

data$s1 <- log2(data$s1+0.1) %>% scale()
data$s2 <- log2(data$s2+0.1) %>% scale()

colnames(data)[c(4, 5)] <- c("s1", "s2")

limitx <- quantile(data$s1, probs = c(0.0001, 0.9999))
limity <- quantile(data$s2, probs = c(0.0001, 0.9999))
axis_max <- max(limitx[2],limity[2])
axis_min <- min(limitx[1],limity[1])

spcor <- cor(data$s1, data$s2, method = "spearman") %>% round(digits = 2)
pscor <- cor(data$s1, data$s2, method = "pearson") %>% round(digits = 2)

pdf(paste0(output_file,"scaled.pdf"),width = 10,height = 8)

ggplot(data)+geom_point(aes(x=s1,y=s2,col=I("DeepSkyBlue"),alpha=I(0.02),size=I(10)))+
  theme_bw()+
  xlab(sample1)+ylab(sample2)+
  xlim(limitx) +
  ylim(limity) +
  annotate("text",x=limitx[2]*0.2+limitx[1],y=limity[2]*0.9,size=10,label=paste0("Pearson Cor=",pscor))+
  annotate("text",x=limitx[2]*0.2+limitx[1],y=limity[2]*0.8,size=10,label=paste0("Spearman Cor=",spcor))+
  theme(text = element_text(size = 25))
  #+ggtitle(paste0(output_file,"_scaled"))

while (!is.null(dev.list()))  dev.off()

cat("Plot saved to", paste0(output_file,"scaled.pdf"), "\n")

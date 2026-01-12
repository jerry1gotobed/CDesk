###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
suppressMessages(library(dplyr))
library(ggplot2) %>% suppressMessages()
###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
sys_argv <- commandArgs(T)
input_file <- sys_argv[1]
out_dir <- sys_argv[2]
###################################################################################################
################################                PLOT               ################################
###################################################################################################
temp = basename(input_file)
name = sub("\\..*$", "", temp)
pdf(paste0(out_dir,'/',name,".FragmentsLength.pdf"),width=5,height=5)
fragLen = read.table(input_file, header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)))

fragLen %>% ggplot(aes(x = fragLen, y = fragCount)) +
  geom_line(linewidth = 1) +
  #scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))
# text(x=50,y=max(data[,2])*1.1,paste(ratio,"%\n(",part,")",sep=""),col=cccol[2])
while (!is.null(dev.list()))  dev.off()

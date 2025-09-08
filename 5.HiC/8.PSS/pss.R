suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

sys_argv <- commandArgs(trailingOnly = TRUE)

# Required parameters
FileList= strsplit(sys_argv[1],',')[[1]]
SampleName = strsplit(sys_argv[2],',')[[1]]
outdir=sys_argv[3]
# Optional parameters
# 每一个数量级之间统计时划分的bin数量
nBin=as.integer(sys_argv[4])
# 绘制曲线时均值下界（不含）
GlobalCurve.Mean=as.numeric(sys_argv[5])
# 绘制曲线时bin索引数下界（不含），即前X个bin不要绘制
GlobalCurve.Bin.Min=as.numeric(sys_argv[6])
# 曲线图X轴截断值（含）
GlobalCurve.X.min=as.numeric(sys_argv[7])
GlobalCurve.X.max=as.numeric(sys_argv[8])
# 绘制与【期望】曲线差异时的Y轴偏移量：
ExpCurve.Y.bais=as.numeric(sys_argv[9])
# 绘制与【期望】曲线差异时的X轴舍弃值（不含）：
ExpCurve.X.min=as.numeric(sys_argv[10])# log(dis)> Curve.X.min ?
ExpCurve.X.max=as.numeric(sys_argv[11])# log(dis)< Curve.X.max ?
# 绘制与马赛克图的X轴范围（不含）：
Masc.X.min=as.numeric(sys_argv[12])# Global.Obs2Exp$dis>4.8 
Masc.X.max=as.numeric(sys_argv[13])# Global.Obs2Exp$dis<7.5
width=as.numeric(sys_argv[14])
height=as.numeric(sys_argv[15])

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)  # 这将创建所有必要的父目录
}

ProbStat<-function(tmp){# 两列的dataframe，第一列为距离(bin)，第二列为count
  tmp=data.frame(dis=tmp$V1, count=tmp$V2)
  tmp=as.data.frame(tapply(tmp$count, tmp$dis, sum))
  colnames(tmp)="sum"
  tmp=data.frame(dis=as.integer(rownames(tmp)), sum=tmp$sum)
  tmp$sum=tmp$sum/sum(tmp$sum)
  return(tmp)
}

OriDataList=replicate(length(SampleName) ,NA, simplify=FALSE)
NorDataList=OriDataList
GlobalData =data.frame(
  dis=as.integer(),
  sum=as.numeric(),
  sample=as.character() )

for(ii in 1:length(SampleName)){
  OriDataList[[ii]]=as.data.frame(
    data.table::fread(FileList[ii],sep="\t"))[,c('dist','contact')]
  colnames(OriDataList[[ii]])= c('V1','V2')
  NorDataList[[ii]]=ProbStat(OriDataList[[ii]])
  GlobalData=rbind.data.frame(GlobalData,
                              data.frame(NorDataList[[ii]],sample=SampleName[ii]))
}

# 1. 统计线图：距离取log10 → 线性区间平均
Global.LogBin=GlobalData
Global.LogBin$dis=log10(Global.LogBin$dis+1)
Global.LogBin$bin=round(Global.LogBin$dis*nBin)
Global.LogBin=Global.LogBin %>% group_by(sample,bin) %>%
  summarise(mean=mean(sum), .groups='drop')

pdf(paste0(outdir,"/1.global_distribution_cruve.pdf"),width=width, height=height)
ggplot(data=
         Global.LogBin[(Global.LogBin$bin>GlobalCurve.Bin.Min)&
                         (Global.LogBin$mean>GlobalCurve.Mean),])+
  geom_line(aes(x=bin/nBin, col=sample, group=sample, y=mean, size=I(1)), se=FALSE)+
  theme_classic()+
  annotation_logticks(sides = 'b',outside = TRUE)+
  theme(text=element_text(size = 25))+
  scale_y_continuous(trans = "log10")+
  scale_color_brewer(palette = "Set3")+
  xlab("Genome distance (log10bp)")+ylab("Contact probability")+
  xlim(GlobalCurve.X.min,GlobalCurve.X.max)
while (!is.null(dev.list()))  dev.off()

# 分布线图vs【标准】线图：
Global.Obs2Exp=GlobalData
Global.Obs2Exp$dis[Global.Obs2Exp$dis==0] = 1
Global.Obs2Exp$dis = log10(Global.Obs2Exp$dis)
Global.Obs2Exp$sum = log10(Global.Obs2Exp$sum)+ExpCurve.Y.bais

# 概率与【期望】概率的差值
Global.Obs2Exp$AbsDiff=(abs(Global.Obs2Exp$sum)-abs(Global.Obs2Exp$dis))
if(sum(Global.Obs2Exp$AbsDiff>0)>nrow(Global.Obs2Exp)*0.5){
  Global.Obs2Exp$AbsDiff=-Global.Obs2Exp$AbsDiff
}
# 概率相比【期望】概率偏差的占比
Global.Obs2Exp$Deltas=(abs(Global.Obs2Exp$dis)/abs(Global.Obs2Exp$sum))
# 删除X轴上较小的区域！！！！！！
Global.Obs2Exp=Global.Obs2Exp[Global.Obs2Exp$dis>GlobalCurve.Bin.Min,]

# 输出前述的两个图：
pdf(paste0(outdir,"/2.obs_to_exp_absdiff.pdf"),width = width,height = height)
ggplot(Global.Obs2Exp[
  (Global.Obs2Exp$dis>ExpCurve.X.min)&
    (Global.Obs2Exp$dis<ExpCurve.X.max),])+
  geom_line(aes(x=dis, y=AbsDiff, col=sample, group=sample))+
  annotation_logticks(sides = 'b',outside = TRUE)+
  scale_color_manual(values =brewer.pal(n = length(SampleName), name = "Set3"))+
  theme_classic()+
  theme(text = element_text(size = 20))+
  xlab("Genome distance (log10bp)")+ylab("Difference to expectation (log10bp)")
while (!is.null(dev.list()))  dev.off()

pdf(paste0(outdir,"/3.obs_to_exp_delta.pdf"), width =width,height = height)
ggplot(Global.Obs2Exp[
  (Global.Obs2Exp$dis>ExpCurve.X.min)&
    (Global.Obs2Exp$dis<ExpCurve.X.max),])+
  geom_smooth(aes(x=dis,y=Deltas,col=sample,group=sample), method=loess,se=FALSE)+
  scale_color_manual(values =brewer.pal(n = length(SampleName), name = "Set3"))+
  theme_classic()+
  theme(text = element_text(size = 25))+
  xlab("Genome distance (log10bp)")+ylab("Deviation ratio to expectation")
while (!is.null(dev.list()))  dev.off()

Global.Sample.Mosaic=Global.Obs2Exp[
  (Global.Obs2Exp$dis>Masc.X.min)&
    (Global.Obs2Exp$dis<Masc.X.max),]
Global.Sample.Mosaic$bin=round(Global.Sample.Mosaic$dis*nBin)
Global.Sample.Mosaic=Global.Sample.Mosaic %>% group_by(sample,bin) %>%
  summarise(mean_AbsDiff=mean(AbsDiff), 
            mean_Deltas=mean(Deltas),
            .groups='drop')
Global.Sample.Mosaic$bin=as.integer(Global.Sample.Mosaic$bin)

Glo.Sam.Mos.Scaled=Global.Sample.Mosaic
# GlobalData =data.frame(
#   dis=as.integer(),
#   sum=as.numeric(),
#   sample=as.character() )
Scale2Zero<-function(x){
  y=scale(x)
  y=y-min(y)}

Glo.Sam.Mos.Scaled=Global.Sample.Mosaic %>% group_by(sample) %>%
  mutate(scaled1 = Scale2Zero(mean_AbsDiff))
colnames(Glo.Sam.Mos.Scaled)[5]="ScaledLog1"
Glo.Sam.Mos.Scaled=Glo.Sam.Mos.Scaled %>% group_by(sample) %>%
  mutate(scaled2 = Scale2Zero(mean_Deltas))
colnames(Glo.Sam.Mos.Scaled)[6]="ScaledLog2"
Glo.Sam.Mos.Scaled=as.data.frame(Glo.Sam.Mos.Scaled)
Glo.Sam.Mos.Scaled$bin=Glo.Sam.Mos.Scaled$bin/nBin
Glo.Sam.Mos.Scaled$sample=
  factor(Glo.Sam.Mos.Scaled$sample,levels = SampleName[length(SampleName):1])

pdf(paste0(outdir,"/4.Heatmap_absdiff.pdf"),width = width,height = height)
ggplot(Glo.Sam.Mos.Scaled)+geom_tile(aes(x=bin,y=sample,fill=ScaledLog1))+
  theme_classic()+
  annotation_logticks(sides = 'b',outside = TRUE)+
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(vjust = 1))+
  scale_fill_distiller(palette = "RdBu",)+
  xlab("Genome distance (log10bp)")+ylab("Sample")
while (!is.null(dev.list()))  dev.off()

pdf(paste0(outdir,"/5.Heatmap_delta.pdf"),width = width,height = height)
ggplot(Glo.Sam.Mos.Scaled)+geom_tile(aes(x=bin,y=sample,fill=ScaledLog2))+
  theme_classic()+
  annotation_logticks(sides = 'b',outside = TRUE)+
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(vjust = 1))+
  scale_fill_distiller(palette = "RdBu",)+
  xlab("Genome distance (log10bp)")+ylab("Sample")
while (!is.null(dev.list()))  dev.off()

cat('\nDone, you can see the results now\n')

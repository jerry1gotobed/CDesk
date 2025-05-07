import argparse
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import gridspec
import omicverse as ov
import subprocess

def deseqMode(args):
    # 获取输入
    count_file = args.input
    output_dir = args.output_dir
    meta_file = args.meta_file
    plot = args.plot
    goi = args.gene_of_interest.split(',')
    goi = [x for x in goi if x != '']
    top_genes = int(args.top_genes)
    fc_threshold = args.fc_threshold
    pval_threshold = args.pval_threshold
    
    # 检查目录是否存在，不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Make folder: {output_dir}")
    else:
        print(f"Folder exists: {output_dir}")
    
    # 把count和fpkm等矩阵存储到输出目录
    count_df = pd.read_csv(count_file, index_col=0, sep=',', header=0)
    count_df.to_csv(os.path.join(output_dir,'count.csv'),sep=',')
    
    # 读取meta文件，里面存放分组信息，并分析
    metas = readMetaFile(meta_file)
    if plot:
        draw_args = vars(args).copy()
        draw_args['input_dir']=args.output_dir

    for meta in metas:
        result = deseqAnalysis(count_file, meta, output_dir, plot=plot,fc_threshold=fc_threshold,pval_threshold=pval_threshold)
        gene_list = list(set(result.sort_values(by="qvalue").index[:top_genes]) | set(goi))
        # 如果要画图，那么就画图
        if plot:
            gene_list = list(set(result.sort_values(by="qvalue").index[:top_genes]) | set(goi))
            drawHeatmapPlot(gene_list=gene_list, meta=meta, **draw_args)
            gene_list = list(set(goi))
            drawvolcanoPlot(gene_list=goi,meta=meta,**draw_args)
    return

def deseqAnalysis(count_file:str, meta:tuple, output_dir:str, plot:bool,fc_threshold=1,pval_threshold:float=0.05,**kwargs):
    
    # ---解析meta文件-----------------------------------------
    (name, treatment_groups, control_groups) = meta
    
    # ---读取count文件，转换为int格式--------------------------
    data=pd.read_csv(count_file,index_col=0,sep=',',header=0)
    data = data.astype(int)
    data_filtered = data[treatment_groups+control_groups]
    data_filtered = data_filtered[data_filtered.sum(axis=1)>=1]
    
    # ---deseq2做分析-----------------------------------------
    
    # 生成deg对象，并除重复index
    dds=ov.bulk.pyDEG(data_filtered)
    dds.drop_duplicates_index()
    print('... drop_duplicates_index success')
    
    # 调用desseq2
    result = dds.deg_analysis(treatment_groups,control_groups,
            method='DEseq2')
    
    # 过滤低表达基因
    #result = result.loc[result['log2(BaseMean)']>1]
    foldchange=float(fc_threshold)
    pval_threshold=float(pval_threshold)
    fc_max,fc_min=foldchange,0-foldchange
    result['sig']='normal'
    result.loc[((result["log2FC"]>fc_max)&(result['qvalue']<pval_threshold)),'sig']='up'
    result.loc[((result["log2FC"]<fc_min)&(result['qvalue']<pval_threshold)),'sig']='down'
    # 结果存储到本地
    result.to_csv(os.path.join(output_dir,name+'.csv'),sep=',')
    
    return result

def TMode(args):
    
    # 获取输入
    count_file = args.input
    output_dir = args.output_dir
    meta_file = args.meta_file
    plot = args.plot
    goi = args.gene_of_interest.split(',')
    goi = [x for x in goi if x != '']
    top_genes = int(args.top_genes)
    fc_threshold = args.fc_threshold
    pval_threshold = args.pval_threshold
    
    # 检查目录是否存在，不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Make folder: {output_dir}")
    else:
        print(f"Folder exists: {output_dir}")
    
    # 把count和fpkm等矩阵存储到输出目录
    count_df = pd.read_csv(count_file, index_col=0, sep=',', header=0)
    count_df.to_csv(os.path.join(output_dir,'count.csv'),sep=',')
 
    # 读取meta文件，里面存放分组信息，并分析
    metas = readMetaFile(meta_file)
    if plot:
        draw_args = vars(args).copy()
        draw_args['input_dir']=args.output_dir

    for meta in metas:
        result = TAnalysis(count_file, meta, output_dir, plot=plot,fc_threshold=fc_threshold,pval_threshold=pval_threshold)
        gene_list = list(set(result.sort_values(by="qvalue").index[:top_genes]) | set(goi))
        # 如果要画图，那么就画图
        if plot:
            gene_list = list(set(result.sort_values(by="qvalue").index[:top_genes]) | set(goi))
            drawHeatmapPlot(gene_list=gene_list, meta=meta, **draw_args)
            gene_list = list(set(goi))
            drawvolcanoPlot(gene_list=gene_list,meta=meta,**draw_args)
    return

def TAnalysis(count_file:str, meta:tuple, output_dir:str, plot:bool,fc_threshold=1,pval_threshold:float=0.05, **kwargs):
    
    # ---解析meta文件-----------------------------------------
    (name, treatment_groups, control_groups) = meta
    
    # ---读取count文件，转换为int格式--------------------------
    data=pd.read_csv(count_file,index_col=0,sep=',',header=0)
    data = data.astype(int)
    data_filtered = data[treatment_groups+control_groups]
    data_filtered = data_filtered[data_filtered.sum(axis=1)>=1]
    
    # ---deseq2做分析-----------------------------------------
    
    # 生成deg对象，并除重复index
    dds=ov.bulk.pyDEG(data_filtered)
    dds.drop_duplicates_index()
    print('... drop_duplicates_index success')
    
    # 调用deseq2
    result = dds.deg_analysis(treatment_groups,control_groups,
            method='ttest',multipletests_method='bonferroni')
    
    # 过滤低表达基因
    #result = result.loc[result['log2(BaseMean)']>1]
    foldchange=float(fc_threshold)
    pval_threshold=float(pval_threshold)
    fc_max,fc_min=foldchange,0-foldchange
    result['sig']='normal'
    result.loc[((result["log2FC"]>fc_max)&(result['qvalue']<pval_threshold)),'sig']='up'
    result.loc[((result["log2FC"]<fc_min)&(result['qvalue']<pval_threshold)),'sig']='down'
    # 结果存储到本地
    result.to_csv(os.path.join(output_dir,name+'.csv'),sep=',')
        
    return result

def gfoldMode(args): 
    # 获取参数
    species = getGTF(args.species,args.config_path)
    input_dir = args.input_dir
    output_dir = args.output_dir
    meta_file = args.meta_file
    plot = args.plot
    goi = args.gene_of_interest.split(',')
    goi = [x for x in goi if x != '']
    top_genes = int(args.top_genes)
    fc_threshold = args.fc_threshold
    pval_threshold = args.pval_threshold
    
    # 将输入目录中的所有bam文件存入list
    bam_files = [os.path.join(input_dir,f) for f in os.listdir(input_dir) if f.endswith('.bam')]

    # 检查输出目录是否存在，不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Make folder: {output_dir}")
    else:
        print(f"Folder already exists: {output_dir}")
    if not os.path.exists(os.path.join(output_dir,"readCnt")):    
        os.makedirs(os.path.join(output_dir,"readCnt"))
    if not os.path.exists(os.path.join(output_dir,"diff")):
        os.makedirs(os.path.join(output_dir,"diff"))

    # 假设 .gtf 文件存在，输出加载信息
    if species:
        print(f"Species is specified as {species}, GTF file will be used for mapping gene IDs to gene names")
    else:
        print(f"Species file not exist!")
        return
    
    # 读取meta文件，里面存放分组信息   
    metas = readMetaFile(meta_file)
    
    # bam文件转count matrix
    gfoldBam2Cnt(bam_files=bam_files, species=species, output_dir=output_dir)
    
    if plot:
        draw_args = vars(args).copy()
        draw_args['input_dir']=args.output_dir
    
    # 做差异分析
    for meta in metas:
        result = gfoldAnalysis(meta=meta, **vars(args))
        gene_list = list(set(result['GFOLD'].abs().nlargest(top_genes).index) | set(goi))
        
    # 画图
        if plot:
            gene_list = list(set(result['GFOLD'].abs().nlargest(top_genes).index) | set(goi))
            drawHeatmapPlot(gene_list=gene_list, meta=meta, **draw_args)
            gene_list = list(set(goi))
            drawvolcanoPlot(gene_list=gene_list,meta=meta,fc_type='GFOLD', **draw_args)    
    
    return

def gfoldBam2Cnt(bam_files:list, species:str, output_dir:str):
    # GFOLD方法需要大量调用命令行中的GFOLD
    
    # ---首先将bam转换为 read_cnt文件----------------------------------------

    processes = []
    for file in bam_files:
        name = file.split("/")[-1].split(".")[0]
        # command = f"samtools view -@ 20 {output_path}/bam/{file} | gfold count -ann ${gtf} -tag stdin -o ${name}.read_cnt"
        command = f"samtools view -@ 20 {file} | gfold count -ann {species} -tag stdin -o {os.path.join(output_dir, 'readCnt', name + '.read_cnt')}"

        processes.append(subprocess.Popen(command, shell=True))
    
    # ---等待所有进程完成，将所有结果整合入count_matrix方便后续分析--------------
    
    for p in processes:
        p.wait()
    
    # 初始化一个空的 DataFrame
    combined_count_df = pd.DataFrame()
    combined_fpkm_df = pd.DataFrame()
    for file in os.listdir(os.path.join(output_dir, 'readCnt')):
        if file.endswith('.read_cnt'):
            name = file[:-9]
            
            # 读取文件
            df = pd.read_csv(os.path.join(output_dir, 'readCnt',file), header=None, index_col=0, sep='\t')
            count_df = df[2]
            count_df.name = name
            fpkm_df = df[4]
            fpkm_df.name = name
            
            # 合并到总的 DataFrame
            combined_count_df = pd.concat([combined_count_df, count_df], axis=1)
            combined_fpkm_df = pd.concat([combined_fpkm_df, fpkm_df], axis=1)
    
    # 存储到output
    combined_count_df.to_csv(os.path.join(output_dir,'count.csv'),sep=',')
    combined_fpkm_df.to_csv(os.path.join(output_dir,'fpkm.csv'),sep=',')

def gfoldAnalysis(meta:tuple, species:str, output_dir:str, fc_threshold=1,pval_threshold:float=0.05,**kwargs):
    
    (name, treatment_groups, control_groups) = meta
    # ---基于meta分组信息进行差异分析------------------------------------------
    original_directory = os.getcwd()
    os.chdir(os.path.join(output_dir,"readCnt"))
    command = "gfold diff \
    -s1 {} \
    -s2 {} \
    -suf .read_cnt \
    -o ../diff/{}.diff \
    > ../diff/{}.gfold.log".format(
        ",".join(control_groups),
        ",".join(treatment_groups),
        name, name
        )
    subprocess.run(command, shell=True)
    
    os.chdir(original_directory)
    # ---统计结果生成result文档，方便后续调用绘图------------------------------
    result = deseqAnalysis(os.path.join(output_dir,'count.csv'), meta, output_dir,plot=False)
    
    df = pd.read_csv("{}/diff/{}.diff".format(output_dir,name), sep="\t", comment="#",header=None, index_col=1)
    df = df[[2,3,4]]
    df.columns = ["GFOLD", "E-FDR", "log2fdc"]
    result = pd.concat([result,df], join='inner', axis=1)
    foldchange=float(fc_threshold)
    pval_threshold=float(pval_threshold)
    fc_max,fc_min=foldchange,0-foldchange
    result['sig']='normal'
    result.loc[((result["GFOLD"]>fc_max)&(result['qvalue']<pval_threshold)),'sig']='up'
    result.loc[((result["GFOLD"]<fc_min)&(result['qvalue']<pval_threshold)),'sig']='down'
    #result = result.reset_index()  # 将索引转化为普通列
    #result.rename(columns={'index': 'gene_name'}, inplace=True)  # 重命名索引列
    result.index.name = 'gene_name'
    result.to_csv(os.path.join(output_dir,name+'.csv'),sep=',')
    
    return result

def drawHeatmapPlot(input_dir:str,gene_list:list, meta:tuple, output_dir:str,**kwargs):
    
    heatmap_types = ['count','fpkm']
    
    # 动态设置子图布局
    num_plots = len(heatmap_types)
    cols = 3  # 设置列数
    rows = (num_plots + cols - 1) // cols  # 计算行数

    fig = plt.figure(figsize=(8 * cols, 16 * rows))
    gs = gridspec.GridSpec(rows, cols)
    
    (name, treatment_groups, control_groups) = meta
    # 绘制每个 heatmap
    for i, heatmap_type in enumerate(heatmap_types):
        
        data_path = os.path.join(input_dir,heatmap_type+'.csv')
        if not os.path.exists(data_path):
            #raise ValueError(f'{heatmap_type} data not exists!')
            continue
        data = pd.read_csv(data_path, sep=",", header=0, index_col=0)
        missing_genes = [gene for gene in gene_list if gene not in data.index]
        if missing_genes:
            print(f"这个{missing_genes}基因不存在，跳过")
            gene_list = [gene for gene in gene_list if gene in data.index]
        data = data[treatment_groups+control_groups].loc[gene_list]
        
        ### log2plus变换
        data = np.log2(data+1) 

        ax = fig.add_subplot(gs[i // cols, i % cols])  # 计算子图的位置
        sns.heatmap(data, ax=ax, annot=False, cmap="RdYlBu_r", fmt="g", square=True)
        ax.set_title(f'{heatmap_type} Heatmap')

    # 显示图形
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, name+"_heatmap.pdf"), format="pdf")

def drawvolcanoPlot(input_dir:str, meta:tuple, output_dir:str,
                   pval_threshold, fc_threshold,  gene_list:list=None,fc_type='padj',**kwargs):
    
    (name, treatment_groups, control_groups) = meta
    # 绘制每个 heatmap
    #result = pd.read_csv(os.path.join(input_dir,name+'.csv'),header=0,index_col=0,sep=',')
    gene_list = ",".join(gene_list)
    R_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'drawVolcano.R')
    cmd = ['Rscript',R_script,
           os.path.join(input_dir,name+'.csv'),str(pval_threshold),str(fc_threshold),gene_list,name,output_dir,fc_type]
    subprocess.run(cmd)
        
 
def readMetaFile(meta):
    if not os.path.exists(meta):
        raise ValueError(f"Meta file not exists: {meta}")
    df = pd.read_excel(meta,header=0,index_col=0)
    result = []
    for column in df.columns[1:]:
        # 获取该列的值
        values = df[column]
        
        # 生成实验组和对照组的样本名列表
        treatment_groups = values.index[values==1].to_list()
        control_groups = values.index[values==-1].to_list()
        
        # 创建元组并添加到结果列表
        result.append((str(column), treatment_groups, control_groups))
    return result

def getGTF(species,config_path):
    with open(config_path, "r") as f:
        config = json.load(f)
    gtf = config['data'][species]['refseq_gtf']
    return gtf

def addCommonArguments(parser, func):
    
    # 添加 --save_dir / -s 参数（用于指定保存目录）
    parser.add_argument(
        '--output_dir', '-o', 
        type=str, 
        required=True,  # 这个参数是必须的
        help='Directory for saving output'
    )
    
    # 添加 --meta 参数，用于提供分组信息
    parser.add_argument(
        '--meta_file', '-m',
        type=str, 
        required=True,  # 必须提供
        help='Group informations for analysis'
    )
    
    # 添加 --plot 参数，决定是否要绘图
    parser.add_argument(
        '--plot', '-p',
        action='store_true', 
        help='Whether draw plots or not'
    )

    # 添加 --gene_of_interest 参数，决定是否要绘图
    parser.add_argument(
        '--gene_of_interest', '--goi', '-g',
        default='', 
        help='Interest gene for plot heatmap'
    )
    # 添加 --top_genes 参数，决定是否要绘图
    parser.add_argument(
        '--top_genes', '-t',
        default=30, 
        help='Number of top differential expression genes for plot heatmap'
    )
    
    # 添加 --fc_thresold 参数，决定是否要绘图
    parser.add_argument(
        '--fc_threshold',
        default=1, 
        help='fc_threshold'
    )

    # 添加 --pval_threshold参数，决定是否要绘图
    parser.add_argument(
        '--pval_threshold',
        default=0.05,
        help='pval_threshold'
    )
    
    parser.add_argument("config_path", help="Path to configuration file")
    # 指定默认分析函数
    parser.set_defaults(func=func)
      
def main():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description="Script for differential expression gene detection & analysis")
    
    # 创建子解析器
    subparsers = parser.add_subparsers(dest='mode', required=True, help="Choose a mode to run the script")


    # -------------------------------------------
    # Mode DESeq2 子解析器
    # -------------------------------------------
    
    parser_deseq2 = subparsers.add_parser('deseq2', help="Run by deseq2")
    
    # 添加 --input / -i 参数（用于接受输入文件）
    parser_deseq2.add_argument(
        '--input', '-i', 
        type=str, 
        required=True,  # 这个参数也是必须的
        help='Input count matrix for DESeq2'
    )
    
    addCommonArguments(parser_deseq2, func=deseqMode)


    # -------------------------------------------
    # Mode GFOLD 子解析器
    # -------------------------------------------
    parser_gfold = subparsers.add_parser('gfold', help="Run by gfold")
    
    parser_gfold.add_argument(
        '--input_dir', '-i', 
        type=str,
        required=True, 
        help="Input bam file directory for GFOLD")
    
    # 添加 --species / -s 参数（用于指定映射文件）
    parser_gfold.add_argument(
        '--species', '-s', 
        type=str, 
        default="mouse",  # 默认为小鼠
        help='Species for choosing .gtf for mapping gene id to gene name'
    )
    addCommonArguments(parser_gfold,func=gfoldMode)
    
    
    # -------------------------------------------
    # Mode adjusted_t 子解析器
    # -------------------------------------------
    
    parser_adjusted_t = subparsers.add_parser('adjusted_t', help="Run by adjusted_t")
    
    # 添加 --input / -i 参数（用于接受输入文件）
    parser_adjusted_t.add_argument(
        '--input', '-i', 
        type=str, 
        required=True,  # 这个参数也是必须的
        help='Input count matrix for adjusted_t'
    )

    addCommonArguments(parser_adjusted_t, func=TMode)
    
    
    # 解析命令行参数
    args = parser.parse_args()

    # 检查是否指定了子命令
    if args.mode is None:
        parser.print_help()  # 打印帮助信息
        print("\nError: Must specify a mode ('deseq2' , 'gfold' or 'adjusted_t').")
    
    # ---调用子func-------------------------------------
    args.func(args)   

if __name__ == '__main__':
    main()
    print('任务运行完成，您可以查看结果了！')

import os
import argparse
import sys
import subprocess
import json
import tempfile
#import psutil

root_dir = os.path.abspath(os.path.dirname(__file__))
config_path = os.path.join(root_dir, 'config.json')
with open(config_path, "r") as f:
    config = json.load(f)
    R_lib = config['R_lib']

# bulkRNA 1. preprocessing function
def bulkRNA_preprocess(args):
    bash_script = os.path.join(root_dir,'1.RNAseq','1.Preprocess','BulkRNAseqData.sh')
    args.i = args.i.rstrip('/');args.o = args.o.rstrip('/')
    cmd = [bash_script,'-s', args.species, '-i', args.i, '-o', args.o, '-p', str(args.thread), '-l', str(args.l),config_path]
    subprocess.run(cmd,check=True)

    if args.bw:
        bash_script = os.path.join(root_dir,'1.RNAseq','1.Preprocess','BulkRNAseqBW.sh')
        cmd = [bash_script, '-i', f"{args.o}/Bam", '-o', f"{args.o}/BW", '-p', str(args.thread),config_path]
        subprocess.run(cmd, check=True)

    if args.te:
        bash_script = os.path.join(root_dir,'1.RNAseq','1.Preprocess','BulkRNAseqTE.sh')
        if args.te_gene == '':
            cmd = [bash_script, '-i', f"{args.o}/Bam", '-o', f"{args.o}/TE", '-p', str(args.thread),'-s',args.species,config_path]
        else:
            cmd = [bash_script, '-i', f"{args.o}/Bam", '-o', f"{args.o}/TE", '-p', str(args.thread),'-s',args.species,'-t',args.te_gene,config_path]
        subprocess.run(cmd, check=True)

# bulkRNA 2. QC
def bulkRNA_QC(args):
    R_script = os.path.join(root_dir,'1.RNAseq','2.QC','QC.R')
    args.o = args.o.rstrip('/')
    cmd = ['Rscript', R_script,args.i,args.o,args.group]
    subprocess.run(cmd, check=True)

# bulkRNA 3. DEG
def bulkRNA_DEG(args):
    python_script = os.path.join(root_dir,'1.RNAseq','3.DEGD','DEGD.py')
    degd_env = config['conda_env']['degd'].rstrip('/')
    args.o = args.o.rstrip('/')
    if args.DEG_mode != 'gfold':
        if args.p:
            cmd = [
                degd_env+'/bin/python', python_script, args.DEG_mode,
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m,"-p",
                "--gene_of_interest",args.g,"--top_genes",str(args.t),
                "--fc_threshold",str(args.fc),
                "--pval_threshold",str(args.pval),config_path
            ]
            subprocess.run(cmd, check=True)
        else:
            cmd = [
                degd_env+'/bin/python', python_script, args.DEG_mode,
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m,config_path
            ]
            subprocess.run(cmd, check=True)
    elif args.DEG_mode == 'gfold':
        if args.p:
            cmd = [
                degd_env+'/bin/python', python_script, "gfold",
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m,"-p","--species",args.s,
                "--gene_of_interest",args.g,"--top_genes",str(args.t),
                "--fc_threshold",str(args.fc),"--pval_threshold",str(args.pval),config_path
            ]
            subprocess.run(cmd, check=True)
        else:
            cmd = [
                degd_env+'/bin/python', python_script, "gfold",
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m,"--species",args.s,config_path
            ]
            subprocess.run(cmd, check=True)

# bulkRNA 4. GO analysis
def bulkRNA_GO(args):
    R_script = os.path.join(root_dir,'1.RNAseq','4.GeneOntologyAnalysis','enrichment_analysis.R')
    args.o = args.o.rstrip('/')
    if args.GO_mode == 'custom':
        R_script = os.path.join(root_dir,'1.RNAseq','4.GeneOntologyAnalysis','enrichment_analysis_customer.R')
        cmd = [
            "Rscript", R_script, 
            args.i, args.custom, args.o,args.color_high, args.color_low,R_lib
        ]
        subprocess.run(cmd, check=True)
    elif args.GO_mode == 'GO':
        cmd = [
            "Rscript", R_script, 
            args.i, 'GO', args.species, str(args.t), args.o,args.color_high, args.color_low,R_lib
        ]
        subprocess.run(cmd, check=True)
    elif args.GO_mode == 'KEGG':
        cmd = [
            "Rscript", R_script, 
            args.i, 'KEGG', args.species, args.o,args.color_high, args.color_low,R_lib
        ]
        subprocess.run(cmd, check=True)

# bulkRNA 5. GSEA
def bulkRNA_GSEA(args):
    python_script = os.path.join(root_dir,'1.RNAseq','5.GSEA','gsea.py')
    if args.gsea_mode == 'analyze':
        args.o = args.o.rstrip('/')
        if args.preprocess:
            cmd = [
                'python', python_script, "gsea","--rnk",args.rnk, "-o",args.o,
                "--species",args.species,"--preprocess","-t",str(args.thread),"--permutation_num",str(args.permutation_num)]
        else:
            cmd = [
                'python', python_script, "gsea","--rnk",args.rnk, "-o",args.o,
                "--species",args.species,"-t",str(args.thread),"--permutation_num",str(args.permutation_num)]
        subprocess.run(cmd, check=True)
    elif args.gsea_mode == 'plot':
        args.i = args.i.rstrip('/')
        cmd = ['python', python_script,'plot','-i',args.i,'--term',args.term]
        subprocess.run(cmd, check=True)

# bulkRNA 6. Similarity Analysis
def bulkRNA_Similarity(args):
    R_script = os.path.join(root_dir,'1.RNAseq','6.SimilarityAnalysis','Similarity.R')
    args.o = args.o.rstrip('/')
    cmd = ['Rscript',R_script,args.sample1,args.sample2,args.group,args.o,args.gene,R_lib]
    subprocess.run(cmd, check=True)

# bulkRNA 7. WGCNA
def bulkRNA_WGCNA(args):
    R_script = os.path.join(root_dir,'1.RNAseq','7.WGCNA','WGCNA.R')
    args.o = args.o.rstrip('/')
    cmd = ['Rscript', R_script,args.count_file,args.pheno_file,args.trait,args.o,R_lib]
    subprocess.run(cmd, check=True)

# bulkRNA 8. AlternativeSplicingDetection
def bulkRNA_AlternativeSpliceDetection(args):
    args.o = args.o.rstrip('/')
    if args.Splice_mode == 'detect':
        bash_script = os.path.join(root_dir,'1.RNAseq','8.AlternativeSplicingDetection','alternative_splicing.sh')
        with tempfile.NamedTemporaryFile(mode="w", delete=True) as b1_temp,tempfile.NamedTemporaryFile(mode="w", delete=True) as b2_temp:
            # 写入内容到临时文件
            b1_temp.write(args.b1)
            b2_temp.write(args.b2)
            # 确保临时文件内容被写入磁盘
            b1_temp.flush()
            b2_temp.flush()
            # 执行脚本
            gtf_file = config['data'][args.species]['refseq_gtf']
            cmd = [bash_script, '--b1', b1_temp.name, '--b2', b2_temp.name, '--gtf', gtf_file,'--type',args.type,
                    '--length',str(args.readLength),'--thread',str(args.thread),'--od', args.o,'--config',config_path]
            subprocess.run(cmd, check=True)
    elif args.Splice_mode == 'draw':
        bash_script = os.path.join(root_dir,'1.RNAseq','8.AlternativeSplicingDetection','draw_sashimiplot.sh')
        gff3_file = config['data'][args.species]['gff3']
        cmd = [bash_script,'--b1',args.b1,'--b2',args.b2,'--gff3',gff3_file,'--r',args.interval,'--g',args.group,'--od',args.o,
                '--intron_scale',str(args.intron_s),'--exon_scale',str(args.exon_s),'--config',config_path]
        subprocess.run(cmd,check=True)

def parse_arguments():
    # CDesk ArgumentParser Object
    parser = argparse.ArgumentParser(prog='CDesk', description="CDesk multiomics pipeline")
    # 创建二级子命令
    subparsers = parser.add_subparsers(dest='omics', required=True, help="Choose an omic to run the script")

    # bulkRNA
    # 创建 bulkRNA 子命令的子解析器
    bulkRNA_parser = subparsers.add_parser('bulkRNA', help='bulkRNA-seq pipeline')
    bulkRNA_subparsers = bulkRNA_parser.add_subparsers(help='bulkRNA-seq pipeline', title='bulkRNA-seq pipeline', dest='bulkRNA_command')

    # bulkRNA: 1. preprocess
    # 添加 bulkRNA preprocessing 子命令
    bulkPreprocessing_parser = bulkRNA_subparsers.add_parser('preprocess', help='bulkRNA-seq preprocessing')
    # Required parameters
    bulkPreprocessing_parser.add_argument('-i', required=True, help='Specify the input directory of input file (fq/bam files)')
    bulkPreprocessing_parser.add_argument('-o', required=True, help='Specify the output directory')
    bulkPreprocessing_parser.add_argument('--species','-s', required=True, help='Specify the species',
                                      choices=[ "mm10","rn7","susScr11",'hg38','galGal6'])
    # Optional parameters
    bulkPreprocessing_parser.add_argument('-bw', action='store_true', help='If specified, perform a bam -> bigwig file transfer (optional)')
    bulkPreprocessing_parser.add_argument('-te', action='store_true', help='If specified, perform a TE expression analysis (optional)')
    bulkPreprocessing_parser.add_argument('--thread','-t', help='Specify number of threads (default is 8)',default=8,type=int)
    bulkPreprocessing_parser.add_argument('-l', help='Single sequencing, 2:Pair sequencing (default is 2)',default=2,type=int)
    bulkPreprocessing_parser.add_argument('--te_gene', help='Specify the TE gene txt file (optional)',default='')

    # bulkRNA: 2. QC
    # 添加 bulkRNA QC 子命令
    bulkQC_parser = bulkRNA_subparsers.add_parser('QC', help='bulkRNA-seq QC')
    # Required parameters
    bulkQC_parser.add_argument('-i', required=True, help='Specify the input gene expression file (.csv), column name as sample, row name as gene')
    bulkQC_parser.add_argument('-o', required=True, help='Specify the output directory')
    bulkQC_parser.add_argument('--group',required=True ,help='The grouping xlsx file')

    # bulkRNA: 3. DEG
    bulkDEGD_parser = bulkRNA_subparsers.add_parser('DEG', help='bulkRNA-seq differential gene expression analysis')
    bulkDEGD_subparser = bulkDEGD_parser.add_subparsers(dest='DEG_mode', required=True, help="Choose a mode to run the script")
    # DEseq2
    bulkDEGD_deseq2 = bulkDEGD_subparser.add_parser('deseq2',help='Perform DEseq2 differential gene expression analysis')
    bulkDEGD_deseq2.add_argument('-i',required=True,help='Input count matrix for DESeq2')
    bulkDEGD_deseq2.add_argument('-o',required=True,help='Output directory')
    bulkDEGD_deseq2.add_argument('-m',required=True,help='Group informations for analysis')
    bulkDEGD_deseq2.add_argument('-p',help='Whether to draw plot or not',action='store_true')
    bulkDEGD_deseq2.add_argument('-g',help='Interest gene for plot heatmap,split by , (optional)',default='')
    bulkDEGD_deseq2.add_argument('-t',help='Number of top differential expression genes for heatmap, default: 30',default=30,type=int)
    bulkDEGD_deseq2.add_argument('-fc',help='fc_threshold for vocano plot, default: 1.5',type=float,default=1.5)
    bulkDEGD_deseq2.add_argument('-pval',help='pval_threshold for vocano plot, default: 0.01',type=float,default=0.01)
    # adjusted_t
    bulkDEGD_adjust_t = bulkDEGD_subparser.add_parser('adjusted_t',help='Perform adjusted_t differential gene expression analysis')
    bulkDEGD_adjust_t.add_argument('-i',required=True,help='Input count matrix for DESeq2')
    bulkDEGD_adjust_t.add_argument('-o',required=True,help='Output directory')
    bulkDEGD_adjust_t.add_argument('-m',required=True,help='Group informations for analysis')
    bulkDEGD_adjust_t.add_argument('-p',help='Whether to draw plot or not',action='store_true')
    bulkDEGD_adjust_t.add_argument('-g',help='Interest gene for plot heatmap,split by , (optional)',default='')
    bulkDEGD_adjust_t.add_argument('-t',help='Number of top differential expression genes for heatmap, default: 30',default=30,type=int)
    bulkDEGD_adjust_t.add_argument('-fc',help='fc_threshold for vocano plot, defaullt: 1.5',type=float,default=1.5)
    bulkDEGD_adjust_t.add_argument('-pval',help='pval_threshold for vocano plot, default: 0.01',type=float,default=0.01)
    # gfold
    bulkDEGD_gfold = bulkDEGD_subparser.add_parser('gfold',help='Perform adjusted_t differential gene expression analysis')
    bulkDEGD_gfold.add_argument('-i',required=True,help='Input count matrix for DESeq2')
    bulkDEGD_gfold.add_argument('-o',required=True,help='Output directory')
    bulkDEGD_gfold.add_argument('-m',required=True,help='Group informations for analysis')
    bulkDEGD_gfold.add_argument('-s',required=True,help='Specify the species',choices=["pig","human","mm10"])
    bulkDEGD_gfold.add_argument('-p',help='Whether to draw plot or not',action='store_true')
    bulkDEGD_gfold.add_argument('-g',help='Interest gene for plot heatmap,split by , (optional)',default='')
    bulkDEGD_gfold.add_argument('-t',help='Number of top differential expression genes for heatmapi, default: 30',default=30,type=int)
    bulkDEGD_gfold.add_argument('-fc',help='fc_threshold for vocano plot, default: 1.5',type=float,default=1.5)
    bulkDEGD_gfold.add_argument('-pval',help='pval_threshold for vocano plot, default: 0.01',type=float,default=0.01)

    # bulkRNA: 4. GeneOntologyAnalysis
    # 添加 bulkRNA GO 子命令
    bulkGO_parser = bulkRNA_subparsers.add_parser('enrichment', help='bulkRNA-seq gene ontology analysis')
    bulkGO_subparser = bulkGO_parser.add_subparsers(dest='GO_mode', required=True, help="Choose a mode to run the script")
    # enrichment_analysis_customer
    bulkGO_custom = bulkGO_subparser.add_parser('custom', help='Perform the enrichment analysis based on the customization gene sets')
    bulkGO_custom.add_argument('-i', required=True, help='path of a single column file containing genes of interests(SYMBOL)')
    bulkGO_custom.add_argument('--custom','-c', required=True, help='A customer file containing two columns, every row is consisted of two factors: the customer term name and gene of interests separated by tab')
    bulkGO_custom.add_argument('-o', required=True, help='Output directory')
    bulkGO_custom.add_argument('--color_high', default='red', help='Assign a color representing the maximum –lg(p.adjust) value, default: red')
    bulkGO_custom.add_argument('--color_low', default='blue', help='Assign a color representing the minimum –lg(p.adjust) value, default: blue')
    # GO enrichment_analysis
    bulkGO = bulkGO_subparser.add_parser('GO', help='Perform the GO enrichment analysis based on the pubsihed GO/KEGG gene sets')
    bulkGO.add_argument('-i', required=True, help='path of a single column file containing genes of interests(SYMBOL)')
    bulkGO.add_argument('--species','-s', required=True, help='Choose the functional annotation database of specific species: human, rat, mouse, chicken, or pig',choices=["human", "rat","mouse","chicken","pig"])
    bulkGO.add_argument('-t', required=True, help='Define a threshold between 0~1, the GO terms will be merged into one once the gene sets matched to them exceed the threshold')
    bulkGO.add_argument('-o', required=True, help='Output directory')
    bulkGO.add_argument('--color_high', default='red', help='Assign a color representing the maximum –lg(p.adjust) value, default: red')
    bulkGO.add_argument('--color_low', default='blue', help='Assign a color representing the minimum –lg(p.adjust) value, default: blue')
    # KEGG enrichment_analysis
    bulkKEGG = bulkGO_subparser.add_parser('KEGG', help='Perform the KEGG enrichment analysis based on the pubsihed GO/KEGG gene sets')
    bulkKEGG.add_argument('-i', required=True, help='Path of a single column file containing genes of interests(SYMBOL)')
    bulkKEGG.add_argument('--species','-s', required=True, help='Choose the functional annotation database of specific species: human, rat, mouse, chicken, or pig',
                        choices=["human", "rat","mouse","chicken","pig"])
    bulkKEGG.add_argument('-o', required=True, help='Output directory')
    bulkKEGG.add_argument('--color_high', default='red', help='Assign a color representing the maximum –lg(p.adjust) value, default: red')
    bulkKEGG.add_argument('--color_low', default='blue', help='Assign a color representing the minimum –lg(p.adjust) value, default: blue')

    # bulkRNA: 5. GSEA
    bulkGSEA_parser = bulkRNA_subparsers.add_parser('GSEA', help='bulkRNA-seq GSEA analysis')
    bulkGSEA_subparser = bulkGSEA_parser.add_subparsers(dest='gsea_mode', required=True, help="Choose to analyze or plot")
    # analysis
    bulkGSEA_analyze_parser = bulkGSEA_subparser.add_parser('analyze', help='Perform the gsea analysis')
    bulkGSEA_analyze_parser.add_argument("--rnk",required=True,type=str, help="The ranked gene list file (e.g. gene_name(upper_case) log2FC(descending)) or the CDesk deg analysis result csv file")
    bulkGSEA_analyze_parser.add_argument("-o",required=True,  type=str, help="Directory to save the results")
    bulkGSEA_analyze_parser.add_argument("-s","--species",default="human",type=str, help="Specify the species, human/mouse, default: human",choices=['human','mouse'])
    bulkGSEA_analyze_parser.add_argument("--preprocess", action="store_true", help="Preprocess CDesk DEG analysis result, add this if you input CDesk DEG analysis result")   
    bulkGSEA_analyze_parser.add_argument("--permutation_num",default=1000 , type=int, help="Relevance of selected trait data, default: 1000")
    bulkGSEA_analyze_parser.add_argument("-t","--thread",default=4 , type=int, help="Relevance of selected trait data, default: 4")
    # plot
    bulkGSEA_plot_parser = bulkGSEA_subparser.add_parser('plot', help='Specify the GSEA term plot')
    bulkGSEA_plot_parser.add_argument("-i",required=True,help="The GSEA analysis output directory")
    bulkGSEA_plot_parser.add_argument("--term",required=True,help="Specify the terms to plot, separate by ,")

    # bulkRNA: 6. Similarity Analysis
    bulkSimilarity_parser = bulkRNA_subparsers.add_parser('similarity', help='bulkRNA-seq similarity analysis')
    # Required parameters
    bulkSimilarity_parser.add_argument('--sample1', required=True, help='Secify the first gene expression csv file')
    bulkSimilarity_parser.add_argument('--sample2', required=True, help='Specify the second gene expression csv file')
    bulkSimilarity_parser.add_argument('-o', required=True, help='Output directory')
    # Optional parameters
    bulkSimilarity_parser.add_argument('--group',default='NO',help='Specify the grouping xlsx file, default: NO')
    bulkSimilarity_parser.add_argument('--gene',default='ALL',help='Specify the gene list, default: ALL')

    # bulkRNA: 7. WGCNA
    bulkWGCNA_parser = bulkRNA_subparsers.add_parser('WGCNA', help='bulkRNA-seq WGCNA')
    bulkWGCNA_parser.add_argument('--count_file', required=True, help='The path to the file containing the gene expression count data.')
    bulkWGCNA_parser.add_argument('--pheno_file', required=True, help='A file path containing phenotypic information that describes the characteristics or experimental conditions of each sample')
    bulkWGCNA_parser.add_argument('--trait',required=True, help='Phenotypic information, calculating the correlation between gene modules and phenotypes')
    bulkWGCNA_parser.add_argument('-o',required=True,help='Output directory')

    # bulkRNA: 8. AlternativeSplicingDetection
    bulkSplice_parser = bulkRNA_subparsers.add_parser('splice', help='bulkRNA-seq alternative splicing detection')
    bulkSplice_subparser = bulkSplice_parser.add_subparsers(dest='Splice_mode', required=True, help="Choose a mode to run the script")
    # enrichment_analysis_customer.
    bulkSplice_enrich = bulkSplice_subparser.add_parser('detect', help='Detect the differentail splicing events by rMATs based on the bam file')
    bulkSplice_enrich.add_argument('--b1', required=True, help='Path to the b1 file containing BAM file paths separated by comma')
    bulkSplice_enrich.add_argument('--b2', required=True, help='Path to the b2 file containing BAM file paths separated by comma')
    bulkSplice_enrich.add_argument('--species','-s', required=True, help='Specify the species',choices=[ "mm10","rn7","susScr11",'hg38','galGal6'])
    bulkSplice_enrich.add_argument('-o', required=True, help='Output directory')
    bulkSplice_enrich.add_argument('--type', default='paired', help='Type of read used in the analysis, default: paired',choices=['paired','single'])
    bulkSplice_enrich.add_argument('--thread','-t', default=20,type=int, help='Number of threads, default: 20')
    bulkSplice_enrich.add_argument('--readLength', default=150,type=int, help='The length of each read, default: 150')
    # draw_sashimiplot
    bulkSplice_draw = bulkSplice_subparser.add_parser('draw', help='Draw the sashimi plot to visualize the differential splicing events based on the bam files')
    bulkSplice_draw.add_argument('--b1', required=True, help='Path to the b1 file containing BAM file paths separated by comma')
    bulkSplice_draw.add_argument('--b2', required=True, help='Path to the b2 file containing BAM file paths separated by comma')
    bulkSplice_draw.add_argument('--species','-s', required=True, help='Specify the species',choices=[ "mm10","rn7","susScr11",'hg38','galGal6'])
    bulkSplice_draw.add_argument('-o', required=True, help='Output directory')
    bulkSplice_draw.add_argument('--interval', required=True, help='Interval of interests,chromosome:orientation:start:end. example: chr11:-:4752513:4854957')
    bulkSplice_draw.add_argument('--group',required=True,help='The path to a *.gf file which groups the replicates.')
    bulkSplice_draw.add_argument('--exon_s',default=1,help='How much to scale down exons, default: 1')
    bulkSplice_draw.add_argument('--intron_s',default=5,help='How much to scale down introns, default: 5')
        
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = parse_arguments()

    # 判断并执行 bulkRNA 相关子命令
    if args.omics == 'bulkRNA':
        if args.bulkRNA_command == 'preprocess':
            bulkRNA_preprocess(args)
        elif args.bulkRNA_command == 'QC':
            bulkRNA_QC(args)
        elif args.bulkRNA_command == 'DEG':
            bulkRNA_DEG(args)
        elif args.bulkRNA_command == 'enrichment':
            bulkRNA_GO(args)
        elif args.bulkRNA_command == 'GSEA':
            bulkRNA_GSEA(args)
        elif args.bulkRNA_command == 'similarity':
            bulkRNA_Similarity(args)
        elif args.bulkRNA_command == 'WGCNA':
            bulkRNA_WGCNA(args)
        elif args.bulkRNA_command == 'splice':
            bulkRNA_AlternativeSpliceDetection(args)


if __name__ == "__main__":
    main()

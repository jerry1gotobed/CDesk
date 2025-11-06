import argparse
import os
import json
import pandas as pd
import numpy as np
import subprocess
import shutil
import chardet
import sys
from datetime import datetime

# Check tools
def is_tool_available(name):
    return shutil.which(name) is not None

tools = ['gfold','samtools']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

config_path = os.environ.get('CDesk_config')

def gfoldMode(args): 
    # Load parameters
    species = getGTF(args.species,config_path)
    output_dir = args.output_dir
    meta_file = args.meta_file
    plot = bool(args.plot)
    goi = args.gene_of_interest
    fc_threshold = float(args.fc_threshold)
    top_num = int(args.top_genes)
    thread = int(args.thread)
    width = float(args.width)
    height = float(args.height)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(os.path.join(output_dir,"readCnt")):    
        os.makedirs(os.path.join(output_dir,"readCnt"))
    if not os.path.exists(os.path.join(output_dir,"diff")):
        os.makedirs(os.path.join(output_dir,"diff"))

    # Check speices
    if not species:
        print(f"Species gtf file not exist!")
        sys.exit(1)
   
    with open(meta_file, 'rb') as f:
        result = chardet.detect(f.read(1000))  # Check encoding
        encoding = result['encoding']
    meta_test = pd.read_csv(meta_file, sep=None, engine='python', encoding=encoding,header=0)
    required = ["sample", "group","bam"]
    missing = set(required) - set(meta_test.columns)
    if missing:
        print('Need columns: sample,group,bam')
        sys.exit(1)
    for bam in meta_test['bam']:
        if not os.path.exists(bam):
            print(f'Can not find {bam}')
            sys.exit(1)
    if 'readcnt' in meta_test.columns:
        for readcnt in meta_test['readcnt']:
            if not os.path.exists(readcnt):
                print(f'Can not find {readcnt}')
                sys.exit(1)

    # Read the meta file 
    metas = readMetaFile(meta_file)
    
    # Bam -> count matrix
    gfoldBam2Cnt(meta=meta_test, species=species, output_dir=output_dir,thread=thread)
    
    # Differential analysis
    print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Start gfold analysis")
    for meta in metas:
        result = gfoldAnalysis(meta=meta,goi=goi,plot=plot,output_dir=output_dir,fc_threshold=fc_threshold,width=width,height=height)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if plot:
        print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Do MD plot and plot heatmap")
        cmd = ['Rscript',os.path.join(script_dir,'gfold_plot.R'),'heatmap',os.path.join(output_dir,'fpkm.csv'),goi,meta_file,str(top_num),output_dir,str(width),str(height)]
        subprocess.run(cmd)
    return

def gfoldBam2Cnt(species:str, output_dir:str,thread:int,meta: pd.DataFrame):
    # bam->read_cnt
    processes = []
    for idx, row in meta.iterrows():
        name = row['sample']
        output_read_cnt = os.path.join(output_dir, 'readCnt', name + '.read_cnt')
        if 'readcnt' in row.index:
            input_read_cnt = row['readcnt']
        else:
            input_read_cnt = ''

        # Check whether .read_cnt file exists，skip if exists
        if os.path.exists(output_read_cnt):
            print(f"Find {output_read_cnt}, skip generating")
            continue
        if os.path.exists(input_read_cnt):
            print(f"Find {input_read_cnt}, copy")
            os.makedirs(os.path.dirname(output_read_cnt), exist_ok=True)
            shutil.copy(input_read_cnt, output_read_cnt)
            continue
        
        file = row['bam']
        print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Tranfer from bam to readcnt format for {name}")
        cmd = f"samtools view -@ {thread} {file} | gfold count -ann {species} -tag stdin -o {output_read_cnt} > /dev/null 2>&1"
        subprocess.run(cmd, shell=True)
    
    # Initialize an empty DataFrame
    combined_count_df = pd.DataFrame()
    combined_fpkm_df = pd.DataFrame()
    for file in os.listdir(os.path.join(output_dir, 'readCnt')):
        if file.endswith('.read_cnt'):
            name = file[:-9]
            
            # Read the file
            df = pd.read_csv(os.path.join(output_dir, 'readCnt',file), header=None, index_col=0, sep='\t')
            count_df = df[2]
            count_df.name = name
            fpkm_df = df[4]
            fpkm_df.name = name
            
            # Merge
            combined_count_df = pd.concat([combined_count_df, count_df], axis=1)
            combined_fpkm_df = pd.concat([combined_fpkm_df, fpkm_df], axis=1)
    
    # Save the result
    combined_count_df.to_csv(os.path.join(output_dir,'count.csv'),sep=',')
    combined_fpkm_df.to_csv(os.path.join(output_dir,'fpkm.csv'),sep=',')

def gfoldAnalysis(meta:tuple, plot:bool,goi:str,output_dir:str,fc_threshold:float,width:float,height:float):
    (name, treatment_groups, control_groups) = meta
    # Differential analysis based on meta file
    original_directory = os.getcwd()
    os.chdir(os.path.join(output_dir,"readCnt"))

    command = "{} diff \
    -s1 {} \
    -s2 {} \
    -suf .read_cnt \
    -o ../diff/{}.diff \
    > /dev/null 2>&1".format(
        'gfold',
        ",".join(control_groups),
        ",".join(treatment_groups),
        name
        )
    subprocess.run(command, shell=True)
    
    os.chdir(original_directory)
    # Save the result
    df = pd.read_csv("{}/diff/{}.diff".format(output_dir,name), sep="\t", comment="#",header=None, index_col=0)
    df = df[[1,2,3,4,5,6]]
    df.columns = ["Gene","GFOLD", "E-FDR", "log2fdc",'RPKM1','RPKM2']
    result=df
    foldchange=float(fc_threshold)
    fc_max,fc_min=foldchange,0-foldchange
    result['sig']='normal'
    result.loc[((result["GFOLD"]>fc_max)),'sig']='up'
    result.loc[((result["GFOLD"]<fc_min)),'sig']='down'
    result.index.name = 'gene_name'
    result.to_csv(os.path.join(output_dir,name+'.csv'),sep=',')

    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Plot
    if plot:
        cmd = ['Rscript',os.path.join(script_dir,'gfold_plot.R'),'MA',os.path.join(output_dir,name+'.csv'),goi,str(fc_threshold),name,output_dir,str(width),str(height)]
        subprocess.run(cmd)
    return result
        
def readMetaFile(meta):
    if not os.path.exists(meta):
        raise ValueError(f"Meta file not exists: {meta}")

    with open(meta, 'rb') as f:
        result = chardet.detect(f.read(1000))  # Read the first 1000 bytes for detection
        encoding = result['encoding']
    # Read the meta file
    df = pd.read_csv(meta, sep=None, engine='python', encoding=encoding,header=0,index_col=0)
    result = []
    for column in df.columns:
        if column not in ['sample','group','readcnt','bam']:
            # Get the value of this column
            values = df[column]
            # Generate the list of sample names for the experimental group and the control group
            treatment_groups = values.index[values==1].to_list()
            control_groups = values.index[values==-1].to_list()
            # Create a tuple and add it to the result list
            result.append((str(column), treatment_groups, control_groups))
    return result

def getGTF(species,config_path):
    with open(config_path, "r") as f:
        config = json.load(f)
    gtf = config['data'][species]['refseq_gtf']
    return gtf


def parse_arguments():
    parser = argparse.ArgumentParser(description="Script for gfold differential expression gene detection & analysis")
    
    parser.add_argument(
        '--species', '-s', 
        type=str, 
        default="mouse", 
        help='Species for choosing .gtf for mapping gene id to gene name'
    )

    parser.add_argument(
        '--output_dir', '-o', 
        type=str, 
        required=True,  
        help='Directory for saving output'
    )
    
    parser.add_argument(
        '--meta_file', '-m',
        type=str, 
        required=True,  
        help='Group informations for analysis'
    )
    
    parser.add_argument(
        '--plot', '-p',
        action='store_true', 
        help='Whether draw plots or not'
    )

    parser.add_argument(
        '--top_genes', '-t',
        default=30,
        help='Number of top differential expression genes for plot heatmap'
    )

    parser.add_argument(
        '--gene_of_interest', '--goi', '-g',
        default='NO', 
        help='Interest gene file for plot heatmap'
    )
    
    parser.add_argument(
        '--fc_threshold',
        default=1, 
        help='fc_threshold'
    )
 
    parser.add_argument(
        '--thread',
        default=20,
        help='Thread number'
    )
    
    parser.add_argument(
        '--width',
        default=5,
        help='Plot Width'
    )

    parser.add_argument(
        '--height',
        default=5,
        help='Plot Height'
    )

    return parser.parse_args()
      
def main():
    args = parse_arguments()
    gfoldMode(args)

if __name__ == '__main__':
    main()
    print('Done, you can check the results now.')

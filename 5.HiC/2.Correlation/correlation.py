# 2. Correlation

import fanc
import subprocess
import pandas as pd
import sys
from datetime import datetime
import os
from functools import reduce
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import shutil
import argparse

# Parameters
parser = argparse.ArgumentParser(description="HiC sample correlation analysis")
# Required parameter
parser.add_argument('-i', type=str, required=True, help="The input sample information file")
parser.add_argument('-o','--output',type=str,required=True,help="Specify the output directory")
# Optional parameters
parser.add_argument('--chr',type=str,default='No',help="Bins divided between each order of magnitude")
parser.add_argument('--inter_chrom',type=str,default='False', help="The lower bound of the mean of the curve")
parser.add_argument('--width',type=float,default=8,help="The width of the plot")
parser.add_argument('--height',type=float,default=6,help="The height of the plot")

args = parser.parse_args()
meta_file = args.i
output_dir = args.output
chr = args.chr
inter_chrom = args.inter_chrom
width = args.width
height = args.height

def is_tool_available(name):
    return shutil.which(name) is not None

tools = ['fanc']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

def get_oe_matrix(hic,chr,output_dir,inter_chrom,prefix):
    # o/e
    if chr == 'No':
        command = f"fanc dump -e --only-intra {hic} {output_dir+'/'+prefix+'.oe.txt'}"
        if inter_chrom == 'True':
            command = f"fanc dump -e {hic} {output_dir+'/'+prefix+'.oe.txt'}"
    else:
        command = f"fanc dump -e --only-intra -s {chr}--{chr} {hic} {output_dir+'/'+prefix+'.oe.txt'}"
        if inter_chrom == 'True':
            command = f"fanc dump -e -s {chr}--{chr} {hic} {output_dir+'/'+prefix+'.oe.txt'}"
    try:
        subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")

meta = pd.read_csv(meta_file)
required_columns = ['file', 'tag','group']
assert all(col in meta.columns for col in required_columns), f"Require column {required_columns} for the input dataframe"
if meta['tag'].duplicated().any():
    print("Duplicates in tag column")
    sys.exit(1)

resolutions = []
for file in list(meta['file']):
    hic_check = fanc.load(file)
    resolutions.append(hic_check.bin_size)
    if len(set(resolutions)) != 1:
        print(f'Detect resolutions:{resolutions}')
        print('The resolutions of the hic files should be the same')
        sys.exit(1)

tmp_dir = os.path.join(output_dir,'tmp')
os.makedirs(output_dir,exist_ok=True)
os.makedirs(tmp_dir,exist_ok=True)

# Get matrix
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extract O/E matrix")
for hic, prefix in zip(meta['file'], meta['tag']):
    get_oe_matrix(hic,chr,tmp_dir,inter_chrom,prefix)
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} O/E matrix extraction done")

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Correlation analysis and plot")
script_dir = os.path.dirname(os.path.abspath(__file__))
cmd = ['Rscript',os.path.join(script_dir,'plot.R'),tmp_dir,meta_file,output_dir,str(width),str(height)]
subprocess.run(cmd)
print('Done, you can see the result now')

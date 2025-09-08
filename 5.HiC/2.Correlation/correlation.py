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
meta = args.i
output_dir = args.output
chr = args.chr
inter_chrom = args.inter_chrom
width = args.width
height = args.height

def is_tool_available(name):
    """检查命令是否存在于 PATH 中"""
    return shutil.which(name) is not None

tools = ['fanc']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

def get_oe_matrix(hic,chr,output_dir,inter_chrom):
    prefix = os.path.basename(hic).split('@')[0]
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

meta = pd.read_csv(meta)
required_columns = ['file', 'tag']
assert all(col in meta.columns for col in required_columns), f"Require column {required_columns} for the input dataframe"

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
for hic in list(meta['file']):
    get_oe_matrix(hic,chr,tmp_dir,inter_chrom)
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} O/E matrix extraction done")

# Combine the common matrix
oe_files = [os.path.join(tmp_dir,os.path.basename(i).split('@')[0] + '.oe.txt') for i in list(meta['file'])] 
tags = list(meta['tag'].astype(str))

dfs_to_merge = []
for tag, file in zip(tags,oe_files):
    temp_df = pd.read_csv(file,sep='\t',names=['chr1','start1','end1','chr2','start2','end2',tag])
    dfs_to_merge.append(temp_df)

combined_oe = reduce(
    lambda left, right: pd.merge(left, right, 
                               on=['chr1','start1','end1','chr2','start2','end2'],
                               how='inner'),
    dfs_to_merge
)

combined_oe = combined_oe.drop(['chr1','start1','end1','chr2','start2','end2'], axis=1)
combined_oe.columns = tags

# Correlation and PCA
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Start correlation analysis")
correlation_matrix = combined_oe.corr()

plt.figure(figsize=(width, height))
sns.heatmap(correlation_matrix, 
            annot=True, 
            cmap='Blues', 
            center=0, 
            square=True, 
            fmt='.3f',
            cbar_kws={'label': 'Pearson Correlation'})
plt.title('Sample Correlation Heatmap')
plt.tight_layout()
plt.savefig(os.path.join(output_dir,'correlation_heatmap.pdf'))

X = combined_oe.T
pca = PCA()
pca_result = pca.fit_transform(X)
pca_df = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(pca_result.shape[1])])
pca_df['Sample'] = combined_oe.columns
unique_samples = pca_df['Sample'].unique()
plt.figure(figsize=(width, height))
for sample in unique_samples:
    sample_data = pca_df[pca_df['Sample'] == sample]
    plt.scatter(sample_data['PC1'], sample_data['PC2'], 
                s=200, 
                alpha=0.8, 
                label=sample,
                edgecolors='black',
                linewidth=0.5)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)', fontsize=12)
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)', fontsize=12)
if chr == 'No':
    title = 'All chromosomes'
else:
    title = chr
plt.title(title, fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, 
           frameon=True, fancybox=True, shadow=True)
plt.tight_layout()
plt.savefig(os.path.join(output_dir,'correlation_pca.pdf'))

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Correlation analysis done")

print('Finished, you can see the results now')

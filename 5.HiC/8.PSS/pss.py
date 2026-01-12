# 8. PSS
import fanc
import subprocess
import pandas as pd
import sys
from datetime import datetime
import os
import shutil
import argparse

# Parameters
parser = argparse.ArgumentParser(description="Distance contact analysis")
# Required parameter
parser.add_argument('-i', type=str, required=True, help="The input sample information file")
parser.add_argument('-o', '--output', type=str, required=True, help="The output directory")
# Optional parameters
parser.add_argument('--bin',type=int,default=20,help="Bins divided between each order of magnitude")
parser.add_argument('--curve_mean_min',type=float,default=1e-6, help="The lower bound of the mean of the curve")
parser.add_argument('--curve_bin_min',type=float,default=0,help="The lower bound of bin index of the curve")
parser.add_argument('--curve_x_min',type=float,default=4,help="The minimal X-axis truncation value of the curve")
parser.add_argument('--curve_x_max',type=float,default=8,help="The maximal X-axis truncation value of the curve")
parser.add_argument('--expcurve_y_bias',type=float,default=-4,help="The Y-axis offset when plotting the difference from the expectation curve")
parser.add_argument('--expcurve_x_min',type=float,default=5,help="The minimal X-axis drop value when plotting the difference from the expectation curve")
parser.add_argument('--expcurve_x_max',type=float,default=8,help="The maximal X-axis drop value when plotting the difference from the expectation curve")
parser.add_argument('--heatmap_x_min',type=float,default=4.8,help="The minimal X-axis truncation value of the heatmap")
parser.add_argument('--heatmap_x_max',type=float,default=7.5,help="The maximal X-axis truncation value of the heatmap")
parser.add_argument('--width',type=float,default=10,help="The width of the plot")
parser.add_argument('--height',type=float,default=6,help="The height of the plot")

args = parser.parse_args()
meta = args.i
output_dir = args.output
bin_size=str(args.bin)
curve_mean_min=str(args.curve_mean_min)
curve_bin_min=str(args.curve_bin_min)
curve_x_min=str(args.curve_x_min)
curve_x_max=str(args.curve_x_max)
expcurve_y_bias=str(args.expcurve_y_bias)
expcurve_x_min=str(args.expcurve_x_min)
expcurve_x_max=str(args.expcurve_x_max)
heatmap_x_min=str(args.heatmap_x_min)
heatmap_x_max=str(args.heatmap_x_max)
width=str(args.width)
height=str(args.height)

current_dir = os.path.dirname(os.path.abspath(__file__))

def is_tool_available(name):
    return shutil.which(name) is not None

tools = ['fanc']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

def get_matrix(hic,output_dir):
    prefix = os.path.basename(hic).split('@')[0]
    command = f"fanc dump -u --only-intra {hic} {output_dir+'/'+prefix+'.contact.txt'}"
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

resolution = resolutions[0]
tmp_dir = os.path.join(output_dir,'tmp')
os.makedirs(output_dir,exist_ok=True)
os.makedirs(tmp_dir,exist_ok=True)

# Get matrix
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extract contact matrix")
for hic in list(meta['file']):
    get_matrix(hic,tmp_dir)
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Contact matrix extraction done")

# Compute distance_contact
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Compute distance_contact")
contact_files = [os.path.join(tmp_dir,os.path.basename(i).split('@')[0] + '.contact.txt') for i in list(meta['file'])] 
for file in contact_files:
    df = pd.read_csv(file,sep='\t',names=['chr','start1','end1','chr_end','start2','end2','contact'])
    df = df[df['chr'] == df['chr_end']]
    df['dist'] = abs(df['end2'] - df['end1'])
    df['contact'] = df['contact'].astype(int)
    df = df[['chr','dist','contact']]
    save_file = file.replace('.contact.txt','.dist.contact.txt')
    result = df.groupby(['chr','dist'])['contact'].sum()
    result.to_csv(save_file, sep='\t', index=True, encoding='utf-8')
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Distance_contact computation done")

dist_contact_files = [os.path.join(tmp_dir,os.path.basename(i).split('@')[0] + '.dist.contact.txt') for i in list(meta['file'])] 
file_list = ','.join(dist_contact_files)
tags = ','.join(list(meta['tag'].astype(str)))

# # PSS plot
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Plot PSS")
command = ['Rscript',os.path.join(current_dir,'pss.R'),file_list,tags,output_dir,str(resolution),
            bin_size,curve_mean_min,curve_bin_min,curve_x_min,curve_x_max,expcurve_y_bias,
            expcurve_x_min,expcurve_x_max,heatmap_x_min,heatmap_x_max,width,height]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
     print(f"Error：{e.stderr.decode()}")
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} PSS plot done")
print("Finished, you can check the results now.")

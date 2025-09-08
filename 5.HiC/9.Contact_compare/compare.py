# 9. Contact compare & delta plot
import fanc
import sys
from datetime import datetime
import argparse
import shutil
import os
import subprocess
import argparse

rscript = os.getenv("RSCRIPT_PATH", "Rscript")
current_dir = os.path.dirname(os.path.abspath(__file__))

# Check the environment
def is_tool_available(name):
    """检查命令是否存在于 PATH 中"""
    return shutil.which(name) is not None

def get_matrix(hic,chr,output_dir):
    prefix = os.path.basename(hic).split('@')[0]
    # count
    if chr == 'no':
        command = f"fanc dump -u --only-intra {hic} {output_dir+'/'+prefix+'.txt'}"
    else:
        command = f"fanc dump -u --only-intra -s {chr}--{chr} {hic} {output_dir+'/'+prefix+'.txt'}"
    try:
        subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")
    # o/e
    if chr == 'no':
        command = f"fanc dump -e --only-intra {hic} {output_dir+'/'+prefix+'.oe.txt'}"
    else:
        command = f"fanc dump -e --only-intra -s {chr}--{chr} {hic} {output_dir+'/'+prefix+'.oe.txt'}"
    try:
        subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error：{e.stderr.decode()}")

tools = ['fanc']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

# Parameters
parser = argparse.ArgumentParser(description="Distance contact analysis")
# Required parameter
parser.add_argument('--hic1', type=str, required=True, help="The first hic file")
parser.add_argument('--hic2', type=str, required=True, help="The second hic file")
parser.add_argument('-o','--output',type=str,required=True,help="Specify the output directory")
parser.add_argument('--chr',type=str,required=True,help="Specify the chromosome")
args = parser.parse_args()

hic1 = args.hic1
hic2 = args.hic2
output_dir = args.output
chr = args.chr

hic1_check = fanc.load(hic1)
hic2_check = fanc.load(hic2)

if hic1 == hic2:
    print('Please provide 2 different hic files')
    sys.exit(1)
if hic1_check.bin_size != hic2_check.bin_size:
    print('The resolutions of the 2 samples should be the same!')
    sys.exit(1)

tmp_dir = os.path.join(output_dir,'tmp')
os.makedirs(output_dir,exist_ok=True)
os.makedirs(tmp_dir,exist_ok=True)

# Get contact and O/E matrix
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Extract matrix")
get_matrix(hic1,chr,tmp_dir)
get_matrix(hic2,chr,tmp_dir)
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Matrix extraction done")

prefix1 = os.path.basename(hic1).split('@')[0]
prefix2 = os.path.basename(hic2).split('@')[0]

# Contact compare and Delta plot
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Do contact comparison and Delta plot")
command = [rscript,os.path.join(current_dir,'compare.R'),output_dir,prefix1,prefix2,chr]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Contact compare and Delta plot done")

print("Finished, you can check the results now.")

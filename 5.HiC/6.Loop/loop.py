import argparse
import os
import subprocess
import fanc
import sys
from datetime import datetime
import shutil

def is_tool_available(name):
    """检查命令是否存在于 PATH 中"""
    return shutil.which(name) is not None

tools = ['fanc']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

parser = argparse.ArgumentParser(description="FanC loop calling process")

# Required parameters
parser.add_argument('-i', '--input', type=str, required=True, help="Input Fanc .hic")
parser.add_argument('-o', '--output', type=str, required=True, help="Output directory")
# Optional parameters
parser.add_argument('-t','--thread', type=int, default=10, help="Number of threads, default: 10")

args = parser.parse_args()

output_directory = args.output
hic_file = args.input
thread = args.thread
tmp_directory = os.path.join(args.output,'tmp')

if not os.path.exists(output_directory):
    os.makedirs(output_directory)
if not os.path.exists(tmp_directory):
    os.makedirs(tmp_directory)

file_name = os.path.basename(hic_file)
last_dot_index = file_name.rfind('.')
sample = file_name[:last_dot_index]

hic_sample = fanc.load(hic_file)
chroms = hic_sample.chromosomes()
chroms = [chrom for chrom in chroms if '_' not in chrom]

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotate pixels for loop calling")
command = ["fanc","loops",hic_file,os.path.join(tmp_directory, sample+".loops"),"-t", str(thread), "-f","-c"]+chroms
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotation done")

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Filter annotated pixels")
command = ["fanc","loops",os.path.join(tmp_directory, sample+".loops"),os.path.join(tmp_directory, sample+"_filtered.loops"),"-t", str(thread), "-f",
           "--rh-filter","-d", "5","-o", "5"]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotation done")

check_loop = fanc.load(os.path.join(tmp_directory, sample+"_filtered.loops"))
if len(check_loop)==0:
    print('No loop after filteration')
    sys.exit(1)

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Merge unfiltered pixels into loops")
command = ["fanc","loops",os.path.join(tmp_directory, sample+"_filtered.loops"), os.path.join(tmp_directory, sample+"_merged.loops"),"-t", str(thread), "-f",
           "-j","--remove-singlets"]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Merge done")

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Export to BEDPE")
command = ["fanc","loops",os.path.join(tmp_directory, sample+"_merged.loops"),"-t", str(thread), "-f",
           "-b",os.path.join(output_directory, sample+".bedpe")]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Exportion done")

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loop aggregate plot")
command = ["fanc","aggregate",hic_file,os.path.join(output_directory, sample + ".bedpe"),os.path.join(tmp_directory, sample + ".agg"),
    "-p",os.path.join(output_directory, "loop_aggregate.pdf"),
    "--loops","--loop-strength",os.path.join(output_directory, "loop_strength.bed"),
    "--labels","Loop anchor","--label-locations","0.5"]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loop aggregate plot done")

input_file = os.path.join(output_directory, "loop_strength.bed")
output_file = os.path.join(output_directory, "arc.bed")

with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    for line in fin:
        cols = line.strip().split()
        selected = cols[:6] if len(cols) >= 6 else cols + ['.'] * (6 - len(cols))
        selected.append('1')
        fout.write('\t'.join(selected) + '\n')
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Done, you can check the results now")

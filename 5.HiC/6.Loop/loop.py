import argparse
import os
import subprocess
import fanc
import fanc.plotting as fancplot
import sys
from datetime import datetime
import shutil
import matplotlib.pyplot as plt
import pandas as pd

def is_tool_available(name):
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
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Filtration done")

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Export to BEDPE")
command = ["fanc","loops",os.path.join(tmp_directory, sample+"_merged.loops"),"-t", str(thread), "-f",
           "-b",os.path.join(tmp_directory, sample+".bedpe")]
try:
    subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"Error：{e.stderr.decode()}")
command = f"awk '$1 == $4 && $1 !~ /_/ && $4 !~ /_/' {os.path.join(tmp_directory, sample+'.bedpe')} > {os.path.join(output_directory, sample+'.loop.bedpe')}"
subprocess.run(command,check=True,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Exportion done")

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loop aggregate plot")
hic_rao = fanc.load(hic_file)
loops = fanc.load(os.path.join(output_directory, sample+".loop.bedpe"))
loops_am = fanc.AggregateMatrix.from_center_pairs(hic_rao, loops.region_pairs())
ax = fancplot.aggregate_plot(loops_am, vmin=-1, vmax=1,
                             relative_label_locations=[0.5],
                             labels=['loop anchor'])
fig = ax.get_figure()
fig.savefig(os.path.join(output_directory, "loop_aggregate.pdf"))
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loop aggregate plot done")

print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Grab O/E matrix and calculate loop strength")
command = ["fanc","dump","-e","--only-intra",hic_file,os.path.join(tmp_directory, sample + '.oe.txt')]
subprocess.run(command,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
loop_bed = pd.read_csv(os.path.join(output_directory, sample+".loop.bedpe"),sep='\t',names=['region1','start1','end1','region2','start2','end2','tag','loop'])
oe_bed =  pd.read_csv(os.path.join(tmp_directory,sample+".oe.txt"),sep='\t',names=['region1','start1','end1','region2','start2','end2','o/e'])
results = []
for chrom_pair in loop_bed[['region1', 'region2']].drop_duplicates().itertuples(index=False):
    chrom1, chrom2 = chrom_pair
    chrom_loops = loop_bed[
        (loop_bed['region1'] == chrom1) & 
        (loop_bed['region2'] == chrom2)
    ]
    chrom_oe = oe_bed[
        (oe_bed['region1'] == chrom1) & 
        (oe_bed['region2'] == chrom2)
    ]
    if len(chrom_oe) == 0:
        continue
    for _, loop_row in chrom_loops.iterrows():
        mask = (
            (chrom_oe['start1'] >= loop_row['start1']) &
            (chrom_oe['end1'] <= loop_row['end1']) &
            (chrom_oe['start2'] >= loop_row['start2']) &
            (chrom_oe['end2'] <= loop_row['end2'])
        )
        matching_oe = chrom_oe[mask]
        if len(matching_oe) > 0:
            avg_oe = matching_oe['o/e'].mean()
            results.append({
                'region1': chrom1,
                'start1': loop_row['start1'],
                'end1': loop_row['end1'],
                'region2': chrom2,
                'start2': loop_row['start2'],
                'end2': loop_row['end2'],
                'avg_o/e': avg_oe
            })
result_df = pd.DataFrame(results)
result_df.to_csv(os.path.join(output_directory, 'loop_strength.bed'),sep='\t')
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} O/E matrix extraction and loop strength calculation done")

input_file = os.path.join(output_directory, "loop_strength.bed")
output_file = os.path.join(output_directory, "arc.bed")

with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    for line in fin:
        cols = line.strip().split()
        selected = cols[:6] if len(cols) >= 6 else cols + ['.'] * (6 - len(cols))
        selected.append('1')
        fout.write('\t'.join(selected) + '\n')
print(f">>>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Done, you can check the results now")

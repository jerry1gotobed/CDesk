import sys
import subprocess
import os
import argparse
import shutil
import json

# Check tools
def is_tool_available(name):
    return shutil.which(name) is not None

tools = ['rmats2sashimiplot']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

config_path = os.environ.get('CDesk_config')

parser = argparse.ArgumentParser(description="CDesk bulkRNA rmats2sashimiplot process")
# Required parameters
parser.add_argument('--b1',type=str,required=True,help="The first txt file containing BAM files")
parser.add_argument('-s','--species',type=str,required=True,help="Specify the species")
parser.add_argument('-o','--output',type=str,required=True,help="The output directory")
parser.add_argument('--region',type=str,required=True,help="The genome region coordinates: format:{chromosome}:{strand}:{start}:{end}")
# Optional parameters
parser.add_argument('--l1',type=str,help="The first sample name",default = 'sample1')
parser.add_argument('--l2',type=str,help="The second sample name",default = 'sample2')
parser.add_argument('--b2',type=str,default='No',help="The second txt file containing BAM files")
parser.add_argument('--exon_s',type=int, default=1, required=True,help="How much to scale down exon, default: 1")
parser.add_argument('--intron_s',type=int, default=1, required=True,help="How much to scale down introns, default: 1")
parser.add_argument('--group',type=str,help="The path to a .gf file which groups the replicates",default='No')

args = parser.parse_args()

b1_file = args.b1
b2_file = args.b2
species = args.species
output = args.output
region = args.region
group = args.group
exon_s = str(args.exon_s)
intron_s = str(args.intron_s)
sample1 = args.l1
sample2 = args.l2

os.makedirs(output,exist_ok=True)

with open(config_path, "r") as f:
        config = json.load(f)
gff = config['data'][species]['gff3']

interval = region + ':' + gff

with open(b1_file, 'r', encoding='utf-8') as file:
    lines = file.readlines()
    lines = [line.strip() for line in lines]
    b1 = ','.join(lines)

command = ['rmats2sashimiplot','--b1',b1,'-c',interval,
           '--intron_s',intron_s,'--exon_s',exon_s,
           '-o',output,'--l1',sample1,'--l2',sample2]

if b2_file != 'No':
    with open(b2_file, 'r', encoding='utf-8') as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]
        b2 = ','.join(lines)
    command += ['--b2',b2]

if group != 'No':
    sample_names = []
    with open(group, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            sample_name = line.split(':', 1)[0].strip()
            sample_names.append(sample_name)

    sample1 = sample_names[0]
    sample2 = sample_names[1]

    command += ['--group-info',group]

#print(' '.join(command))
subprocess.run(command,check=True)

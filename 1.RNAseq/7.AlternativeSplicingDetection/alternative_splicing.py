import sys
import subprocess
import os
import json
import argparse
import shutil

# Check tools
def is_tool_available(name):
    return shutil.which(name) is not None

tools = ['rmats.py']
for tool in tools:
    if not is_tool_available(tool):
        print(f"❌ {tool} not available, please check the envionment variables")
        sys.exit(1)

config_path = os.environ.get('CDesk_config')
if os.environ.get('CDesk'):
    python = os.path.join(os.environ.get('CDesk'),'python')
else:
    python = 'python'

parser = argparse.ArgumentParser(description="CDesk bulkRNA rmats process")
# Required parameters
parser.add_argument('--b1',type=str,required=True,help="The first txt file containing BAM files")
parser.add_argument('-s','--species',type=str,required=True,help="Specify the species")
parser.add_argument('-o','--output',type=str,required=True,help="The output directory")
# Optional parameters
parser.add_argument('--b2',type=str,default='No',help="The second txt file containing BAM files")
parser.add_argument('--type',choices=['paired','single'],default='paired',required=True,help="Type of read: paired/single, default:paired")
parser.add_argument('-t','--thread', type=int, default=10, help="Number of threads, default: 10")
parser.add_argument('--length',required=int,default=150,help="The length of each read, default: 150")
parser.add_argument('--variable_read_length',choices=['True','False'],default='False',help="Allow reads with lengths that differ from read Length")
parser.add_argument('--allow_clipping',choices=['True','False'],default='False',help="Allow alignments with soft or hard clipping to be used")

args = parser.parse_args()
b1 = args.b1
b2 = args.b2
species = args.species
output = args.output
type = args.type
thread = args.thread
length = args.length
variable_read_length = args.variable_read_length
allow_clipping = args.allow_clipping

tmp = os.path.join(output,'tmp')
os.makedirs(output,exist_ok=True)
os.makedirs(tmp,exist_ok=True)

with open(config_path, "r") as f:
        config = json.load(f)
gtf = config['data'][species]['refseq_gtf']

command = [python,os.path.join(os.environ.get('CDesk'),'rmats.py'),'--b1',b1,'--gtf',gtf,'-t',type,'--readLength',str(length),'--nthread',str(thread),'--od',output,'--tmp',tmp]

with open(b1, 'r') as f:
    lines = f.read().strip().splitlines()
b1 = os.path.join(output,'b1.txt')
with open(b1, 'w') as f:
    f.write(','.join(lines))

if b2 != 'No':
    with open(b2, 'r') as f:
        lines = f.read().strip().splitlines()
    b2 = os.path.join(output,'b2.txt')
    with open(b2, 'w') as f:
        f.write(','.join(lines))
    command += ['--b2',b2]

if variable_read_length == 'True':
    command = command + ['--variable-read-length']

if allow_clipping == 'True':
    command = command + ['--allow-clipping']

subprocess.run(command,check=True)

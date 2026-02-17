import os 
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(description="Prepare csv input")
parser.add_argument("-i",required=True,type=str, help="Input directory")
parser.add_argument("-o",required=True,type=str, help="Output directory")
parser.add_argument("--pair",default='',type=str, help="Paire suffix")
parser.add_argument("--single",default='',type=str, help="Single suffix")

args = parser.parse_args()
inputdir = args.i
outroot = args.o
pair = args.pair
single = args.single

os.makedirs(outroot, exist_ok=True)

# Exception handling
if pair == '' and single == '':
    print('At least provide pair or single postfix')
    sys.exit(1)

pairs = pair.split(',')
if len(pairs) != 2:
    print('Pair postfix: comma separate, 2 postfixes!')
    sys.exit(1)

os.makedirs(outroot, exist_ok=True)

# Run
pair_samples = [];single_samples = []
fq_1 = []; fq_2 = [];ports = []
pair1_files = [];pair2_files=[];single_files=[]

if pair != '':
    temp_1 = [i.split(pairs[0])[0] for i in os.listdir(inputdir) if str.endswith(i,pairs[0])]
    temp_2 = [i.split(pairs[1])[0]  for i in os.listdir(inputdir) if str.endswith(i,pairs[1])]
    pair_samples = list(set(temp_1).intersection(set(temp_2)))
    pair1_files = [i + pairs[0] for i in pair_samples]
    pair2_files = [i + pairs[1] for i in pair_samples]

if single != '':
    single_samples = [i.split(single)[0]  for i in os.listdir(inputdir) if str.endswith(i,single)]
    single_files = [i + single for i in single_samples]

samples = single_samples + pair_samples

if len(samples) == 0:
    print('No sample found')
    sys.exit(1)

fq_1 = [os.path.join(inputdir,i) for i in single_files] + [os.path.join(inputdir,i) for i in pair1_files]
fq_2 = [os.path.join(inputdir,i) for i in single_files] + [os.path.join(inputdir,i) for i in pair2_files]
ports = [1] * len(single_samples) + [2] * len(pair_samples)

df = pd.DataFrame({
    'sample': samples,
    'fq1': fq_1,
    'fq2': fq_2,
    'ports': ports
})

df.to_csv(os.path.join(outroot,'Input.csv'),index=None)
print('Done, you can check the results now')

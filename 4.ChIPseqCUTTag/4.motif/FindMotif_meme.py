import os
import subprocess
import sys
import pandas as pd

input = sys.argv[1]
database = sys.argv[2]
output_dir = sys.argv[3]
thread = int(sys.argv[4]) #30
nmotifs = int(sys.argv[5]) # 3

os.makedirs(output_dir, exist_ok=True)

input = pd.read_csv(input)
if 'fasta' not in input.columns or 'sample' not in input.columns:
    print('Need columns: peak,sample')
    sys.exit(1)

fasta = list(input['fasta'])
sample = list(input['sample'])

with open(database, 'r') as file:
        databases = [line.strip() for line in file if line.strip()]

db = []
for i in databases:
    if not os.path.exists(os.path.abspath(i)):
        print(f'Can not find {i}')
        sys.exit(1)
    db = db + ['-db',i]

for fasta_file, sample_name in zip(fasta, sample):
    print(f'Processing {fasta_file} ......')
    cmd = [
    "meme-chip",
    "-meme-p", str(thread),
    "-oc",os.path.join(output_dir,sample_name),
    fasta_file,
    "-meme-nmotifs", str(nmotifs)
    ] + db
    #print(cmd) 
    with open(os.path.join(output_dir,sample_name+'.log'), 'w') as log:
        subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True
        )

print('Done, you can see the results now')

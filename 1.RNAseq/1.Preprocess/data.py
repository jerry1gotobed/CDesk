import os 
import pandas as pd
import argparse
import sys
import shutil

parser = argparse.ArgumentParser(description="Run CStreet analysis")
parser.add_argument("-i",required=True,type=str, help="Input expression matrixes list file")
parser.add_argument("-o",required=True,type=str, help="Input expression matrixes list file")
args = parser.parse_args()
input_file = args.i
outroot = args.o
os.makedirs(outroot, exist_ok=True)

data_dir = os.path.join(outroot,'data')
if os.path.exists(data_dir):
    if os.path.isdir(data_dir):
        shutil.rmtree(data_dir)
os.makedirs(data_dir, exist_ok=True)

bam_dir = os.path.join(outroot,'Bam')
os.makedirs(bam_dir, exist_ok=True)

paired_dir = data_dir
single_dir = data_dir

test = pd.read_csv(input_file)
# Check columns
for col in ("sample","fq1","fq2","ports"):
    if col not in test.columns:
        print(f"Error: Can not find column {col} in the input csv")
        sys.exit(1)

if test.duplicated(subset=['sample', 'ports']).any():
    print("Find duplications!")
    print(test[test.duplicated(subset=['sample', 'ports'], keep=False)])
    sys.exit(1)

if (~test['ports'].isin([1, 2])).any():
    invalid_values = test[~test['ports'].isin([1, 2])]['ports'].unique()
    print(f"Error ports column value: {invalid_values}. Only allow: 1,2")
    sys.exit(1)

check_cols = ["fq1", "fq2"]
if 'bam' in test.columns:
    check_cols = ["fq1", "fq2", "bam"]

for idx, row in test.iterrows():
    sample = row["sample"]
    # Check files exist or not
    for col in check_cols:
        raw = row[col]
        if pd.notna(raw) and raw.strip():
            src = raw
            if not os.path.exists(src):
                print(f"Can not find {sample} {col}: {src}")
                sys.exit(1)
    # Link files
    if row['ports'] == 2:
        if 'bam' in test.columns and pd.notna(row['bam']):
            if os.path.exists(os.path.join(bam_dir,sample+'.bam')):
                os.remove(os.path.join(bam_dir,sample+'.bam')) 
            os.symlink(os.path.abspath(row['bam']),os.path.join(bam_dir,sample+'.bam'))
        else:
            if pd.isna(row['fq1']) or pd.isna(row['fq2']):
                print(f"Must provide 2 fastq files for paired sequencing data {sample}")
                sys.exit(1)  
            if os.path.exists(os.path.join(paired_dir,sample+'_1.fq.gz')):
                os.remove(os.path.join(paired_dir,sample+'_1.fq.gz')) 
            os.symlink(os.path.abspath(row['fq1']),os.path.join(paired_dir,sample+'_1.fq.gz'))
            if os.path.exists(os.path.join(paired_dir,sample+'_2.fq.gz')):
                os.remove(os.path.join(paired_dir,sample+'_2.fq.gz')) 
            os.symlink(os.path.abspath(row['fq2']),os.path.join(paired_dir,sample+'_2.fq.gz'))

    if row['ports'] == 1:
        if 'bam' in test.columns and pd.notna(row['bam']):
            if os.path.exists(os.path.join(bam_dir,sample+'.bam')):
                os.remove(os.path.join(bam_dir,sample+'.bam')) 
            os.symlink(os.path.abspath(row['bam']),os.path.join(bam_dir,sample+'.bam'))
        else:
            if pd.isna(row['fq1']) and pd.isna(row['fq2']):
                print(f"At least provide one fastq file for single sequenceing data {sample}")
                sys.exit(1)

            if pd.notna(row['fq1']):
                if os.path.exists(os.path.join(single_dir,sample+'_single.fq.gz')):
                    os.remove(os.path.join(single_dir,sample+'_single.fq.gz')) 
                os.symlink(os.path.abspath(row['fq1']),os.path.join(single_dir,sample+'_single.fq.gz'))
            else:
                if os.path.exists(os.path.join(single_dir,sample+'_single.fq.gz')):
                    os.remove(os.path.join(single_dir,sample+'_single.fq.gz')) 
                os.symlink(os.path.abspath(row['fq2']),os.path.join(single_dir,sample+'_single.fq.gz'))  

if len(os.listdir(data_dir)) == 0:
    print('No fq files found')
    sys.exit(1)

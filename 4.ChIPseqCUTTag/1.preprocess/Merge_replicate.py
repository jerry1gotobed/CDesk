import os
import subprocess
import sys
import pandas as pd

def list_all_files_walk(folder_path):
    file_list = [];check=[]
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            file_list.append(file_path)
            check.append(file)
    return file_list,check

group = sys.argv[1]
bam_directory = sys.argv[2]

bams,check = list_all_files_walk(bam_directory)
group = pd.read_csv(group)

if 'group' not in group.columns:
    print('Skip merging bam files')
    sys.exit(0)

for i in group['sample']:
    if i+'.bam' not in check:
        print(f'{i}.bam does not exit')
        sys.exit(1)


# Single 
bam_dir = os.path.join(bam_directory,'single')
df = group[group['ports'] == 1]
if len(df) > 0:
    grouped_samples = df.groupby('group')['sample'].apply(list).to_dict()
    for group_name, samples in grouped_samples.items():
        if len(samples) > 1:  
            bam_files = [bams[check.index(sample+'.bam')] for sample in samples]  
            merged_bam_path = os.path.join(bam_dir, f"{group_name}.bam") 
            # Merge BAM files using samtools
            try:
                #print(f"Merging BAM files for group: {group_name}")
                subprocess.run(['samtools', "merge", merged_bam_path] + bam_files, check=True)
                #print(f"Successfully merged BAM files into {merged_bam_path}")
                for sample in samples:
                    subprocess.run(["rm", os.path.join(bam_dir, f"{sample}.bam")], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error during merging BAM files for group {group_name}: {e}")
        else:
            subprocess.run(["mv", os.path.join(bam_dir, f"{samples[0]}.bam"), os.path.join(bam_dir, f"{group_name}.bam")], check=True)

# Paired
bam_dir = os.path.join(bam_directory,'pair')
df = group[group['ports'] == 2]
if len(df) > 0:
    grouped_samples = df.groupby('group')['sample'].apply(list).to_dict()
    for group_name, samples in grouped_samples.items():
        if len(samples) > 1:  
            bam_files = [bams[check.index(sample+'.bam')] for sample in samples]  
            merged_bam_path = os.path.join(bam_dir, f"{group_name}.bam") 
            # Merge BAM files using samtools
            try:
                #print(f"Merging BAM files for group: {group_name}")
                subprocess.run(['samtools', "merge", merged_bam_path] + bam_files, check=True)
                #print(f"Successfully merged BAM files into {merged_bam_path}")
                for sample in samples:
                    subprocess.run(["rm", os.path.join(bam_dir, f"{sample}.bam")], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error during merging BAM files for group {group_name}: {e}")
        else:
            subprocess.run(["mv", os.path.join(bam_dir, f"{samples[0]}.bam"), os.path.join(bam_dir, f"{group_name}.bam")], check=True)

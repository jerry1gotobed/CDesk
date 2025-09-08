import os
import subprocess
import sys
import pandas as pd

group = sys.argv[1]
bam_dir = sys.argv[2]
samtools = sys.argv[3]

group = pd.read_csv(group)

for i in group['sample']:
    if not os.path.exists(os.path.join(bam_dir,f"{i}.bam")):
        print('{i}.bam不存在')
        sys.exit(1)

grouped_samples = group.groupby('group')['sample'].apply(list).to_dict()

for group_name, samples in grouped_samples.items():
    if len(samples) > 1:  
        bam_files = [os.path.join(bam_dir, f"{sample}.bam") for sample in samples]  
        merged_bam_path = os.path.join(bam_dir, f"{group_name}.bam") 
        # Merge BAM files using samtools
        try:
            print(f"Merging BAM files for group: {group_name}")
            subprocess.run([samtools, "merge", merged_bam_path] + bam_files, check=True)
            print(f"Successfully merged BAM files into {merged_bam_path}")
            for sample in samples:
                subprocess.run(["rm", os.path.join(bam_dir, f"{sample}.bam")], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error during merging BAM files for group {group_name}: {e}")
    else:
         subprocess.run(["mv", os.path.join(bam_dir, f"{samples[0]}.bam"), os.path.join(bam_dir, f"{group_name}.bam")], check=True)

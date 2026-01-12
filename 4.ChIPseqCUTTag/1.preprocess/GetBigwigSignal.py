import pyBigWig
import numpy as np
import pandas as pd
from collections import defaultdict
import sys

def calculate_gene_signal(bw_path, bed_path):
    """
    Caculate the average signal of each gene in the bigwig file

    Parameters:
    bw_path: bigWig file 
    bed_path: BED file
    """
    # Read from bed file
    with open(bed_path, 'r') as f:
        bed_lines = f.readlines()
    
    bed_data = []
    for line in bed_lines:
        # Skip empty and annotation line
        if line.startswith('#') or line.strip() == '':
            continue
            
        parts = line.strip().split()
        if len(parts) < 5:
            print('promoter.bed needs to have at least 5 columns: chrom start end transcript gene!')
            sys.exit(1)
        chrom = parts[0]
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError:
            continue  # Skip invalid coordinate rows
        transcript = parts[3]
        gene = parts[4]
        #strand = parts[5]
        bed_data.append((chrom, start, end, transcript, gene))
    
    # Store all transcript intervals grouped by gene
    transcript_intervals = defaultdict(list)
    for chrom, start, end, transcript, gene in bed_data:
        transcript_intervals[transcript].append((chrom, start, end,gene))
    
    # Calculate each gene signal
    results = [];wrong = 0
    try:
        bw = pyBigWig.open(bw_path)
        bw_chroms = set(bw.chroms().keys())  # Get the chrom list of the bigWig
        
        for transcript, intervals in transcript_intervals.items():
            total_signal = 0
            total_length = 0
            
            for chrom, start, end,gene in intervals:
                # Check whether the chrom is in bigWig file
                if chrom not in bw_chroms:
                    continue
                
                # Get the interval value
                try:
                    values = np.array(bw.values(chrom, start, end))
                    # Remove NaN
                    valid_values = values[~np.isnan(values)]
                    
                    if len(valid_values) > 0:
                        # Calculate the weighted signal value for this transcript interval
                        interval_signal = np.sum(valid_values)
                        interval_length = len(valid_values)
                        
                        total_signal += interval_signal
                        total_length += interval_length
                except Exception as e:
                    wrong += 1
            
            # Calculate the average signal of the gene
            if total_length > 0:
                mean_signal = total_signal / total_length
            else:
                mean_signal = np.nan
            
            results.append({
                "gene": gene,
                "transcript": transcript,
                "mean_signal": mean_signal
                #"num_transcripts": len(intervals)
                #"covered_bases": total_length
            })
    finally:
        if 'bw' in locals() and bw is not None:
            bw.close()
    
    # Transfer to DataFrame and sort
    df = pd.DataFrame(results)
    df = df.sort_values(by="transcript")
    print(f'Skip {wrong} transcripts')
    return df

if __name__ == "__main__":
    bed_file = sys.argv[1]
    bw_file = sys.argv[2]
    out_file = sys.argv[3]

    # Calculate gene signal
    gene_signal_df = calculate_gene_signal(bw_file, bed_file)
     
    # Save
    gene_signal_df.to_csv(out_file, index=False, sep=",")
    print(f"Results saved to {out_file}")

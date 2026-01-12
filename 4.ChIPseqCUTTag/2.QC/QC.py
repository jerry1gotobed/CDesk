import argparse
import subprocess
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
import pyBigWig
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run ATAC sample replicate correlation")
    parser.add_argument("-i", required=True,type=str, help="The input ATAC sample information file")
    parser.add_argument("--species", required=True,type=str, help="Species")
    parser.add_argument("--bin",required=True, type=str, help="Genome partition bin")
    parser.add_argument("--step",required=True, type=str, help="BW file partition step")
    parser.add_argument("-o",required=True, type=str, help="Output directory")
    parser.add_argument("-t",required=True, type=str, help="Threads to run")
    parser.add_argument('--width',default=12,type=int, help='The plot width, default: 12')
    parser.add_argument('--height',default=8,type=int, help='The plot height, default: 8')
    return parser.parse_args()

def Sample_Correlation_Genome(meta, chrom_size, bin_size, output_dir,width,height,dataframe):
    # Check files
    if not os.path.isfile(chrom_size):
        print(f"Error：chromosome size file {chrom_size} not valid")
        sys.exit(1)

    # Grab all bw files
    bw_files = list(meta['bw'])

    # Read chrom size
    with open(chrom_size, "r") as f:
        len_info = f.readlines()

    # Process each bw file
    for file in bw_files:
        prefix = os.path.basename(file).split('.')[0]
        with open(os.path.join(output_dir,'tmp',prefix + '.SignalAcrossGenome.txt'), "w") as output_info:
            for each in len_info:
                each = each.strip().split()
                tmp_chr = each[0]
                end = int(each[1])
                steps = int(end / bin_size)
                if "_" not in tmp_chr:
                    cmd = ['bigWigSummary', file, tmp_chr, "1", str(end), str(steps)]
                    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                    result, err = pipe.communicate()
                    if not result:
                        output_info.write("\n".join(["NA"] * steps) + "\n")
                    else:
                        if isinstance(result, bytes):
                            try:
                                result = result.decode("utf-8")
                            except UnicodeDecodeError:
                                print('Error: The output cannot be decoded. Please check the result format of bigWigSummary')
                                exit()
                        output_info.write(result.replace("\t", "\n").replace("n/a", "NA"))
    
    # with open(os.path.join(output_dir,'tmp','SignalAcrossGenome_index.txt'), "w") as output_info:
    #     for each in len_info:
    #         each = each.strip().split()
    #         tmp_chr = each[0]
    #         end = int(each[1])
    #         steps = int(end / bin_size)
    #         if "_" not in tmp_chr:
    #             output_info.write(tmp_chr + "\t" + str(steps) + "\n")

    # Plot
    cmd = ['Rscript',os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SampleCorrelation.Genome.R'), output_dir, dataframe,str(width),str(height)]
    subprocess.run(cmd, check=True)


def process_peaks(bw_file, bed_file, step, out_file):
    with pyBigWig.open(bw_file) as bw, open(bed_file, "r") as bed_info, open(out_file, "w") as out_info:
        for each in bed_info:
            each = each.strip().split()
            chr = each[0]
            start = int(each[1])
            end = int(each[2])
            signal = bw.stats(chr, start, end, nBins=step)
            if signal is None:
                out_info.write("\t".join(["NA"] * step) + "\n")
            else:
                signal = [str(x) if x is not None else "NA" for x in signal]
                out_info.write("\t".join(signal) + "\n")

def Sample_Correlation_Peak(meta,output_dir, step,threads,width,height,dataframe):
    # Merge peaks
    peak_files = list(meta['peak'])
    peak_files_str = " ".join(peak_files)
    output_file = os.path.join(output_dir,'tmp',"merged_peaks.bed")
    command = f"cat {peak_files_str} | cut -f 1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > {output_file}"
    subprocess.run(command, shell=True, check=True)

    # Grab all bw files
    bw_files = list(meta['bw'])
    tasks = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for file in bw_files:
            prefix = os.path.basename(file).split('.')[0]
            out_file = os.path.join(output_dir,'tmp',prefix + '.SignalOnPeaks.txt')
            tasks.append(executor.submit(process_peaks, file, output_file, step, out_file))
        for future in tasks:
            future.result()
    # plot
    r_cmd = ['Rscript',os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SampleCorrelation.Peaks.R'), output_dir, dataframe,str(width),str(height)]
    subprocess.run(r_cmd, check=True)

def main():
    args = parse_arguments()
    
    output_dir = args.o
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir,'tmp'),exist_ok=True)

    meta = args.i
    meta = pd.read_csv(meta)

    required = ['tag', 'group','bw']
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"Miss required columns: {missing}")
   
    if meta['tag'].duplicated().any():
        print("Duplicates in tag column")
        sys.exit(1)

    bin_size = int(args.bin) 
    config_path = os.environ.get('CDesk_config')
    with open(config_path, "r") as f:
        config = json.load(f)
    chrom_size = config['data'][args.species]['chromInfo']
    
    step = int(args.step)
    thread = int(args.t)

    width = args.width
    height = args.height
    # Sample Correlation in Genome
    for file in meta['bw']:
        if not os.path.exists(file):
            print(f'Can not find {file}')
            sys.exit(1)

    if 'peak' in meta.columns:
        for file in meta['peak']:
            if not os.path.exists(file):
                print(f'Can not find {file}')
                sys.exit(1)
    print('Process on genome scale')
    Sample_Correlation_Genome(meta,chrom_size,bin_size,output_dir,width,height,args.i)
    if 'peak' in meta.columns:
        print('Process on merged peaks')
        Sample_Correlation_Peak(meta,output_dir,step,thread,width,height,args.i)

if __name__ == "__main__":
    main()
    print('Done, you can check the results now.')

import argparse
import subprocess
import json
import os
import math
import sys
import glob
import shutil
from concurrent.futures import ThreadPoolExecutor
import pyBigWig

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run ATAC sample replicate correlation")
    parser.add_argument("--bw",required=True, type=str, help="Input bw file directory")
    parser.add_argument("--species", required=True,type=str, help="Species")
    parser.add_argument("--bin",required=True, type=str, help="Genome partition bin")
    parser.add_argument("--group", required=True,type=str, help="Grouping file")

    parser.add_argument("--peak",required=True, type=str, help="Input peak file directory")
    parser.add_argument("--suffix",type=str,default='_peaks.bed', help="Peak file suffix, default: _peaks.bed")
    parser.add_argument("--step",required=True, type=str, help="BW file partition step")
    parser.add_argument("-o",required=True, type=str, help="Output directory")
    parser.add_argument("-t",required=True, type=str, help="Threads to run")
    parser.add_argument('--width',default=12,type=int, help='The plot width, default: 12')
    parser.add_argument('--height',default=8,type=int, help='The plot height, default: 8')
    parser.add_argument('--batch', default='no',choices=['no','removeBatchEffect','combat'],help='Choose the method to remove batch effect by the batch column in the grouping file') 
    parser.add_argument("--config_path", type=str,help="Path to configuration file")
    return parser.parse_args()

def Sample_Correlation_Genome(bw, chrom_size, bin_size, output_dir, meta,width,height,config):
    # 检查输入文件夹和文件
    if not os.path.isdir(bw):
        print(f"错误：{bw} 不是有效的文件夹路径")
        exit()
    if not os.path.isfile(chrom_size):
        print(f"错误：{chrom_size} 不是有效的文件路径")
        exit()

    # 获取所有 .bw 文件
    bw_files = [f for f in os.listdir(bw) if f.endswith(".bw")]
    if len(bw_files) == 0:
        print("输入文件夹无 .bw 文件")
        exit()

    # 将染色体信息文件读取到内存
    with open(chrom_size, "r") as f:
        len_info = f.readlines()

    # 处理每个 .bw 文件
    for file in bw_files:
        input_bw = os.path.join(bw, file)
        prefix = file.replace('.bw', '')
        with open(os.path.join(output_dir, prefix + '.SignalAcrossGenome.txt'), "w") as output_info:
            for each in len_info:
                each = each.strip().split()
                tmp_chr = each[0]
                end = int(each[1])
                steps = int(end / bin_size)
                if "_" not in tmp_chr:
                    cmd = [config['software']['bigWigSummary'], input_bw, tmp_chr, "1", str(end), str(steps)]
                    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                    result, err = pipe.communicate()
                    if not result:
                        output_info.write("\n".join(["NA"] * steps) + "\n")
                    else:
                        if isinstance(result, bytes):
                            try:
                                result = result.decode("utf-8")
                            except UnicodeDecodeError:
                                print("错误：无法解码输出，请检查 bigWigSummary 的结果格式")
                                exit()
                        output_info.write(result.replace("\t", "\n").replace("n/a", "NA"))
    
    with open(os.path.join(output_dir, 'SignalAcrossGenome_index.txt'), "w") as output_info:
        for each in len_info:
            each = each.strip().split()
            tmp_chr = each[0]
            end = int(each[1])
            steps = int(end / bin_size)
            if "_" not in tmp_chr:
                output_info.write(tmp_chr + "\t" + str(steps) + "\n")

    # 调用 R 脚本
    cmd = ['Rscript',os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SampleCorrelation.Genome.R'), output_dir, meta,str(width),str(height)]
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

def Sample_Correlation_Peak(peak, bw, output_dir, step, meta,threads,suffix,width,height,config):
    # 检查输入文件夹

    if not os.path.isdir(peak):
        print(f"错误：{peak} 不是有效的文件夹路径")
        exit()
    if not os.path.isdir(bw):
        print(f"错误：{bw} 不是有效的文件夹路径")
        exit()

    # 合并峰文件

    peak_files = glob.glob(os.path.join(peak, "*"+suffix))
    if not peak_files:
        print(f"输入文件夹无{suffix}文件.")
        exit()
    peak_files_str = " ".join(peak_files)
    output_file = os.path.join(output_dir, "merged_peaks.bed")
    command = f"cat {peak_files_str} | cut -f 1-3 | sort -k1,1 -k2,2n | {config['software']['bedtools']} merge -i - > {output_file}"
    subprocess.run(command, shell=True, check=True)

    # 获取所有 .bw 文件

    bw_files = [f for f in os.listdir(bw) if f.endswith(".bw")]
    tasks = []
    with ThreadPoolExecutor(max_workers=threads) as executor:  # 使用多线程同时处理

        for file in bw_files:
            prefix = file.replace('.bw', '')
            bw_file = os.path.join(bw, file)
            out_file = os.path.join(output_dir, prefix + '.SignalOnPeaks.txt')
            tasks.append(executor.submit(process_peaks, bw_file, output_file, step, out_file))
        for future in tasks:
            future.result()  # 等待所有任务完成
    # 调用 R 脚本
    r_cmd = ['Rscript',os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SampleCorrelation.Peaks.R'), output_dir, meta,str(width),str(height)]
    subprocess.run(r_cmd, check=True)

def main():
    # 解析命令行参数
    args = parse_arguments()
    
    output_dir = args.o
    os.makedirs(output_dir, exist_ok=True)

    bw = args.bw
    bin_size = int(args.bin) #500
    with open(args.config_path, "r") as f:
        config = json.load(f)
    chrom_size = config['data'][args.species]['chromInfo']
    meta = args.group

    peak = args.peak
    step = int(args.step)
    thread = int(args.t)
    suffix = args.suffix
    width = args.width
    height = args.height
    batch = args.batch
    # Sample Correlation in Genome
    Sample_Correlation_Genome(bw,chrom_size,bin_size,output_dir,meta,width,height,config)
    if args.peak != '':
        Sample_Correlation_Peak(peak,bw,output_dir,step,meta,thread,suffix,width,height,config)

if __name__ == "__main__":
    main()
    print('运行完成')

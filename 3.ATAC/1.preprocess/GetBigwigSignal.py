import pyBigWig
import numpy as np
import pandas as pd
from collections import defaultdict
import sys

def calculate_gene_signal(bw_path, bed_path):
    """
    计算每个基因在bigWig文件上的平均信号值
    
    参数:
    bw_path: bigWig文件路径
    bed_path: BED文件路径
    """
    # 从BED文件中读取内容
    with open(bed_path, 'r') as f:
        bed_lines = f.readlines()
    
    bed_data = []
    for line in bed_lines:
        # 跳过注释行和空行
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
            continue  # 跳过无效坐标行
        transcript = parts[3]
        gene = parts[4]
        #strand = parts[5]
        bed_data.append((chrom, start, end, transcript, gene))
    
    # 按基因分组，存储所有转录本区间
    transcript_intervals = defaultdict(list)
    for chrom, start, end, transcript, gene in bed_data:
        transcript_intervals[transcript].append((chrom, start, end,gene))
    
    # 计算每个基因的信号值
    results = [];wrong = 0
    try:
        bw = pyBigWig.open(bw_path)
        bw_chroms = set(bw.chroms().keys())  # 获取bigWig中存在的染色体列表
        
        for transcript, intervals in transcript_intervals.items():
            total_signal = 0
            total_length = 0
            
            for chrom, start, end,gene in intervals:
                # 检查染色体是否在bigWig文件中
                if chrom not in bw_chroms:
                    continue
                
                # 获取区间信号值
                try:
                    values = np.array(bw.values(chrom, start, end))
                    # 移除NaN值
                    valid_values = values[~np.isnan(values)]
                    
                    if len(valid_values) > 0:
                        # 计算该转录本区间的加权信号值
                        interval_signal = np.sum(valid_values)
                        interval_length = len(valid_values)
                        
                        total_signal += interval_signal
                        total_length += interval_length
                except Exception as e:
                    wrong += 1
            
            # 计算该基因的平均信号值
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
    
    # 转换为DataFrame并排序
    df = pd.DataFrame(results)
    df = df.sort_values(by="transcript")
    print(f'Skip {wrong} transcripts')
    return df

# 示例用法
if __name__ == "__main__":
    bed_file = sys.argv[1]
    bw_file = sys.argv[2]
    out_file = sys.argv[3]

    # 计算基因信号值
    gene_signal_df = calculate_gene_signal(bw_file, bed_file)
     
    # 保存到CSV文件
    gene_signal_df.to_csv(out_file, index=False, sep=",")
    print(f"Results saved to {out_file}")

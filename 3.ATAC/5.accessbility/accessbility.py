#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
import pyBigWig
import tempfile
import subprocess
import re
import json

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="分析基因组区域可访问性并绘制profile和boxplot")
    
    # 必需参数
    #parser.add_argument("-b", "--bed", required=True, help="输入bed文件列表，用逗号分隔")
    parser.add_argument("-bw", "--bigwig", required=True, help="bigwig文件列表文件")
    parser.add_argument("-o", "--output", required=True, help="输出目录")
    parser.add_argument("-m", "--mode", required=True, choices=["region", "gene"], help="分析模式：region或gene")
    parser.add_argument("-i", "--input", required=True, help="输入区域列表文件")
    
    # 可选参数
    parser.add_argument("-s", "--species", help="物种信息")
    parser.add_argument("--method", choices=["sample", "region_list"], default="sample", 
                        help="boxplot绘制方法：sample（样本为列）或region_list（区域列表为列）")
    parser.add_argument("--region_list", help="region_list模式下的区域列表文件（每行包含多个基因/区域）")
    parser.add_argument("--center", help="average_profile对齐方式：region模式可选peak/position，gene模式可选TSS/TTS/whole/peak")
    parser.add_argument("--promoter_size", type=int, default=5000, help="基因启动子区域大小（默认：TSS上下游各5kb）")
    #parser.add_argument("--threads","-t", type=int, default=4, help="使用的CPU线程数（默认：4）")
    parser.add_argument("config_path", help="配置信息文件")
    
    args = parser.parse_args()
    
    # 参数验证
    if args.mode == "gene" and not args.species:
        parser.error("gene模式下需要提供物种信息")
    
    if args.method == "region_list" and not args.region_list:
        parser.error("region_list模式下需要提供--region_list参数")
    
    for bw in args.bigwig.split(","):
        if not os.path.exists(bw):
            parser.error(f"BigWig文件不存在: {bw}")
    
    if not os.path.exists(args.input):
        parser.error(f"输入区域列表文件不存在: {args.input}")
    
    if args.region_list and not os.path.exists(args.region_list):
        parser.error(f"区域列表文件不存在: {args.region_list}")
    
    # 设置center参数的默认值和验证
    if args.center is None:
        args.center = "position" if args.mode == "region" else "whole"
    else:
        if args.mode == "region" and args.center not in ["peak", "position"]:
            parser.error("region模式下--center参数只能是peak或position")
        elif args.mode == "gene" and args.center not in ["TSS", "TTS", "whole", "peak"]:
            parser.error("gene模式下--center参数只能是TSS、TTS、whole或peak")
    
    # 确保输出目录存在
    os.makedirs(args.output, exist_ok=True)
    
    return args

# 基因模式下先读取基因，然后去获得promoter区域
def read_gene_list(gene_list_file):
    """读取基因列表文件"""
    try:
        with open(gene_list_file, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        print(f"成功读取 {len(genes)} 个基因")
        return genes
    except Exception as e:
        print(f"读取基因列表文件时出错: {e}", file=sys.stderr)
        sys.exit(1)

# 区域模式下读取区域
def read_region_list(region_file):
    """读取区域列表文件"""
    try:
        regions_df = pd.read_csv(region_file, sep='\t', header=None)
        if regions_df.shape[1] < 3:
            raise ValueError("区域文件必须至少包含3列：chr, start, end")
        
        # 仅使用前三列（chr, start, end）
        regions_df = regions_df.iloc[:, 0:3]
        regions_df.columns = ['chr', 'start', 'end']
        
        print(f"成功读取 {len(regions_df)} 个区域")
        return regions_df
    except Exception as e:
        print(f"读取区域列表文件时出错: {e}", file=sys.stderr)
        sys.exit(1)

# 读取region_list文件，每行包含多个基因/区域
def read_region_list_file(region_list_file):
    """读取region_list文件，每行包含多个基因/区域"""
    try:
        region_groups = []
        with open(region_list_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line:
                    # 分割每行的基因/区域
                    regions = [r.strip() for r in line.split() if r.strip()]
                    if regions:
                        region_groups.append({
                            'group_id': f"Group_{line_num}",
                            'regions': regions
                        })
        
        print(f"成功读取 {len(region_groups)} 个区域组")
        return region_groups
    except Exception as e:
        print(f"读取region_list文件时出错: {e}", file=sys.stderr)
        sys.exit(1)

# 基因模式下获得promoter区域
def get_gene_promoter_regions(genes, gtf_file, promoter_size=5000):
    """获取基因启动子区域（TSS上下游各5kb）从GTF文件中"""
    # 使用临时文件存储从基因组文件中提取的基因坐标
    with tempfile.NamedTemporaryFile(mode='w+t', suffix='.bed', delete=False) as temp_file:
        temp_genes_bed = temp_file.name
    
    try:
        # 从GTF文件中获取基因信息
        print(f"从GTF文件中提取基因信息...")
        
        gene_regions = []
        for gene in genes:
            try:
                # 使用grep搜索GTF文件中的基因信息，筛选transcript行
                cmd = f"grep -w 'transcript' {gtf_file} | grep 'gene_name \"{gene}\"'"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                
                if result.returncode != 0 or not result.stdout:
                    print(f"警告: 未找到基因 {gene} 的信息", file=sys.stderr)
                    continue
                
                # 解析基因信息（GTF格式）
                # GTF格式：chr source feature start end score strand frame attributes
                transcript_lines = result.stdout.strip().split('\n')
                
                # 解析第一个匹配的转录本
                # 如果有多个转录本，可以选择最长的或者指定规则的一个
                line = transcript_lines[0].strip()
                fields = line.split('\t')
                
                if len(fields) >= 9:
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    
                    # 确定TSS和TTS位置（根据基因方向）
                    if strand == '+':
                        tss = start  # 正链：TSS是起始位置
                        tts = end    # 正链：TTS是结束位置
                    else:
                        tss = end    # 负链：TSS是结束位置
                        tts = start  # 负链：TTS是起始位置
                    
                    # 计算启动子区域（TSS上下游各promoter_size）
                    promoter_start = max(0, tss - promoter_size)
                    promoter_end = tss + promoter_size
                    
                    gene_regions.append({
                        'chr': chrom,
                        'start': promoter_start,
                        'end': promoter_end,
                        'name': gene,
                        'score': 0,
                        'strand': strand,
                        'tss': tss,
                        'tts': tts,
                        'gene_start': start,  # 保存原始基因起始位置
                        'gene_end': end       # 保存原始基因结束位置
                    })
                    
                    print(f"找到基因 {gene}: {chrom}:{promoter_start}-{promoter_end} ({strand})")
            except Exception as e:
                print(f"处理基因 {gene} 时出错: {e}", file=sys.stderr)
        
        # 将基因区域写入临时BED文件
        if gene_regions:
            gene_df = pd.DataFrame(gene_regions)
            gene_df.to_csv(temp_genes_bed, sep='\t', header=False, index=False)
            print(f"成功生成 {len(gene_regions)} 个基因的启动子区域")
        else:
            print("错误: 未找到有效的基因区域", file=sys.stderr)
            sys.exit(1)
        
        return gene_df
        
    except Exception as e:
        print(f"获取基因启动子区域时出错: {e}", file=sys.stderr)
        sys.exit(1)

# 计算区域矩阵并生成平均profile图
def compute_matrix_and_profile(regions_df, bigwig_files, output_dir, mode, center):#, threads=4):
    """计算区域矩阵并生成平均profile图"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 结果存储
    profile_results = {}  # 存储每个区域的profile数据
    region_ids = []       # 区域ID列表
    
    # 对每个区域处理
    for idx, row in regions_df.iterrows():
        chrom, start, end = row['chr'], row['start'], row['end']
        region_id = f"{chrom}:{start}-{end}"
        region_ids.append(region_id)
        print(f"处理区域: {region_id}")
        
        # 确定TSS和TTS位置
        if 'strand' in row and 'name' in row:  # gene模式
            tss = row['tss'] if 'tss' in row else (start if row['strand'] == '+' else end)
            tts = row['tts'] if 'tts' in row else (end if row['strand'] == '+' else start)
            gene_name = row['name']
        else:  # region模式
            tss = (start + end) // 2
            tts = (start + end) // 2  # region模式下TSS和TTS相同
            gene_name = region_id
        
        # 根据center参数定义区域边界
        if center == "peak":
            # peak模式：需要先获取完整的gene/region信号来找peak summit
            if 'gene_start' in row and 'gene_end' in row:
                # gene模式：使用完整基因区域
                full_start = row['gene_start']
                full_end = row['gene_end']
            else:
                # region模式：使用给定的区域
                full_start = start
                full_end = end
            region_start = max(0, full_start)
            region_end = full_end
            center_pos = None  # 暂时设为None，后面会根据peak summit确定
        elif center == "TTS" and 'tts' in row:
            # TTS模式：TTS上下游各5000bp
            center_pos = tts
            region_start = max(0, center_pos - 5000)
            region_end = center_pos + 5000
        else:
            # 默认TSS模式：TSS上下游各5000bp
            center_pos = tss
            region_start = max(0, center_pos - 5000)
            region_end = center_pos + 5000
        
        # 存储每个样本在该区域的信号
        sample_signals = []
        
        # 对每个样本获取信号
        for bw_file in bigwig_files:
            try:
                bw = pyBigWig.open(bw_file)
                
                # 确保染色体存在于bigwig文件中
                if chrom in bw.chroms():
                    # 获取区域内每个位置的信号值
                    try:
                        # 调整区域边界确保在染色体范围内
                        valid_start = max(0, region_start)
                        valid_end = min(bw.chroms()[chrom], region_end)
                        
                        # 获取信号值
                        values = bw.values(chrom, valid_start, valid_end)
                        
                        # 处理可能的NaN值
                        values = [0 if v is None or np.isnan(v) else v for v in values]
                        
                        # 如果区域开始位置被调整，填充0
                        if valid_start > region_start:
                            pad_left = valid_start - region_start
                            values = [0] * pad_left + values
                        
                        # 如果区域结束位置被调整，填充0
                        if valid_end < region_end:
                            pad_right = region_end - valid_end
                            values = values + [0] * pad_right
                        
                        # 确保长度正确
                        if len(values) != (region_end - region_start):
                            # 调整为正确长度
                            values = values[:region_end - region_start]
                            if len(values) < (region_end - region_start):
                                values = values + [0] * ((region_end - region_start) - len(values))
                        
                    except Exception as e:
                        print(f"  警告: 获取区域 {region_id} 的信号时出错: {e}")
                        values = [0] * (region_end - region_start)
                else:
                    print(f"  警告: 染色体 {chrom} 不存在于文件 {bw_file} 中")
                    values = [0] * (region_end - region_start)
                
                sample_signals.append(values)
                bw.close()
                
            except Exception as e:
                print(f"  警告: 处理文件 {bw_file} 时出错: {e}")
                # 出错时填充0
                sample_signals.append([0] * (region_end - region_start))
        
        # 计算所有样本的平均信号
        if sample_signals:
            # 确保所有样本信号长度一致
            min_length = min(len(s) for s in sample_signals)
            sample_signals = [s[:min_length] for s in sample_signals]
            
            # 计算平均值
            avg_signal = np.mean(sample_signals, axis=0)
            
            # peak模式特殊处理：找到peak summit并重新提取信号
            if center == "peak":
                # 找到平均信号中的peak summit位置
                peak_idx = np.argmax(avg_signal)
                peak_summit_pos = region_start + peak_idx
                
                print(f"  找到peak summit: {chrom}:{peak_summit_pos}")
                
                # 重新提取peak summit周围5kb的信号
                summit_start = max(0, peak_summit_pos - 5000)
                summit_end = peak_summit_pos + 5000
                
                # 重新获取所有样本在summit周围的信号
                summit_signals = []
                for bw_file in bigwig_files:
                    try:
                        bw = pyBigWig.open(bw_file)
                        if chrom in bw.chroms():
                            try:
                                valid_start = max(0, summit_start)
                                valid_end = min(bw.chroms()[chrom], summit_end)
                                values = bw.values(chrom, valid_start, valid_end)
                                values = [0 if v is None or np.isnan(v) else v for v in values]
                                
                                # 填充到正确长度
                                if valid_start > summit_start:
                                    pad_left = valid_start - summit_start
                                    values = [0] * pad_left + values
                                if valid_end < summit_end:
                                    pad_right = summit_end - valid_end
                                    values = values + [0] * pad_right
                                
                                # 确保长度正确
                                target_length = summit_end - summit_start
                                if len(values) != target_length:
                                    values = values[:target_length]
                                    if len(values) < target_length:
                                        values = values + [0] * (target_length - len(values))
                                        
                            except Exception as e:
                                print(f"    警告: 获取summit区域信号时出错: {e}")
                                values = [0] * (summit_end - summit_start)
                        else:
                            values = [0] * (summit_end - summit_start)
                        summit_signals.append(values)
                        bw.close()
                    except Exception as e:
                        print(f"    警告: 处理summit信号时出错: {e}")
                        summit_signals.append([0] * (summit_end - summit_start))
                
                # 重新计算summit周围的平均信号
                if summit_signals:
                    min_length = min(len(s) for s in summit_signals)
                    summit_signals = [s[:min_length] for s in summit_signals]
                    avg_signal = np.mean(summit_signals, axis=0)
                    # 更新区域边界为summit周围的实际范围
                    region_start = summit_start
                    region_end = summit_start + min_length
                    center_pos = peak_summit_pos
                    print(f"    重新提取信号: {region_start}-{region_end}, summit: {peak_summit_pos}")
                else:
                    print(f"    警告: 无法获取summit周围信号，使用原始信号")
                    center_pos = peak_summit_pos
            
            # 重新采样到500个bin以获得更平滑的profile并保持peak位置精度
            n_bins = 500
            original_length = region_end - region_start
            if len(avg_signal) > n_bins:
                # 计算每个bin的大小
                bin_size = len(avg_signal) // n_bins
                # 重新采样
                resampled_signal = []
                for i in range(n_bins):
                    start_idx = i * bin_size
                    end_idx = min((i + 1) * bin_size, len(avg_signal))
                    bin_avg = np.mean(avg_signal[start_idx:end_idx])
                    resampled_signal.append(bin_avg)
                avg_signal = resampled_signal
                # 更新位置信息：保持原始的区域长度比例
                positions = np.linspace(region_start, region_end, n_bins)
            else:
                positions = list(range(region_start, region_start + len(avg_signal)))
            
            # 应用高斯平滑进一步提高平滑度
            from scipy import ndimage
            try:
                # 高斯平滑，sigma=1.5 提供适度的平滑
                avg_signal = ndimage.gaussian_filter1d(avg_signal, sigma=1.5)
            except ImportError:
                # 如果没有scipy，使用简单的移动平均
                window_size = 5
                if len(avg_signal) >= window_size:
                    smoothed_signal = []
                    for i in range(len(avg_signal)):
                        start_idx = max(0, i - window_size // 2)
                        end_idx = min(len(avg_signal), i + window_size // 2 + 1)
                        smoothed_signal.append(np.mean(avg_signal[start_idx:end_idx]))
                    avg_signal = smoothed_signal
            
            # 存储结果，包括TSS、TTS位置和基因名称
            profile_results[region_id] = {
                'positions': list(positions),
                'signal': list(avg_signal),
                'tss': tss,
                'tts': tts,
                'center_pos': center_pos,  # 记录实际使用的中心位置
                'gene_name': gene_name
            }
        else:
            print(f"  警告: 区域 {region_id} 没有有效的信号数据")
    
    # 绘制每个区域的平均profile图
    profile_pdf = os.path.join(output_dir, "average_profile.pdf")
    
    if profile_results:
        # 新的绘图逻辑：所有曲线在一张图中
        plt.figure(figsize=(12, 8))
        
        # 使用不同颜色表示不同的基因/区域
        colors = plt.cm.Set3(np.linspace(0, 1, len(profile_results)))
        
        if center == "peak":
            # peak模式：region和gene通用，以peak summit为中心对齐
            for i, (region_id, color) in enumerate(zip(profile_results.keys(), colors)):
                data = profile_results[region_id]
                center_pos = data['center_pos']  # peak summit位置
                gene_name = data['gene_name']
                
                # 计算相对于peak summit的位置
                rel_positions = [p - center_pos for p in data['positions']]
                
                plt.plot(rel_positions, data['signal'], color=color, linewidth=2, 
                        label=f"{gene_name}", alpha=0.8)
            
            # 添加peak summit标记线
            plt.axvline(x=0, color='red', linestyle='--', alpha=0.7, linewidth=2)
            plt.text(0, plt.gca().get_ylim()[1]*0.95, 'Peak Summit', 
                    color='red', fontweight='bold', ha='center')
            plt.xlabel('Distance from Peak Summit (bp)')
            plt.xlim(-5000, 5000)
            plt.title('Average Profile - Peak Summit Aligned')
            
        elif mode == "region":
            if center == "position":
                # region + position: 将所有region的横坐标范围标准化
                for i, (region_id, color) in enumerate(zip(profile_results.keys(), colors)):
                    data = profile_results[region_id]
                    signal = data['signal']
                    
                    # 标准化横坐标到0-100%
                    x_normalized = np.linspace(0, 100, len(signal))
                    
                    plt.plot(x_normalized, signal, color=color, linewidth=2, 
                            label=f"{data['gene_name']}", alpha=0.8)
                
                plt.xlabel('Relative Position (%)')
                plt.title('Average Profile - Position Normalized')
        
        elif mode == "gene":
            if center == "whole":
                # gene + whole: TSS前5kb + 全基因(缩放) + TTS后5kb，比例1:5:1
                total_width = 700  # 总宽度点数
                tss_width = 100    # TSS前5kb部分
                tts_width = 100    # TTS后5kb部分
                gene_width = 500   # 基因主体部分
                
                for i, (region_id, color) in enumerate(zip(profile_results.keys(), colors)):
                    data = profile_results[region_id]
                    gene_name = data['gene_name']
                    
                    signal = data['signal']
                    
                    # 创建三段式坐标
                    tss_x = np.arange(0, tss_width)
                    gene_x = np.arange(tss_width, tss_width + gene_width)
                    tts_x = np.arange(tss_width + gene_width, total_width)
                    
                    # 简化处理：将现有信号分为三段
                    third = len(signal) // 3
                    tss_signal = signal[:third]
                    gene_signal = signal[third:2*third]
                    tts_signal = signal[2*third:]
                    
                    # 重新采样到目标长度
                    from scipy.interpolate import interp1d
                    try:
                        if len(tss_signal) > 1:
                            f_tss = interp1d(np.arange(len(tss_signal)), tss_signal, kind='linear')
                            tss_resampled = f_tss(np.linspace(0, len(tss_signal)-1, tss_width))
                        else:
                            tss_resampled = np.full(tss_width, tss_signal[0] if tss_signal else 0)
                        
                        if len(gene_signal) > 1:
                            f_gene = interp1d(np.arange(len(gene_signal)), gene_signal, kind='linear')
                            gene_resampled = f_gene(np.linspace(0, len(gene_signal)-1, gene_width))
                        else:
                            gene_resampled = np.full(gene_width, gene_signal[0] if gene_signal else 0)
                        
                        if len(tts_signal) > 1:
                            f_tts = interp1d(np.arange(len(tts_signal)), tts_signal, kind='linear')
                            tts_resampled = f_tts(np.linspace(0, len(tts_signal)-1, tts_width))
                        else:
                            tts_resampled = np.full(tts_width, tts_signal[0] if tts_signal else 0)
                    except ImportError:
                        # 如果没有scipy，使用简单的重采样
                        tss_resampled = np.interp(np.linspace(0, len(tss_signal)-1, tss_width), 
                                                np.arange(len(tss_signal)), tss_signal)
                        gene_resampled = np.interp(np.linspace(0, len(gene_signal)-1, gene_width), 
                                                 np.arange(len(gene_signal)), gene_signal)
                        tts_resampled = np.interp(np.linspace(0, len(tts_signal)-1, tts_width), 
                                                np.arange(len(tts_signal)), tts_signal)
                    
                    # 合并三段信号
                    full_signal = np.concatenate([tss_resampled, gene_resampled, tts_resampled])
                    full_x = np.arange(total_width)
                    
                    plt.plot(full_x, full_signal, color=color, linewidth=2, 
                            label=f"{gene_name}", alpha=0.8)
                
                # 添加TSS和TTS标记线
                plt.axvline(x=tss_width, color='red', linestyle='--', alpha=0.7, linewidth=2)
                plt.axvline(x=tss_width + gene_width, color='blue', linestyle='--', alpha=0.7, linewidth=2)
                plt.text(tss_width, plt.gca().get_ylim()[1]*0.95, 'TSS', 
                        color='red', fontweight='bold', ha='center')
                plt.text(tss_width + gene_width, plt.gca().get_ylim()[1]*0.90, 'TTS', 
                        color='blue', fontweight='bold', ha='center')
                plt.xlabel('Relative Position (TSS-5kb | Gene Body | TTS+5kb)')
                plt.title('Average Profile - Whole Gene (1:5:1)')
                
            elif center == "TSS":
                # gene + TSS: TSS前后5kb
                for i, (region_id, color) in enumerate(zip(profile_results.keys(), colors)):
                    data = profile_results[region_id]
                    center_pos = data['center_pos']  # 使用实际的中心位置
                    gene_name = data['gene_name']
                    
                    # 计算相对于中心位置的位置
                    rel_positions = [p - center_pos for p in data['positions']]
                    
                    plt.plot(rel_positions, data['signal'], color=color, linewidth=2, 
                            label=f"{gene_name}", alpha=0.8)
                
                # 添加TSS标记线
                plt.axvline(x=0, color='red', linestyle='--', alpha=0.7, linewidth=2)
                plt.text(0, plt.gca().get_ylim()[1]*0.95, 'TSS', 
                        color='red', fontweight='bold', ha='center')
                plt.xlabel('Distance from TSS (bp)')
                plt.xlim(-5000, 5000)
                plt.title('Average Profile - TSS Centered')
                
            else:  # center == "TTS"
                # gene + TTS: TTS前后5kb  
                for i, (region_id, color) in enumerate(zip(profile_results.keys(), colors)):
                    data = profile_results[region_id]
                    center_pos = data['center_pos']  # 使用实际的中心位置
                    gene_name = data['gene_name']
                    
                    # 计算相对于中心位置的位置
                    rel_positions = [p - center_pos for p in data['positions']]
                    
                    plt.plot(rel_positions, data['signal'], color=color, linewidth=2, 
                            label=f"{gene_name}", alpha=0.8)
                
                # 添加TTS标记线
                plt.axvline(x=0, color='blue', linestyle='--', alpha=0.7, linewidth=2)
                plt.text(0, plt.gca().get_ylim()[1]*0.95, 'TTS', 
                        color='blue', fontweight='bold', ha='center')
                plt.xlabel('Distance from TTS (bp)')
                plt.xlim(-5000, 5000)
                plt.title('Average Profile - TTS Centered')
        
        # 通用图形设置
        plt.ylabel('Average Signal')
        plt.grid(True, alpha=0.3)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(profile_pdf, bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"平均Profile图已生成: {profile_pdf}")
    else:
        print("警告: 没有生成平均Profile图，因为没有有效的区域数据")
    
    # 创建并保存矩阵文件（包含所有区域和样本的信号）
    matrix_tab = os.path.join(output_dir, "matrix.tab")
    
    try:
        # 保存一个简单的矩阵文件，以便后续分析
        with open(matrix_tab, 'w') as f:
            # 写入标题行
            f.write("region_id\tregion_start\tregion_end\tsample_count\tavg_signal\n")
            
            # 写入每个区域的数据
            for region_id in region_ids:
                if region_id in profile_results:
                    chrom, pos = region_id.split(':')
                    start, end = pos.split('-')
                    avg_signal = np.mean(profile_results[region_id]['signal'])
                    sample_count = len(bigwig_files)
                    
                    f.write(f"{region_id}\t{start}\t{end}\t{sample_count}\t{avg_signal}\n")
        
        print(f"矩阵文件已保存: {matrix_tab}")
    except Exception as e:
        print(f"保存矩阵文件时出错: {e}", file=sys.stderr)
    
    return matrix_tab, profile_pdf

# 计算每个区域在每个样本中的平均信号值
def compute_region_signals(regions_df, bigwig_files, sample_names=None):
    """计算每个区域在每个样本中的平均信号值"""
    if not sample_names:
        sample_names = [f"Sample_{i+1}" for i in range(len(bigwig_files))]
    
    # 创建结果DataFrame
    result_df = regions_df.copy()
    result_df['region_id'] = result_df.apply(lambda row: f"{row['chr']}:{row['start']}-{row['end']}", axis=1)
    
    # 为每个样本计算信号
    for i, bw_file in enumerate(bigwig_files):
        sample_name = sample_names[i]
        print(f"计算样本 {sample_name} 的信号...")
        
        try:
            # 打开BigWig文件
            bw = pyBigWig.open(bw_file)
            
            # 计算每个区域的平均信号
            signal_values = []
            for _, row in regions_df.iterrows():
                chrom, start, end = row['chr'], row['start'], row['end']
                
                # 确保区域在有效范围内
                if chrom in bw.chroms() and start < bw.chroms()[chrom] and end > 0:
                    # 调整区域边界
                    valid_start = max(0, start)
                    valid_end = min(bw.chroms()[chrom], end)
                    
                    try:
                        # 获取区域内的信号值
                        values = bw.values(chrom, valid_start, valid_end)
                        # 过滤掉NaN值并计算平均值
                        valid_values = [v for v in values if v is not None and not np.isnan(v)]
                        avg_signal = np.mean(valid_values) if valid_values else 0
                    except:
                        avg_signal = 0
                else:
                    avg_signal = 0
                
                signal_values.append(avg_signal)
            
            # 添加到结果DataFrame
            result_df[sample_name] = signal_values
            
            # 关闭BigWig文件
            bw.close()
            
        except Exception as e:
            print(f"处理样本 {sample_name} 时出错: {e}", file=sys.stderr)
            result_df[sample_name] = 0
    
    return result_df

# 画箱线图
def plot_boxplot(signal_df, output_file, region_id_col='region_id', region_id_mappings=None):
    """绘制boxplot，每个样本对应一个box，包含该样本在所有region/gene中的signal值
    
    Args:
        signal_df: 包含信号值的DataFrame
        output_file: 输出文件路径
        region_id_col: region_id列名
        region_id_mappings: region_id到显示名称的映射字典（可选）
    """
    try:
        # 获取样本列（除了chr、start、end和region_id等非信号列之外的列）
        non_signal_cols = ['chr', 'start', 'end', region_id_col, 'name', 'score', 'strand','tss','tts','gene_start','gene_end']
        sample_cols = [col for col in signal_df.columns if col not in non_signal_cols]
        
        if not sample_cols:
            raise ValueError("未找到样本列进行绘图")
        
        print(f"绘制boxplot，使用样本列: {sample_cols}")
        print(f"区域数量: {len(signal_df)}")
        
        # 准备绘图数据：每个样本对应一个box，包含该样本在所有region中的signal值
        data_list = []
        for sample_col in sample_cols:
            for i, row in signal_df.iterrows():
                if pd.notna(row[sample_col]):  # 排除NaN值
                    try:
                        signal_value = float(row[sample_col])
                        data_list.append({
                            'Sample': sample_col,
                            'Signal': signal_value
                        })
                    except (ValueError, TypeError):
                        print(f"警告: 无法将列 {sample_col} 的值 '{row[sample_col]}' 转换为浮点数，已跳过")
        
        # 创建DataFrame
        plot_data = pd.DataFrame(data_list)
        
        if plot_data.empty:
            print("警告: 没有有效数据用于绘制boxplot")
            plt.figure(figsize=(8, 6))
            plt.text(0.5, 0.5, "No Data Available", 
                     horizontalalignment='center', 
                     verticalalignment='center',
                     fontsize=20)
            plt.savefig(output_file)
            plt.close()
            return
        
        # 设置seaborn的美观风格
        sns.set_style("whitegrid")
        sns.set_context("notebook", font_scale=1.2)
        
        # 使用自定义颜色
        samples = plot_data['Sample'].unique()
        palette = sns.color_palette("Set2", n_colors=len(samples))
        
        # 创建图形和轴
        plt.figure(figsize=(max(10, len(samples) * 1.5), 8))
        
        # 绘制boxplot：横坐标为样本，每个样本的box包含该样本在所有region中的signal值
        ax = sns.boxplot(
            x='Sample', 
            y='Signal',
            data=plot_data, 
            palette=palette,
            width=0.6,               # 箱子宽度
            fliersize=5,             # 异常值点大小
            linewidth=1.5            # 线宽
        )
        
        # 添加数据点以更好地显示分布
        sns.stripplot(
            x='Sample',
            y='Signal',
            data=plot_data, 
            color='black',
            size=3, 
            alpha=0.3,
            jitter=True
        )
        
        # 优化标题和标签
        plt.title('Signal Distribution Across Samples', fontsize=16, pad=20)
        plt.xlabel('Samples', fontsize=14, labelpad=10)
        plt.ylabel('Signal Value', fontsize=14, labelpad=10)
        
        # 优化x轴标签
        plt.xticks(rotation=45, ha='right')
        
        # 添加网格线
        ax.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # 添加边框
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
            spine.set_color('#2F2F2F')
        
        # 保存图形（高DPI以获得更好的质量）
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Boxplot 已保存到: {output_file}")
        
    except Exception as e:
        print(f"绘制 Boxplot 时出错: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

# 为region_list模式绘制boxplot
def plot_boxplot_region_list(signal_df, region_groups, output_file, region_id_col='region_id', sample_cols=None):
    """为region_list模式绘制boxplot
    
    Args:
        signal_df: 包含信号值的DataFrame
        region_groups: 区域组列表，每组包含多个区域/基因
        output_file: 输出文件路径
        region_id_col: region_id列名
        sample_cols: 样本列名列表
    """
    try:
        # 获取样本列（除了chr、start、end和region_id等非信号列之外的列）
        if sample_cols is None:
            non_signal_cols = ['chr', 'start', 'end', region_id_col, 'name', 'score', 'strand','tss','tts','gene_start','gene_end']
            sample_cols = [col for col in signal_df.columns if col not in non_signal_cols]
        
        if not sample_cols:
            raise ValueError("未找到样本列进行绘图")
        
        print(f"绘制region_list模式的boxplot，使用样本列: {sample_cols}")
        print(f"区域组数量: {len(region_groups)}")
        
        # 创建名称到region_id的映射
        name_to_region_id = {}
        if 'name' in signal_df.columns:
            for _, row in signal_df.iterrows():
                name_to_region_id[row['name']] = row[region_id_col]
        
        # 准备绘图数据
        data_list = []
        
        for group in region_groups:
            group_id = group['group_id']
            regions_in_group = group['regions']
            
            # 对每个样本
            for sample_col in sample_cols:
                group_signals = []
                
                # 收集该组内所有区域在当前样本中的信号值
                for region_name in regions_in_group:
                    # 查找对应的region_id
                    region_id = name_to_region_id.get(region_name, region_name)
                    
                    # 在signal_df中查找匹配的行
                    matching_rows = signal_df[
                        (signal_df[region_id_col] == region_id) |
                        (signal_df.get('name', '') == region_name)
                    ]
                    
                    if not matching_rows.empty:
                        signal_value = matching_rows.iloc[0][sample_col]
                        if pd.notna(signal_value):
                            try:
                                group_signals.append(float(signal_value))
                            except (ValueError, TypeError):
                                pass
                
                # 将该组该样本的所有信号值加入绘图数据
                for signal_value in group_signals:
                    data_list.append({
                        'Group': group_id,
                        'Sample': sample_col,
                        'Signal': signal_value
                    })
        
        # 创建DataFrame
        plot_data = pd.DataFrame(data_list)
        
        if plot_data.empty:
            print("警告: 没有有效数据用于绘制region_list模式的boxplot")
            plt.figure(figsize=(8, 6))
            plt.text(0.5, 0.5, "No Data Available", 
                     horizontalalignment='center', verticalalignment='center', fontsize=20)
            plt.savefig(output_file)
            plt.close()
            return
        
        # 设置绘图风格
        sns.set_style("whitegrid")
        sns.set_context("notebook", font_scale=1.2)
        
        # 创建图形
        g = sns.FacetGrid(plot_data, col="Group", col_wrap=2, height=8, aspect=1.2, sharex=True, sharey=True)
        # 绘制 boxplot
        g.map(sns.boxplot, 'Sample', 'Signal', palette="Set2",order=sorted(plot_data['Sample'].unique()), width=0.6, linewidth=1.5)

        # 绘制 stripplot
        g.map(sns.stripplot, 'Sample', 'Signal', palette="Set2",order=sorted(plot_data['Sample'].unique()), size=7, jitter=True)

        # 优化标题和标签
        g.set_axis_labels("Samples", "Signal Value")
        g.set_titles("{col_name}")

        # 只保留最后一行的 x 轴标签

        for ax in g.axes.flat:
            ax.xaxis.set_tick_params(labelbottom=False)  # 默认隐藏所有 x 轴 label

        # 找到最后一行的 axes
        last_row_axes = g.axes[-(len(g.col_names) % g._col_wrap):] if len(g.col_names) % g._col_wrap != 0 else g.axes[-g._col_wrap:]

        for ax in last_row_axes:
            ax.xaxis.set_tick_params(labelbottom=True)  # 显示最后一行的 x 轴 label
            ax.set_xticklabels(sorted(plot_data['Sample'].unique()), rotation=90, ha='right')

        samples = sorted(plot_data['Sample'].unique())
        n_samples = len(samples)

        # 设置颜色（可根据你的 palette 修改）
        colors = sns.color_palette("Set2", n_samples)
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], 
                marker='o', 
                color='w', 
                label=sample,
                markerfacecolor=colors[i], 
                markersize=8,
                linestyle='') 
            for i, sample in enumerate(samples)
        ]

        # 添加图例到画布右侧
        g.figure.legend(
            handles=legend_elements,
            title='Samples',
            loc='center left',
            bbox_to_anchor=(1.02, 0.5)
        ) 
        g.figure.suptitle('Signal Distribution Across Samples by Group', fontsize=16, y=1.05)

        # 保存图形
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Region_list模式Boxplot 已保存到: {output_file}")
        
    except Exception as e:
        print(f"绘制 Region_list模式Boxplot 时出错: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

def main():
    # 解析命令行参数
    args = parse_arguments()
    with open(args.config_path, "r") as f:
        config = json.load(f)

    if args.species!= None:
        gtf = config['data'][args.species]['refseq_gtf']
        if args.mode == "gene" and not os.path.exists(gtf):
            print(f"基因注释GTF文件不存在: {gtf}")
            sys.exit(1)

    # 获取输入文件列表
    with open(args.bigwig, 'r') as file:
        bigwig_files = [line.strip() for line in file if line.strip()]

   # 检查文件是否存在
    nonexistent_files = [file_path for file_path in bigwig_files if not os.path.exists(file_path)]
    if nonexistent_files:
        print("以下文件不存在:")
        for file_path in nonexistent_files:
            print(file_path)
    bigwig_files = [file_path for file_path in bigwig_files if os.path.exists(file_path)]
    if not bigwig_files:
        print("错误：输入文件列表为空或所有文件不存在！请检查输入文件。")
        sys.exit(1) 

    # 生成样本名称
    sample_names = [os.path.basename(bw).split('.')[0] for bw in bigwig_files]
    
    # 根据模式获取分析区域
    if args.mode == "region":
        # 区域模式：直接读取区域文件
        regions_df = read_region_list(args.input)
    else:
        # 基因模式：提取启动子区域
        genes = read_gene_list(args.input)
        regions_df = get_gene_promoter_regions(genes, gtf, args.promoter_size)
    
    # 创建输出目录
    profile_dir = os.path.join(args.output, "profile")
    os.makedirs(profile_dir, exist_ok=True)
    
    # 计算矩阵并生成average profile图
    matrix_file, profile_pdf = compute_matrix_and_profile(regions_df, bigwig_files, profile_dir, args.mode, args.center)
    print(f"平均Profile图已生成: {profile_pdf}")
    
    # 计算每个区域在每个样本中的信号值
    signal_df = compute_region_signals(regions_df, bigwig_files, sample_names)
    
    # 保存信号值表格
    signal_file = os.path.join(args.output, "region_signals.tsv")
    signal_df.to_csv(signal_file, sep='\t', index=False)
    print(f"区域信号值已保存到: {signal_file}")
    
    # 绘制boxplot
    boxplot_file = os.path.join(args.output, "region_signals_boxplot.png")
    
    if args.method == "sample":
        # sample模式：每个样本一个box，包含该样本在所有region中的signal值
        if args.mode == "gene" and 'name' in regions_df.columns:
            region_to_gene = {}
            for _, row in regions_df.iterrows():
                region_id = f"{row['chr']}:{row['start']}-{row['end']}"
                gene_name = row['name']
                region_to_gene[region_id] = gene_name
            plot_boxplot(signal_df, boxplot_file, region_id_mappings=region_to_gene)
        else:
            plot_boxplot(signal_df, boxplot_file)
    else:
        # region_list模式：新逻辑
        region_groups = read_region_list_file(args.region_list)
        # 获取样本列
        non_signal_cols = ['chr', 'start', 'end', 'region_id', 'name', 'score', 'strand','tss','tts','gene_start','gene_end']
        sample_cols = [col for col in signal_df.columns if col not in non_signal_cols]
        plot_boxplot_region_list(signal_df, region_groups, boxplot_file, sample_cols=sample_cols)
    
    print("分析完成!")

if __name__ == "__main__":
    main() 

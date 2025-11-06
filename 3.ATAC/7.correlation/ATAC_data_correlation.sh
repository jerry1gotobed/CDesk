#!/bin/bash
# 检查参数数量
if [ "$#" -lt 2 ] || [ "$#" -gt 6 ]; then
    echo "Usage: $0 <bw_file1> <bw_file2> <output_directory> <threads> [region_file] <config>"
    echo "  <bw_file1>: Path to the first bw file (e.g., GSM7789744.bw)"
    echo "  <bw_file2>: Path to the second bw file (e.g., GSM7789752.bw)"
    echo "  <output_directory>: Output directory"
    echo "  <threads>: Number of threads (default:30)"
    echo "  [region_file]: Optional BED file for specific regions (e.g., GSM7789744_e5_peaks.bed)"
    exit 1
fi

# 获取文件前缀
bw1=$(echo $1|sed -e "s/.*\///g" -e "s/\.bw$//g")
bw2=$(echo $2|sed -e "s/.*\///g" -e "s/\.bw$//g")
OUTPUT_DIRECTORY=$(echo $3)
THREAD=$(echo $4)
if [ "$#" -gt 5 ]; then
    region=$(echo $5|sed -e "s/.*\///g" -e "s/\.bed$//g" -e "s/\.bedgraph//g")
fi

# 定义输出文件路径
out_npz="$OUTPUT_DIRECTORY/${bw1}_vs_${bw2}.npz"
out_txt="$OUTPUT_DIRECTORY/${bw1}_vs_${bw2}.txt"
out_peak_npz="$OUTPUT_DIRECTORY/${bw1}_vs_${bw2}_${region}.npz"
out_peak_txt="$OUTPUT_DIRECTORY/${bw1}_vs_${bw2}_${region}.txt"

# 检查 bw 文件是否存在
if [ ! -f "$1" ] || [ ! -f "$2" ]; then
    echo "Error: One or both of the bw files do not exist."
    exit 1
fi

# 检查 region 文件是否存在（如果指定）
if [ "$#" -gt 5 ]; then
    if [ -n "$5" ] && [ ! -f "$5" ]; then
        echo "Error: The specified region file does not exist."
        exit 1
    fi
fi

# 检查 OUTPUT_DIRECTORY 是否存在
if [ ! -d "$OUTPUT_DIRECTORY" ]; then
    mkdir -p "$OUTPUT_DIRECTORY"
    echo "Directory $OUTPUT_DIRECTORY created."
fi

SCRIPT_DIR=$(dirname $(realpath $0))

# 执行 multiBigwigSummary 命令
if [ -z "$region" ]; then
    # 默认 bins 模式
    multiBigwigSummary bins -b $1 $2 --labels $bw1 $bw2  -p $THREAD -o $out_npz --outRawCounts $out_txt
    sed -i "s/#//g;s/'//g" $out_txt
    Rscript ${SCRIPT_DIR}/atac.correlationplot.R $out_txt
    echo "Bins mode completed. Output files: $out_npz and $out_txt"
else
    # 使用指定的 region 文件
    multiBigwigSummary bins -b $1 $2 --labels $bw1 $bw2 -p $THREAD -o $out_peak_npz --BED $5  --outRawCounts $out_peak_txt
    sed -i "s/#//g;s/\\'//g" $out_peak_txt
    Rscript ${SCRIPT_DIR}/atac.correlationplot.R $out_peak_txt
    echo "Region mode completed. Output files: $out_peak_npz and $out_peak_txt"
fi

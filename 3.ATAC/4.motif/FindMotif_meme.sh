#!/bin/bash

# 默认参数
config_json=${!#}

default_meme_p=10          # 默认并行线程数
default_meme_nmotifs=3     # 默认提取的 motif 数量

if [ -z "$config_json" ]; then
    echo "错误: 必须提供 JSON 配置文件作为最后一个参数！"
    exit 1

fi

if [ ! -f "$config_json" ]; then
    echo "错误: JSON 配置文件不存在！路径: $config_json"
    exit 1
fi

memechip=$(jq -r '.software.memechip' $config_json)
if [ $? -ne 0 ] || [ -z "$memechip" ]; then
    echo "错误: 无法从 JSON 文件中解析 'software.memechip'！请检查文件内容。"
    exit 1
fi

if [ ! -x "$memechip" ]; then
    echo "错误: 找不到 meme-chip 可执行文件！路径: $memechip"
    exit 1
fi


# 解析输入参数

while getopts "i:o:s:p:n:" opt; do
    case $opt in
        i) input_dir=$OPTARG ;;       # 输入目录路径
        o) output_dir=$OPTARG ;;      # 输出目录路径
        s) species=$OPTARG ;;         # 选择物种：human 或 mouse
        p) meme_p=$OPTARG ;;          # 用户指定的并行线程数
        n) meme_nmotifs=$OPTARG ;;    # 用户指定的提取 motif 数量
        *) echo "Usage: $0 -i input_directory -o output_directory -s human|mouse -p meme_p -n meme_nmotifs"
           exit 1
           ;;
    esac
done

# 检查输入参数

if [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$species" ]; then
    echo "错误: 参数缺失！"
    echo "正确用法: $0 -i input_directory -o output_directory -s human|mouse -p meme_p -n meme_nmotifs"
    exit 1
fi

# 设置 motif 数据库目录
motif_database_dir=$(jq -r '.motif_database' "$config_json")

# 根据物种选择数据库

case $species in

    human)
        motif_file="${motif_database_dir}/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant.meme"
        extra_motif_file="${motif_database_dir}/HOCOMOCO/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
        ;;
    mouse)
        motif_file="${motif_database_dir}/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant.meme"
        extra_motif_file="${motif_database_dir}/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme"
        ;;
    *)
        echo "错误: 物种参数必须是 human 或 mouse"
        exit 1

        ;;
esac

# 检查输入输出目录

if [ ! -d "$input_dir" ]; then

    echo "输入目录 $input_dir 不存在！"
    exit 1

fi

if [ ! -d "$output_dir" ]; then

    echo "输出目录 $output_dir 不存在，正在创建..."
    mkdir -p "$output_dir"
fi

# 如果未指定 meme-p 和 meme-nmotifs，则使用默认值
meme_p=${meme_p:-$default_meme_p}
meme_nmotifs=${meme_nmotifs:-$default_meme_nmotifs}
echo $meme_nmotifs

# 遍历并处理每个 .fa 或 .fasta 文件
find "$input_dir" -type f \( -name "*.fa" -o -name "*.fasta" \) | while read id; do

    file=$(basename "$id" | sed 's/\.[^.]*$//')  # 去除路径和后缀（适配 .fa 和 .fasta）
    echo "正在处理样本: $file"

    $memechip -meme-p "$meme_p" \
        -oc "${output_dir}/${file}.results/" \
        -db "$motif_file" \
        -db "$extra_motif_file" \
        "$id" -meme-nmotifs "$meme_nmotifs"

done

echo "所有运行已完成！"

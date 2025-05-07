#!/bin/bash

# Description : generate TE
# Author      : WEI SHI
# Version     : 1.0
# Time        : 2024-10-20 12:28:00

show_help() {
    echo "Usage: bash BulkRNAseqTE.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (.bam), the ending does not require a '/'"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -p INT		Specify number of threads (default is 8)"
    echo "  -s SPECIES          Specify the species"
    echo "  -t TE_NAME_FILE     Specify the path to TE name reference file (optional)"
}

config_json=${!#}
THREAD=8
SCRIPT_DIR=$(dirname $(realpath $0))
ALLOWED_SPECIES=("mm10" "hg38" "rn7" "susScr11" "galGal6")

SCTE=$(jq -r '.software.scTE' $config_json)
for tool in "$SCTE"; do
    if [ ! -x "$tool" ]; then
        echo "Error: Tool not found or not executable: $tool"
        exit 1
    fi
done

while getopts ":hi:s:o:t:p:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        i )
            INPUT_DIRECTORY=$(realpath "${OPTARG}")
            ;;
        o)
            OUTPUT_DIRECTORY=$(realpath "${OPTARG}")
            ;;
        p )
            THREAD=${OPTARG}
            ;;
	s )
            SPECIES=${OPTARG}
	    if [[ ! " ${ALLOWED_SPECIES[@]} " =~ " ${SPECIES} " ]]; then
                echo "Error: Invalid species '${SPECIES}'. Allowed species are: ${ALLOWED_SPECIES[*]}"
                exit 1
            fi
            ;;
	t )
            TE_NAME_FILE=$(realpath "${OPTARG}")
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            show_help
            exit 1
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" 1>&2
            show_help
            exit 1
            ;;
    esac
done

# 检查是否提供了必要的参数
if [ -z "$INPUT_DIRECTORY" ] || [ -z "$OUTPUT_DIRECTORY" ]; then
    echo "Error: INPUT_DIRECTORY and OUTPUT_DIRECTORY must be provided."
    show_help
    exit 1
fi

if [ ! -d "${OUTPUT_DIRECTORY}" ]; then
    mkdir ${OUTPUT_DIRECTORY}
fi

CONFIG=$(cat $config_json | jq -r --arg species "$SPECIES" '.data[$species]')
REF_IDX=$(echo "$CONFIG" | jq -r '.TE_idx')
if [ -z "$REF_IDX" ]; then
    echo "Error: 'REF_IDX' is not set or empty in the CONFIG."
    exit 1
fi
if [ ! -f "$REF_IDX" ]; then
    echo "Error: File '$REF_IDX' does not exist."
    exit 1
fi

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

echo "---------------------------------- 参考数据使用 ${REF_IDX} --------------------------------"

counter=0
for FILE in ${INPUT_DIRECTORY}/*.bam; do
    ((counter++))
    echo "----------------------------- 这是第 ${counter} 个样本~~~~~~~ ---------------------------"

    echo "132 `get_time` scTE ..."

    cd ${OUTPUT_DIRECTORY}
    SAMPLE=$(basename ${FILE} .bam)
    $SCTE -i ${FILE} -o ${SAMPLE} -p ${THREAD} -x ${REF_IDX} --hdf5 False -CB False -UMI False
    cd -
done

echo "--------------------------开始合并scTE结果-------------------------------"
# 合并所有scTE结果
mkdir ${OUTPUT_DIRECTORY}/merge

# 获取第一个csv文件名
first_file=$(ls ${OUTPUT_DIRECTORY}/*.csv | head -n 1)
# 复制第一个文件
cp "$first_file" ${OUTPUT_DIRECTORY}/merge/combined.csv

# 循环处理其余csv文件
for file in ${OUTPUT_DIRECTORY}/*.csv; do
  if [ "$file" != "$first_file" ]; then
    sed -n '2p' "$file" >> ${OUTPUT_DIRECTORY}/merge/combined.csv
  fi
done

# 转置
python ${SCRIPT_DIR}/transpose_csv.py ${OUTPUT_DIRECTORY}/merge/combined.csv ${OUTPUT_DIRECTORY}/merged_scTE.csv
rm -rf ${OUTPUT_DIRECTORY}/merge

# select TE data
if [ -n "${TE_NAME_FILE}" ] && [ -f "${TE_NAME_FILE}" ]; then
  head -n 1 ${OUTPUT_DIRECTORY}/merged_scTE.csv >> ${OUTPUT_DIRECTORY}/TE_data.csv
  awk -F ',' 'NR==FNR {te[$1]; next} $1 in te' "${TE_NAME_FILE}" ${OUTPUT_DIRECTORY}/merged_scTE.csv > ${OUTPUT_DIRECTORY}/TE_tem.csv
  cat ${OUTPUT_DIRECTORY}/TE_tem.csv >> ${OUTPUT_DIRECTORY}/TE_filter.csv
  rm -rf ${OUTPUT_DIRECTORY}/TE_tem.csv
fi


echo "--------------------------合并完成，可以查看结果了--------------------------"


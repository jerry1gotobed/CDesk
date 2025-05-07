#!/bin/bash

# Description : generate BW
# Author      : WEI SHI
# Version     : 1.0
# Time        : 2024-10-20 12:01:00

show_help() {
    echo "Usage: bash BulkRNAseqBW.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (.bam), the ending does not require a '/'"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -p INT		Specify number of threads (default is 8)"
}

config_json=${!#}
THREAD=8
BAMCOVERAGE=$(jq -r '.software.bamCoverage' $config_json)
for tool in "$BAMCOVERAGE"; do
    if [ ! -x "$tool" ]; then
        echo "Error: Tool not found or not executable: $tool"
        exit 1
    fi
done

while getopts ":hi:o:p:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        i )
            INPUT_DIRECTORY=${OPTARG}
            ;;
        o)
            OUTPUT_DIRECTORY=${OPTARG}
            ;;
        p )
            THREAD=${OPTARG}
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

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

counter=0
for FILE in ${INPUT_DIRECTORY}/*.bam; do
    ((counter++))

    echo "这是第 ${counter} 个样本~~~~~~~~~"
    echo "----------------------------- Start generating BW files for ${FILE} ---------------------------"
    echo "128 `get_time` bamCoverage ..."

    SAMPLE=$(basename ${FILE} .bam)
    $BAMCOVERAGE -p ${THREAD} -b ${FILE} -o ${OUTPUT_DIRECTORY}/${SAMPLE}.bw
done


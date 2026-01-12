#!/usr/bin/env bash
show_help() {
    echo "Usage: bash BulkRNAseqTE.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (.bam), the ending does not require a '/'"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -p INT		Specify number of threads (default is 8)"
    echo "  -s SPECIES          Specify the species"
}

config_json=$CDesk_config
THREAD=8
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")

if command -v scTE > /dev/null 2>&1; then
    SCTE_CMD="scTE"
else
    if [[ -f $config_json ]]; then
        SCTE_PATH=$(jq -r '.software.scTE' $config_json 2>/dev/null)
        if [[ -n "$SCTE_PATH" && "$SCTE_PATH" != "null" ]]; then
            if [[ -x "$SCTE_PATH" ]]; then
                SCTE_CMD="$SCTE_PATH"
            else
                echo "❌ Error: scTE path from config.json is not executable: $SCTE_PATH" >&2
                exit 1
            fi
        else
            echo "❌ Error: 'scTE' not found or empty in config.json" >&2
            exit 1
        fi
    else
        echo "❌ Error: config.json not found" >&2
        exit 1
    fi
fi

while getopts ":hi:s:o:p:" opt; do
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
	s )
            SPECIES=${OPTARG}
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

# Check parameters
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

QC_path=$(dirname "${OUTPUT_DIRECTORY}")/Log

echo " "
echo "----------------------------- Start scTE analysis---------------------------"
echo " "

echo "---------------------------------- Reference data ${REF_IDX} --------------------------------"

counter=0
for FILE in ${INPUT_DIRECTORY}/*.bam; do
    ((counter++))
    echo "`get_time` scTE ..."

    cd ${OUTPUT_DIRECTORY}
    SAMPLE=$(basename ${FILE} .bam)
    echo "----------------------------- Number ${counter} sample: ${SAMPLE} ---------------------------"
    scTE -i ${FILE} -o ${SAMPLE} -p ${THREAD} -x ${REF_IDX} --hdf5 False -CB False -UMI False > ${QC_path}/${SAMPLE}.scTE.log 2>&1
    cd - > /dev/null
done

echo "--------------------------Merge scTE result-------------------------------"
mkdir ${OUTPUT_DIRECTORY}/merge

# The first csv
first_file=$(ls ${OUTPUT_DIRECTORY}/*.csv | head -n 1)
cp "$first_file" ${OUTPUT_DIRECTORY}/merge/combined.csv

# Traverse remaining csv files
for file in ${OUTPUT_DIRECTORY}/*.csv; do
  if [ "$file" != "$first_file" ]; then
    sed -n '2p' "$file" >> ${OUTPUT_DIRECTORY}/merge/combined.csv
  fi
done

# Transpose
$python3_7 ${SCRIPT_DIR}/transpose_csv.py ${OUTPUT_DIRECTORY}/merge/combined.csv ${OUTPUT_DIRECTORY}/merged_scTE.csv
rm -rf ${OUTPUT_DIRECTORY}/merge

echo "--------------------------scTE done--------------------------"

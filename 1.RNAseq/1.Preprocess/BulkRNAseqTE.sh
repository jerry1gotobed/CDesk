#!/usr/bin/env bash
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

config_json=$CDesk_config
THREAD=8
SCRIPT_DIR=$(dirname $(realpath $0))

tools=("scTE")
missing_tools=()
echo "Checking required tools..."
for tool in "${tools[@]}"; do
    if ! command -v "$tool" &>/dev/null; then
        missing_tools+=("$tool")
    fi
done
if [ ${#missing_tools[@]} -gt 0 ]; then
    echo "Error: The following required tools are not installed or not in PATH:" >&2
    printf ' - %s\n' "${missing_tools[@]}" >&2
    exit 1
fi

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

echo "---------------------------------- Reference data ${REF_IDX} --------------------------------"

counter=0
for FILE in ${INPUT_DIRECTORY}/*.bam; do
    ((counter++))
    echo "----------------------------- Number ${counter} sample~~~~~~~ ---------------------------"

    echo "132 `get_time` scTE ..."

    cd ${OUTPUT_DIRECTORY}
    SAMPLE=$(basename ${FILE} .bam)
    scTE -i ${FILE} -o ${SAMPLE} -p ${THREAD} -x ${REF_IDX} --hdf5 False -CB False -UMI False
    cd -
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
python3.7 ${SCRIPT_DIR}/transpose_csv.py ${OUTPUT_DIRECTORY}/merge/combined.csv ${OUTPUT_DIRECTORY}/merged_scTE.csv
rm -rf ${OUTPUT_DIRECTORY}/merge

# select TE data
if [ -n "${TE_NAME_FILE}" ] && [ -f "${TE_NAME_FILE}" ]; then
  head -n 1 ${OUTPUT_DIRECTORY}/merged_scTE.csv >> ${OUTPUT_DIRECTORY}/TE_tem.csv
  awk -F ',' 'NR==FNR {te[$1]; next} $1 in te' "${TE_NAME_FILE}" ${OUTPUT_DIRECTORY}/merged_scTE.csv > ${OUTPUT_DIRECTORY}/TE_tem.csv
  cat ${OUTPUT_DIRECTORY}/TE_tem.csv >> ${OUTPUT_DIRECTORY}/TE_filter.csv
  rm -rf ${OUTPUT_DIRECTORY}/TE_tem.csv
fi


echo "--------------------------scTE done--------------------------"

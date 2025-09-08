#!/usr/bin/env bash
show_help() {
    echo "Usage: bash BulkRNAseqBW.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (.bam), the ending does not require a '/'"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -p INT		Specify number of threads (default is 8)"
}

config_json=$CDesk_config
THREAD=8

tools=("bamCoverage")
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

# Parameters check
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

    echo "Number ${counter} sample~~~~~~~~~"
    echo "----------------------------- Start generating BW files for ${FILE} ---------------------------"
    echo "128 `get_time` bamCoverage ..."

    SAMPLE=$(basename ${FILE} .bam)
    bamCoverage -p ${THREAD} -b ${FILE} -o ${OUTPUT_DIRECTORY}/${SAMPLE}.bw
done


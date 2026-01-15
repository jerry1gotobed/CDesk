#!/usr/bin/env bash
config_json=$CDesk_config

# Check whether the tools are available
# Software required
tools=("ascp")
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

show_help() {
    echo "Usage: Download_SRR.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -k Key_file		ASCP Key Path"
    echo "  -i INPUT_DIRECTORY	Specify the input file"
    echo "  -o OUTPUT_DIRECTORY	Specify the output directory"
    echo "  -p INT		Specify number of threads (default is 4)"
}

MAX_CONCURRENT=4
while getopts ":hk:i:o:p:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        k )
            KEY=${OPTARG}
            ;;
        i )
            INPUT_FILE=${OPTARG}
            ;;
        o)
            OUTPUT_DIRECTORY=${OPTARG}
            ;;
        p )
            MAX_CONCURRENT=${OPTARG}
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

mkdir -p $OUTPUT_DIRECTORY
openssh=$KEY

num_srr_err=$(grep -c "^\(SRR\|ERR\)" "$INPUT_FILE" 2>/dev/null || echo 0)
echo "$num_srr_err SRR/ERR samples in total"

if [ "$num_srr_err" -eq 0 ]; then
    echo "❌ Error: No SRR/ERR IDs found in $INPUT_FILE. Please check the input file." >&2
    exit 1
fi

while IFS= read -r id || [ -n "$id" ]; do
    
    # Skip empty line
    [[ -z "$id" ]] && continue
    # Delete spaces
    id="${id// /}"
    # SRR/ERR starts only
    if [[ ! "$id" =~ ^(SRR|ERR) ]]; then
        echo "⚠️  Skipping non-SRR/ERR ID: $id"
        continue
    fi

    num=$(echo -n "$id" | wc -m)
    echo "📥 Downloading $id"

    while [ $(jobs -r | wc -l) -ge "$MAX_CONCURRENT" ]; do
        sleep 5
    done

    if [ $num -eq 12 ]; then
        x=$(echo "$id" | cut -b 1-6)
        y=$(echo "$id" | cut -b 10-11)
        ascp -QT -l 500m -P33001 -k 1 -i "$openssh" \
            era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/0$y/$id/ $OUTPUT_DIRECTORY &

    elif [ $num -eq 11 ]; then
        x=$(echo "$id" | cut -b 1-6)
        y=$(echo "$id" | cut -b 10-11)
        ascp -QT -l 500m -P33001 -k 1 -i "$openssh" \
            era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/0$y/$id/ $OUTPUT_DIRECTORY &

    elif [ $num -eq 10 ]; then
        x=$(echo "$id" | cut -b 1-6)
        y=$(echo "$id" | cut -b 10-10)
        ascp -QT -l 500m -P33001 -k 1 -i "$openssh" \
            era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/00$y/$id/ $OUTPUT_DIRECTORY &

    elif [ $num -eq 9 ]; then
        x=$(echo "$id" | cut -b 1-6)
        ascp -QT -l 500m -P33001 -k 1 -i "$openssh" \
            era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/$id/ $OUTPUT_DIRECTORY &

    else
        echo "❌ Unsupported format: $id (length = $num)"
    fi

done < "$INPUT_FILE"

echo "⏳ Waiting for all downloads to complete..."
wait
echo "✅ All downloads completed!"

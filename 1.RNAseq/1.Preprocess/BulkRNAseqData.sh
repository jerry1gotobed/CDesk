#!/usr/bin/env bash
config_json=$CDesk_config
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")

show_help() {
    echo "Usage: BulkRNAseqData.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -s SPECIES		Specify the species"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (fastq.gz), the ending does not require a '/'"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -p INT		Specify number of threads (default is 8)"
}

# Check whether the tools are available
# Software required
tools=("jq" "hisat2" "samtools" "trim_galore" "fastqc" "multiqc" "gfold")
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
else
    echo "All required tools are available."
fi

THREAD=8

while getopts ":hs:i:o:p:l:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        i )
            INPUT_DIRECTORY=${OPTARG}
            ;;
        s )
            SPECIES=${OPTARG}
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
if [ -z "$INPUT_DIRECTORY" ] || [ -z "$SPECIES" ] || [ -z "$OUTPUT_DIRECTORY" ]; then
    echo "Error: INPUT_DIRECTORY, OUTPUT_DIRECTORY and SPECIES must be provided."
    show_help
    exit 1
fi

CONFIG=$(cat $config_json | jq -r --arg species "$SPECIES" '.data[$species]')
if [ -z "$CONFIG" ]; then
    echo "Error: No configuration found for species '$SPECIES' in data.json"
    exit 1
fi

# Assign the reference data
mapping_index=$(echo "$CONFIG" | jq -r '.mapping_index')
mapping_gtf=$(echo "$CONFIG" | jq -r '.refseq_gtf')
gfold_gtf=$(echo "$CONFIG" | jq -r '.refseq_gtf') 
RSeQC_bed=$(echo "$CONFIG" | jq -r '.refseq_bed')
# variables check
if [ -z "$mapping_index" ]; then
    echo "Error: 'mapping_index' is not set or empty in the CONFIG."
    exit 1
fi
if [ -z "$mapping_gtf" ]; then
    echo "Error: 'refseq_gtf' is not set or empty in the CONFIG."
    exit 1
fi
if [ -z "$gfold_gtf" ]; then
    echo "Error: 'gfold_gtf' is not set or empty in the CONFIG."
    exit 1
fi
if [ -z "$RSeQC_bed" ]; then
    echo "Error: 'refseq_bed' is not set or empty in the CONFIG."
    exit 1
fi

if ! ls "${mapping_index}"*.ht* 1> /dev/null 2>&1; then
  echo "HISAT2 index files are missing. Please check the index path."
fi
if [ ! -f "$mapping_gtf" ]; then
    echo "Error: File '$mapping_gtf' does not exist."
    exit 1
fi
if [ ! -f "$gfold_gtf" ]; then
    echo "Error: File '$gfold_gtf' does not exist."
    exit 1
fi
if [ ! -f "$RSeQC_bed" ]; then
    echo "Error: File '$RSeQC_bed' does not exist."
    exit 1
fi

# Grab the fq files name
file_list_single=$(ls "$INPUT_DIRECTORY"/*_single.{fastq,fq}.gz 2>/dev/null | sed 's/_single.fastq.gz//;s/_single.fq.gz//' | sort -u)
file_list=$(ls "$INPUT_DIRECTORY"/*_{1,2}.{fastq,fq}.gz 2>/dev/null | sed 's/_1.fastq.gz//;s/_2.fastq.gz//;s/_1.fq.gz//;s/_2.fq.gz//' | sort -u)

log_path="${OUTPUT_DIRECTORY}/Log"
bam_path="${OUTPUT_DIRECTORY}/Bam"
QC_path="${OUTPUT_DIRECTORY}/QC"
Expression_path="${OUTPUT_DIRECTORY}/Expression"
gfold_path="${Expression_path}/gfold"

if [ ! -d "${OUTPUT_DIRECTORY}" ]; then
    mkdir -p ${OUTPUT_DIRECTORY}
fi
if [ ! -d "${log_path}" ]; then
    mkdir -p ${log_path}
fi
if [ ! -d "${bam_path}" ]; then
    mkdir -p ${bam_path}
fi
if [ ! -d "${QC_path}" ]; then
    mkdir -p ${QC_path}
fi
if [ ! -d "${Expression_path}" ]; then
    mkdir -p ${Expression_path}
fi
if [ ! -d "${gfold_path}" ]; then
    mkdir -p ${gfold_path}
fi

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

cp ${INPUT_DIRECTORY}/*.bam ${bam_path}

echo "--------------------------------------------INITIALIZING----------------------------------------------" 
echo "RNA-seq data analysis pipeline is now running..."
echo "Number of threads ---------- ${THREAD}" 
echo "Directory of data ---------- ${INPUT_DIRECTORY}"
echo "Directory of result ---------- ${OUTPUT_DIRECTORY}"
echo "Mapping index ---------- ${mapping_index}"
echo "Mapping gtf ---------- ${mapping_gtf}"
echo "RSeQC bed ---------- ${RSeQC_bed}"

####################################################################################################
#################################          ANALYSIS STEPS           ################################
####################################################################################################
echo "---------------------------------------Process fq files-----------------------------------------------------"
counter=0

touch ${log_path}/fastqc_multiqc.log 
echo '' > ${log_path}/fastqc_multiqc.log 

# Paired
for FILE in $file_list; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename ${FILE})
    BAM_FILE="${bam_path}/${SAMPLE_PREFIX}.bam"
    echo " "
    echo "----------------------------Number ${counter} fq sample: ${SAMPLE_PREFIX}--------------------------"

    if [ -f "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" ] && samtools quickcheck -q "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam"; then
        echo "Available BAM file checked，skip mapping"
        if ln "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" "${bam_path}/${SAMPLE_PREFIX}.bam" 2>/dev/null; then
            echo "Linked ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam -> ${bam_path}/${SAMPLE_PREFIX}.bam"
        else
            cp "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" "${bam_path}/${SAMPLE_PREFIX}.bam"
            echo "Copy ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam -> ${bam_path}/${SAMPLE_PREFIX}.bam"
        fi
        SKIP_MAPPING=true
        BAM_FILE="${bam_path}/${SAMPLE_PREFIX}.bam"
    else
        echo "No available BAM file checked，do mapping"
        SKIP_MAPPING=false
    fi

    FILE_NAME=""
    R1_FILE_NAME=""
    R2_FILE_NAME=""

    # Check according fq file
    if [ -e "${FILE}_1.fastq.gz" ] && [ -e "${FILE}_2.fastq.gz" ]; then
        R1_FILE_NAME="${FILE}_1.fastq.gz"
        R2_FILE_NAME="${FILE}_2.fastq.gz"
    FILE_TYPE="fastq"
    elif [ -e "${FILE}_1.fq.gz" ] && [ -e "${FILE}_2.fq.gz" ]; then
        R1_FILE_NAME="${FILE}_1.fq.gz"
        R2_FILE_NAME="${FILE}_2.fq.gz"
    FILE_TYPE="fq"
    else
        echo "Error: Can not find ${SAMPLE_REFIX} 2 fq files, skip"
        continue
    fi

    echo " "
    echo "`get_time` fastqc ..."
    if [[ "${FILE_TYPE}" == "fastq" ]]; then
        fastqc -o ${QC_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fastq.gz >> ${log_path}/fastqc_multiqc.log 2>&1
    elif [[ "${FILE_TYPE}" == "fq" ]]; then
        fastqc -o ${QC_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fq.gz >> ${log_path}/fastqc_multiqc.log 2>&1
    fi

    if [ "$SKIP_MAPPING" = false ]; then
        # 1.mapping ----------------------------------------------------------------------------------------
        echo " "
        echo "Mapping ..."
        echo " "

        echo "`get_time` trim_galore ..."
        echo " "
        trim_galore --gzip -j ${THREAD} --paired ${R1_FILE_NAME} ${R2_FILE_NAME} --trim-n -o ${OUTPUT_DIRECTORY} --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1

        echo "`get_time` hisat2 ..."
        echo " "
	hisat2 -p ${THREAD} --dta -x ${mapping_index} -1 ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_R1_val_1.fq.gz -2 ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_R2_val_2.fq.gz -S ${bam_path}/${SAMPLE_PREFIX}.align.sam --summary-file ${QC_path}/${SAMPLE_PREFIX}_mapping_summary.txt > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1
        rm -rf ${OUTPUT_DIRECTORY}/*.fq.gz

        echo "`get_time` samtools view ..."
        echo " "
        samtools view -@ ${THREAD} -bS ${bam_path}/${SAMPLE_PREFIX}.align.sam > ${bam_path}/${SAMPLE_PREFIX}.bam
	rm ${bam_path}/${SAMPLE_PREFIX}.align.sam

    else
        echo "Use available BAM: $BAM_FILE"
    	echo " "
    fi

done

# Single
for FILE in $file_list_single; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename ${FILE})
    BAM_FILE="${bam_path}/${SAMPLE_PREFIX}.bam"
    echo " "
    echo "----------------------------Number ${counter} fq sample: ${SAMPLE_PREFIX}--------------------------"

    if [ -f "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" ] && samtools quickcheck -q "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam"; then
        echo "Available BAM file checked，skip mapping"
        if ln "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" "${bam_path}/${SAMPLE_PREFIX}.bam" 2>/dev/null; then
            echo "Linked ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam -> ${bam_path}/${SAMPLE_PREFIX}.bam"
        else
            cp "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" "${bam_path}/${SAMPLE_PREFIX}.bam"
            echo "Copy ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam -> ${bam_path}/${SAMPLE_PREFIX}.bam"
        fi
        SKIP_MAPPING=true
        BAM_FILE="${bam_path}/${SAMPLE_PREFIX}.bam"
    else
        echo "No available BAM file checked，do mapping"
        SKIP_MAPPING=false
    fi

    FILE_NAME=""
    R1_FILE_NAME=""
    R2_FILE_NAME=""

    # Check according fq file
    if [ -e "${FILE}.fastq.gz" ] || [ -e "${FILE}.fq.gz" ]; then
        if [ -e "${FILE}.fastq.gz" ]; then
            FILE_NAME="${FILE}.fastq.gz"
        else
            FILE_NAME="${FILE}.fq.gz"
        fi
    else
        echo "Error：Can not find ${SAMPLE_REFIX} sample fq fiile, skip"
        continue
    fi

    echo " "
    echo "`get_time` fastqc ..."
    if [[ "${FILE_TYPE}" == "fastq" ]]; then
        fastqc -o ${QC_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fastq.gz >> ${log_path}/fastqc_multiqc.log 2>&1
    elif [[ "${FILE_TYPE}" == "fq" ]]; then
        fastqc -o ${QC_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fq.gz >> ${log_path}/fastqc_multiqc.log 2>&1
    fi

    if [ "$SKIP_MAPPING" = false ]; then
        # 1.mapping ----------------------------------------------------------------------------------------
        echo " "
        echo "Mapping ..."
        echo " "

        echo "`get_time` trim_galore ..."
        echo " "
        trim_galore --gzip -j ${THREAD} ${FILE_NAME} --trim-n -o ${OUTPUT_DIRECTORY} --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1

        echo "`get_time` hisat2 ..."
        echo " "
        hisat2 -p ${THREAD} --dta -x ${mapping_index} -U ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_trimmed.fq.gz -S ${bam_path}/${SAMPLE_PREFIX}.align.sam --summary-file ${QC_path}/${SAMPLE_PREFIX}_mapping_summary.txt > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1 
        rm -rf ${OUTPUT_DIRECTORY}/*.fq.gz

        echo "`get_time` samtools view ..."
        echo " "
        samtools view -@ ${THREAD} -bS ${bam_path}/${SAMPLE_PREFIX}.align.sam > ${bam_path}/${SAMPLE_PREFIX}.bam
	    rm ${bam_path}/${SAMPLE_PREFIX}.align.sam

    else
        echo "Use available BAM: $BAM_FILE"
    	echo " "
    fi
done

multiqc $QC_path --outdir $QC_path >> ${log_path}/fastqc_multiqc.log 2>&1

echo "---------------------------------------Process bam file----------------------------------------------------"
counter=0
file_list=$(ls "$bam_path"/*.bam 2>/dev/null | sort -u)
for FILE in $file_list; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename "$FILE" .bam)
    echo " "
    echo "----------------------------Number ${counter} bam sample: ${SAMPLE_PREFIX}--------------------------"
    echo " "
    # Check sort status
    if ! samtools view -H "${bam_path}/${SAMPLE_PREFIX}.bam" | grep -q "SO:coordinate"; then
        echo "Sort BAM..."
        tmp_sorted=$(mktemp)
        samtools sort -@ $THREAD "${bam_path}/${SAMPLE_PREFIX}.bam" -o "$tmp_sorted"
        mv "$tmp_sorted" "${bam_path}/${SAMPLE_PREFIX}.bam"
    fi

    # Check index status
    if [ ! -f "${bam_path}/${SAMPLE_PREFIX}.bam.bai" ]; then
        echo "BAM index..."
        samtools index -@ $THREAD "${bam_path}/${SAMPLE_PREFIX}.bam"
    fi
    
    echo "`get_time` RSeQC ..."
    echo " "
    $python3_6 ${SCRIPT_DIR}/read_distribution.py -i ${bam_path}/${SAMPLE_PREFIX}.bam -r ${RSeQC_bed} > ${QC_path}/${SAMPLE_PREFIX}_dirstribution.txt

    echo "Calculate gene expression levels"
    echo "`get_time` gfold count ..."
    echo " "
    samtools view -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam | gfold count -ann ${gfold_gtf} -tag stdin -o ${gfold_path}/${SAMPLE_PREFIX}.read_cnt > ${log_path}/${SAMPLE_PREFIX}.gfold.log 2>&1
done

echo "`get_time` Merge expression matrix ..."
# -------------------------------------------Merge fpkm----------------------------------------------
# Gram the first column of the first .read_cnt to a csv file
awk '{print $1}' $(ls ${gfold_path}/*.read_cnt | head -n 1) >> ${Expression_path}/output1.csv

# Traverse all .read_cnt files
for file in ${gfold_path}/*.read_cnt; do
    paste -d ',' ${Expression_path}/output1.csv <(awk '{print $5}' "$file") > ${Expression_path}/temp1.csv
    mv ${Expression_path}/temp1.csv ${Expression_path}/output1.csv
done

# Add header
echo "gene_name,$(ls ${gfold_path}/*.read_cnt | xargs -n 1 basename | sed 's/.read_cnt//;s/^/,/' | tr -d '\n' | sed 's/^,//')" > ${Expression_path}/temp1.csv
cat ${Expression_path}/output1.csv >> ${Expression_path}/temp1.csv

# Get merged_fpkm
mv ${Expression_path}/temp1.csv ${Expression_path}/merged_fpkm.csv
rm -rf ${Expression_path}/output1.csv

# -------------------------------------------Merge count----------------------------------------------
# Gram the first column of the first .read_cnt to a csv file
awk '{print $1}' $(ls ${gfold_path}/*.read_cnt | head -n 1) >> ${Expression_path}/output1.csv

# Traverse all .read_cnt files
for file in ${gfold_path}/*.read_cnt; do
    paste -d ',' ${Expression_path}/output1.csv <(awk '{print $3}' "$file") > ${Expression_path}/temp1.csv
    mv ${Expression_path}/temp1.csv ${Expression_path}/output1.csv
done

# Add header
echo "gene_name,$(ls ${gfold_path}/*.read_cnt | xargs -n 1 basename | sed 's/.read_cnt//;s/^/,/' | tr -d '\n' | sed 's/^,//')" > ${Expression_path}/temp1.csv
cat ${Expression_path}/output1.csv >> ${Expression_path}/temp1.csv

# Get merged_count
mv ${Expression_path}/temp1.csv ${Expression_path}/merged_count.csv
rm -rf ${Expression_path}/output1.csv
echo " "
echo "Expression matrix has been merged..."
echo "Preprocess Done"

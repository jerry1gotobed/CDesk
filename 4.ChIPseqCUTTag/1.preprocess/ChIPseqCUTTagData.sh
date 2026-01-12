#!/usr/bin/env bash
config_json=$CDesk_config
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")

# Check whether the tools are available
# Software required
tools=("jq" "bowtie2" "samtools" "trim_galore" "fastqc" "multiqc" "macs" "bedClip" "bedGraphToBigWig" "macs3")
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

show_help() {
    echo "Usage: ATACseqData.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -s SPECIES		Specify the species"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (fastq.gz)"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file"
    echo "  -p INT		Specify number of threads (default is 8)"
    echo "  -g GROUP		Grouping file"
}

THREAD=8

while getopts ":hs:i:o:p:g:" opt; do
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
        o )
            OUTPUT_DIRECTORY=${OPTARG}
            ;;
        p )
            THREAD=${OPTARG}
            ;;
	    g )
            Group=${OPTARG}
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

# Check necessary parameters
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
mapping_index=$(echo "$CONFIG" | jq -r '.bowtie2_mapindex')
species=$(echo "$CONFIG" | jq -r '.effective_genome_size')
chromInfo=$(echo "$CONFIG" | jq -r '.chromInfo')
promoter_file=$(echo "$CONFIG" | jq -r '.promoter_file')

if ! ls "${mapping_index}"*.bt* 1> /dev/null 2>&1; then
  echo "bowtie2 index files are missing. Please check the index path ${mapping_index}."
fi
if [ -z "$species" ]; then
    echo "Error: Effective genome size of {$SPECIES} is not set or empty in the CONFIG."
    exit 1
fi
if [ ! -f "$chromInfo" ]; then
    echo "Error: chromInfo file ${chromInfo} does not exist."
    exit 1
fi
if [ ! -f "$promoter_file" ]; then
    echo "Error: promter file ${promoter_file} does not exist."
    exit 1
fi

# Grab the file names
file_list_single=$(ls "$INPUT_DIRECTORY"/*_single.{fastq,fq}.gz 2>/dev/null | sed 's/_single.fastq.gz//;s/_single.fq.gz//' | sort -u)
file_list=$(ls "$INPUT_DIRECTORY"/*_{1,2}.{fastq,fq}.gz 2>/dev/null | sed 's/_1.fastq.gz//;s/_2.fastq.gz//;s/_1.fq.gz//;s/_2.fq.gz//' | sort -u)

log_path="${OUTPUT_DIRECTORY}/Log"
bam_path="${OUTPUT_DIRECTORY}/Bam"
bw_path="${OUTPUT_DIRECTORY}/Bw"
peak_path="${OUTPUT_DIRECTORY}/Peak"
qc_path="${OUTPUT_DIRECTORY}/QC"
signal_path="${OUTPUT_DIRECTORY}/Signal"
fastqc_path="${OUTPUT_DIRECTORY}/QC/fastqc"
temp_path="${OUTPUT_DIRECTORY}/Temp"
bam_path_single="${OUTPUT_DIRECTORY}/Bam/single"
bam_path_paired="${OUTPUT_DIRECTORY}/Bam/pair"
fragment_path="${OUTPUT_DIRECTORY}/QC/fragment"

if [ ! -d "${OUTPUT_DIRECTORY}" ]; then
    mkdir -p ${OUTPUT_DIRECTORY}
fi
if [ ! -d "${log_path}" ]; then
    mkdir -p ${log_path}
fi
if [ ! -d "${bam_path}" ]; then
    mkdir -p ${bam_path}
fi
if [ ! -d "${bam_path_single}" ]; then
    mkdir -p ${bam_path_single}
fi
if [ ! -d "${bam_path_paired}" ]; then
    mkdir -p ${bam_path_paired}
fi
if [ ! -d "${bw_path}" ]; then
    mkdir -p ${bw_path}
fi
if [ ! -d "${peak_path}" ]; then
    mkdir -p ${peak_path}
fi
if [ ! -d "${qc_path}" ]; then
    mkdir -p ${qc_path}
fi
if [ ! -d "${signal_path}" ]; then
    mkdir -p ${signal_path}
fi
if [ ! -d "${fastqc_path}" ]; then
    mkdir -p ${fastqc_path}
fi
if [ ! -d "${temp_path}" ]; then
    mkdir -p ${temp_path}
fi
if [ ! -d "${fragment_path}" ]; then
    mkdir -p ${fragment_path}
fi

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

echo "--------------------------------------------INITIALIZING----------------------------------------------" 
echo "ChIPseqCUTTag data analysis pipeline is now running..."
echo "Number of threads ---------- ${THREAD}" 
echo "Directory of data ---------- ${INPUT_DIRECTORY}"
echo "Directory of result ---------- ${OUTPUT_DIRECTORY}"
echo "Mapping index ---------- ${mapping_index}"
echo "ChromInfo ---------- ${chromInfo}"
echo "PromoterInfo ---------- ${promoter_file}"

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
    
    echo "-----------------------------------------------------------------------------------------------------"
    echo "----------------------------Number ${counter} fq sample: ${SAMPLE_PREFIX}--------------------------"
    echo "-----------------------------------------------------------------------------------------------------"

    if [ -f "${bam_path_paired}/${SAMPLE_PREFIX}.bam" ] && samtools quickcheck -q "${bam_path_paired}/${SAMPLE_PREFIX}.bam"; then
        echo "Available BAM file checked，skip mapping"
        SKIP_MAPPING=true
        BAM_FILE="${bam_path_paired}/${SAMPLE_PREFIX}.bam"
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

    echo "fastqc `get_time` fastqc ..."
    if [[ "${FILE_TYPE}" == "fastq" ]]; then
        fastqc -o ${qc_path}/fastqc ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fastq.gz >> ${log_path}/fastqc_multiqc.log 2>&1
    elif [[ "${FILE_TYPE}" == "fq" ]]; then
        fastqc -o ${qc_path}/fastqc ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fq.gz >> ${log_path}/fastqc_multiqc.log 2>&1
    fi
    
    if [ "$SKIP_MAPPING" = false ]; then
        # 1.mapping ----------------------------------------------------------------------------------------
        echo " "
        echo "1.mapping..."
        echo " "

        # trim
        echo "paired data `get_time` trim_galore ..."
        echo " "
        trim_galore --gzip -j ${THREAD} --paired ${R1_FILE_NAME} ${R2_FILE_NAME} --trim-n -o ${temp_path} --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1

        # map
        echo "paired data `get_time` bowtie2 ..."
        echo " "
	bowtie2 -q -p ${THREAD} -x ${mapping_index} -1 ${temp_path}/${SAMPLE_PREFIX}_R1_val_1.fq.gz -2 ${temp_path}/${SAMPLE_PREFIX}_R2_val_2.fq.gz -S ${bam_path_paired}/${SAMPLE_PREFIX}.align.sam > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1
        rm -rf ${temp_path}/${SAMPLE_PREFIX}_R1_val_1.fq.gz ${temp_path}/${SAMPLE_PREFIX}_R2_val_2.fq.gz

        # sam2bam
        echo "sam2bam `get_time` samtools view ..."
        echo " "
        samtools view -@ ${THREAD} -bS ${bam_path_paired}/${SAMPLE_PREFIX}.align.sam > ${bam_path_paired}/${SAMPLE_PREFIX}.bam
    else
        echo "Use available BAM: $BAM_FILE"
    fi

done

# Single
for FILE in $file_list_single; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename ${FILE})
    
    echo "-----------------------------------------------------------------------------------------------------"
    echo "----------------------------Number ${counter} fq sample: ${SAMPLE_PREFIX}--------------------------"
    echo "-----------------------------------------------------------------------------------------------------"

    if [ -f "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" ] && samtools quickcheck -q "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam"; then
        echo "Available BAM file checked，skip mapping"
        if ln "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" "${bam_path_single}/${SAMPLE_PREFIX}.bam" 2>/dev/null; then
            echo "Linked ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam -> ${bam_path_single}/${SAMPLE_PREFIX}.bam"
        else
            cp "${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam" "${bam_path_single}/${SAMPLE_PREFIX}.bam"
            echo "Copy ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}.bam -> ${bam_path_single}/${SAMPLE_PREFIX}.bam"
        fi
        SKIP_MAPPING=true
        BAM_FILE="${bam_path_single}/${SAMPLE_PREFIX}.bam"
    else
        echo "No available BAM file checked，do mapping"
        SKIP_MAPPING=false
    fi

    FILE_NAME=""
    R1_FILE_NAME=""
    R2_FILE_NAME=""

    # Check according fq file
    if [ -e "${FILE}_single.fastq.gz" ] || [ -e "${FILE}_single.fq.gz" ]; then
        if [ -e "${FILE}_single.fastq.gz" ]; then
            FILE_NAME="${FILE}_single.fastq.gz"
	    FILE_TYPE="fastq"
        else
            FILE_NAME="${FILE}_single.fq.gz"
	    FILE_TYPE="fq"
        fi
    else
        echo "Error：Can not find ${SAMPLE_PREFIX} sample fq fiile, skip"
        continue
    fi

    echo "fastqc `get_time` fastqc ..."
    if [[ "${FILE_TYPE}" == "fastq" ]]; then
        fastqc -o ${qc_path}/fastqc ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fastq.gz >> ${log_path}/fastqc_multiqc.log 2>&1
    elif [[ "${FILE_TYPE}" == "fq" ]]; then
        fastqc -o ${qc_path}/fastqc ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fq.gz >> ${log_path}/fastqc_multiqc.log 2>&1
    fi
    
    if [ "$SKIP_MAPPING" = false ]; then
        # 1.mapping ----------------------------------------------------------------------------------------
        echo " "
        echo "1.mapping..."
        echo " "

        # trim
        echo "single data `get_time` trim_galore ..."
        echo " "
        trim_galore --gzip -j ${THREAD} ${FILE_NAME} --trim-n -o ${temp_path}/ --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1

        # map
        echo "single data `get_time` bowtie2 ..."
        echo " "
        bowtie2 -q -p ${THREAD} -x ${mapping_index} -U ${temp_path}/${SAMPLE_PREFIX}_trimmed.fq.gz -S ${bam_path_single}/${SAMPLE_PREFIX}.align.sam > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1
        rm -rf ${temp_path}/${SAMPLE_PREFIX}_trimmed.fq.gz

        # sam2bam
        echo "sam2bam `get_time` samtools view ..."
        echo " "
        samtools view -@ ${THREAD} -bS ${bam_path_single}/${SAMPLE_PREFIX}.align.sam > ${bam_path_single}/${SAMPLE_PREFIX}.bam
    else
        echo "Use available BAM: $BAM_FILE"
    fi

done

multiqc ${qc_path}/fastqc --outdir ${qc_path}/fastqc >> ${log_path}/fastqc_multiqc.log 2>&1

echo "---------------------------------------Process bam file----------------------------------------------------"
# Single
file_list=$(ls "$bam_path_single"/*.bam 2>/dev/null | sort -u)
for FILE in $file_list; do
    SAMPLE_PREFIX=$(basename "$FILE" .bam)
    echo "----------------------------Process ${SAMPLE_PREFIX}--------------------------"
    echo " "

    # Sort and index
    echo "Sort BAM..."
    tmp_sorted=$(mktemp)
    samtools sort -@ $THREAD "${bam_path_single}/${SAMPLE_PREFIX}.bam" -o "$tmp_sorted"
    mv "$tmp_sorted" "${bam_path_single}/${SAMPLE_PREFIX}.bam"
    echo "Index BAM..."
    samtools index -@ $THREAD "${bam_path_single}/${SAMPLE_PREFIX}.bam"
done

# Paired
file_list=$(ls "$bam_path_paired"/*.bam 2>/dev/null | sort -u)
for FILE in $file_list; do
    SAMPLE_PREFIX=$(basename "$FILE" .bam)
    echo "----------------------------Process ${SAMPLE_PREFIX}--------------------------"
    echo " "
    # sort bam
    echo "bam2sortbam `get_time` samtools sort ..."
    echo " "
    
    # Sort and index
    echo "Sort BAM..."
    tmp_sorted=$(mktemp)
    samtools sort -@ $THREAD "${bam_path_paired}/${SAMPLE_PREFIX}.bam" -o "$tmp_sorted"
    mv "$tmp_sorted" "${bam_path_paired}/${SAMPLE_PREFIX}.bam"
    echo "Index BAM..."
    samtools index -@ $THREAD "${bam_path_paired}/${SAMPLE_PREFIX}.bam"
done

if [ -f "$Group" ]; then
    echo " "
    echo "Merge replicates..."
    echo " "
    $python3_7 ${SCRIPT_DIR}/Merge_replicate.py $Group $bam_path
fi

# Single
file_list=$(ls "$bam_path_single"/*.bam 2>/dev/null | sort -u)
counter=0
file_list=$(ls "$bam_path_single"/*.bam 2>/dev/null | sort -u)
for FILE in $file_list; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename "$FILE" .bam)
    samtools view -@ ${THREAD} ${bam_path_single}/${SAMPLE_PREFIX}.bam | awk '/^[^#]/ {print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${fragment_path}/${SAMPLE_PREFIX}.ChromosomeDistribution.txt
    samtools view -@ ${THREAD} ${bam_path_single}/${SAMPLE_PREFIX}.bam | awk '$5 >= 30 && substr($3,1,3)=="chr" && $3!="chrM" && $3!="chrEBV" && $3 !~ /_/ {print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${fragment_path}/${SAMPLE_PREFIX}.Filtered.ChromosomeDistribution.txt
    echo "----------------------------Number ${counter} bam sample: ${SAMPLE_PREFIX}--------------------------"
    echo "single data `get_time` Remove duplicate reads ..."
    # Remove duplicate reads
    samtools view -@ ${THREAD} -H ${bam_path_single}/${SAMPLE_PREFIX}.bam > ${bam_path_single}/${SAMPLE_PREFIX}.Filtered.sam
    samtools view -@ ${THREAD} ${bam_path_single}/${SAMPLE_PREFIX}.bam | awk '$5 >= 30 && substr($3,1,3)=="chr" && $3!="chrM" && $3!="chrEBV" && $3 !~ /_/' >> ${bam_path_single}/${SAMPLE_PREFIX}.Filtered.sam
    samtools view -@ ${THREAD} -bS ${bam_path_single}/${SAMPLE_PREFIX}.Filtered.sam > ${bam_path_single}/${SAMPLE_PREFIX}.Filtered.bam
    rm ${bam_path_single}/${SAMPLE_PREFIX}.Filtered.sam

    samtools sort -@ ${THREAD} -o ${bam_path_single}/${SAMPLE_PREFIX}.FilteredSorted.bam ${bam_path_single}/${SAMPLE_PREFIX}.Filtered.bam
    rm ${bam_path_single}/${SAMPLE_PREFIX}.Filtered.bam
    samtools markdup -@ ${THREAD} -r ${bam_path_single}/${SAMPLE_PREFIX}.FilteredSorted.bam ${bam_path_single}/${SAMPLE_PREFIX}.bam
    rm ${bam_path_single}/${SAMPLE_PREFIX}.FilteredSorted.bam
    
    echo "single data `get_time` Generate bigwig file ..."
    # Generate bigwig file
    macs3 callpeak -f BAM -g ${species} -B -q 0.001 --keep-dup all --SPMR -t ${bam_path_single}/${SAMPLE_PREFIX}.bam -n ${temp_path}/${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.macs3.log 2>&1
    bedClip ${temp_path}/${SAMPLE_PREFIX}_treat_pileup.bdg ${chromInfo} ${temp_path}/${SAMPLE_PREFIX}.tmp.bdg
    grep -v _ ${temp_path}/${SAMPLE_PREFIX}.tmp.bdg | sort -k1,1 -k2,2n  > ${temp_path}/${SAMPLE_PREFIX}.clean.bdg
    bedGraphToBigWig ${temp_path}/${SAMPLE_PREFIX}.clean.bdg ${chromInfo} ${bw_path}/${SAMPLE_PREFIX}.bw
    
    echo "single data `get_time` Peak calling ..."
    # Peak calling
    macs callpeak -g ${species} -p 1e-5 --keep-dup 1 --nomodel --shiftsize 73 -t ${bam_path_single}/${SAMPLE_PREFIX}.bam -n ${peak_path}/${SAMPLE_PREFIX}_e5 > ${log_path}/${SAMPLE_PREFIX}.macs_e5.log 2>&1

    echo "single data `get_time` Promoter signal calculation ..."
    # Promoter signal calculation
    $python3_7 ${SCRIPT_DIR}/GetBigwigSignal.py ${promoter_file} ${bw_path}/${SAMPLE_PREFIX}.bw ${signal_path}/${SAMPLE_PREFIX}.PromoterSignal.csv
done

# paired
file_list=$(ls "$bam_path_paired"/*.bam 2>/dev/null | sort -u)
for FILE in $file_list; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename "$FILE" .bam)
    samtools view -@ ${THREAD} ${bam_path_paired}/${SAMPLE_PREFIX}.bam | awk '/^[^#]/ {print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${fragment_path}/${SAMPLE_PREFIX}.ChromosomeDistribution.txt
    samtools view -@ ${THREAD} -f 0x2 ${bam_path_paired}/${SAMPLE_PREFIX}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3 !~ /_/) print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${fragment_path}/${SAMPLE_PREFIX}.Filtered.ChromosomeDistribution.txt
    echo "----------------------------Number ${counter} bam sample: ${SAMPLE_PREFIX}--------------------------"
    echo "paired data `get_time` Remove duplicate reads ..."
    # Remove duplicate reads
    samtools view -@ ${THREAD} -H ${bam_path_paired}/${SAMPLE_PREFIX}.bam > ${bam_path_paired}/${SAMPLE_PREFIX}.Filtered.sam
    samtools view -@ ${THREAD} -f 0x2 ${bam_path_paired}/${SAMPLE_PREFIX}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3!="chrM" && $3!="chrEBV" && $3 !~ /_/) print forward"\n"$0}' >> ${bam_path_paired}/${SAMPLE_PREFIX}.Filtered.sam
    samtools view -@ ${THREAD} -bS ${bam_path_paired}/${SAMPLE_PREFIX}.Filtered.sam > ${bam_path_paired}/${SAMPLE_PREFIX}.Filtered.bam
    rm ${bam_path_paired}/${SAMPLE_PREFIX}.Filtered.sam

    samtools sort -n -@ ${THREAD} ${bam_path_paired}/${SAMPLE_PREFIX}.Filtered.bam > ${bam_path_paired}/${SAMPLE_PREFIX}.FilteredSorted.bam
    rm ${bam_path_paired}/${SAMPLE_PREFIX}.Filtered.bam
    samtools fixmate -@ ${THREAD} -m ${bam_path_paired}/${SAMPLE_PREFIX}.FilteredSorted.bam ${bam_path_paired}/${SAMPLE_PREFIX}.FixMate.bam
    rm ${bam_path_paired}/${SAMPLE_PREFIX}.FilteredSorted.bam
    samtools sort -@ ${THREAD} -o ${bam_path_paired}/${SAMPLE_PREFIX}.FixMateSorted.bam ${bam_path_paired}/${SAMPLE_PREFIX}.FixMate.bam
    rm ${bam_path_paired}/${SAMPLE_PREFIX}.FixMate.bam
    samtools markdup -@ ${THREAD} -r ${bam_path_paired}/${SAMPLE_PREFIX}.FixMateSorted.bam ${bam_path_paired}/${SAMPLE_PREFIX}.bam
    rm ${bam_path_paired}/${SAMPLE_PREFIX}.FixMateSorted.bam
    samtools view -F 0x04 ${bam_path_paired}/${SAMPLE_PREFIX}.bam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > ${fragment_path}/${SAMPLE_PREFIX}.FragmentsLength.txt
    Rscript ${SCRIPT_DIR}/FragmentLengthDistribution.R ${fragment_path}/${SAMPLE_PREFIX}.FragmentsLength.txt ${fragment_path}

    echo "paired data `get_time` Generate bigwig file ..."
    # Generate bigwig file
    macs3 callpeak -f BAMPE -g ${species} -B -q 0.001 --keep-dup all --SPMR -t ${bam_path_paired}/${SAMPLE_PREFIX}.bam -n ${temp_path}/${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.macs3.log 2>&1
    bedClip ${temp_path}/${SAMPLE_PREFIX}_treat_pileup.bdg ${chromInfo} ${temp_path}/${SAMPLE_PREFIX}.tmp.bdg
    grep -v _ ${temp_path}/${SAMPLE_PREFIX}.tmp.bdg | sort -k1,1 -k2,2n  > ${temp_path}/${SAMPLE_PREFIX}.clean.bdg
    bedGraphToBigWig ${temp_path}/${SAMPLE_PREFIX}.clean.bdg ${chromInfo} ${bw_path}/${SAMPLE_PREFIX}.bw

    echo "paired data `get_time` Peak calling ..."
    # Peak calling
    macs callpeak -g ${species} -p 1e-5 --keep-dup 1 --nomodel --shiftsize 73 -t ${bam_path_paired}/${SAMPLE_PREFIX}.bam -n ${peak_path}/${SAMPLE_PREFIX}_e5 > ${log_path}/${SAMPLE_PREFIX}.macs_e5.log 2>&1

    echo "paired data `get_time` Promoter signal calculation ..."
    # Promoter signal calculation
    $python3_7 ${SCRIPT_DIR}/GetBigwigSignal.py ${promoter_file} ${bw_path}/${SAMPLE_PREFIX}.bw ${signal_path}/${SAMPLE_PREFIX}.PromoterSignal.csv
done

rm -rf ${temp_path}
rm tmp*
echo "ChIPseqCUTTag preprocessing has been finished!!!"

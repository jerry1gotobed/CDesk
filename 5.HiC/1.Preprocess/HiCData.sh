#!/usr/bin/env bash
config_json=$CDesk_config

show_help() {
    echo "Usage: HiCData.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -s SPECIES		Specify the species"
    echo "  -i INPUT_DIRECTORY	Specify the input directory of input file (fastq.gz)"
    echo "  -o OUTPUT_DIRECTORY	Specify the output directory of output file"
    echo "  -t INT		Specify number of threads (default is 8)"
    echo "  -l 1 or 2		1:Single sequencing, 2:Pair sequencing (default is 2)"
    echo "  -f		Use strict filtering"
    echo "  -r		Specify the resolutions"
}

# Default parameters
THREAD=8
LAYOUT=2
Dofilter=False
RESOLUTIONS=10000,100000,1000000

# Software required
tools=("jq" "cooler" "pairtools" "bgzip" "pairix" "bwa" "samtools")
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

# Parameters
while getopts ":hs:i:o:t:l:f:r:z:" opt; do
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
        t )
            THREAD=${OPTARG}
            ;;
        l )
            LAYOUT=${OPTARG}
            if [[ "$LAYOUT" != "1" && "$LAYOUT" != "2" ]]; then
                echo "Invalid value for -l: $LAYOUT. It must be 1 or 2." 1>&2
                show_help
                exit 1
            fi
            ;;
        f )
            Dofilter=True
            ;;
        r )
            RESOLUTIONS=${OPTARG}
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

# Check the required resolutions
if [ -z "$INPUT_DIRECTORY" ] || [ -z "$SPECIES" ] || [ -z "$OUTPUT_DIRECTORY" ]; then
    echo "Error: INPUT_DIRECTORY, OUTPUT_DIRECTORY and SPECIES must be provided."
    show_help
    exit 1
fi

CONFIG=$(cat $config_json | jq -r --arg species "$SPECIES" '.data[$species]')
mapping_index=$(echo "$CONFIG" | jq -r '.fasta')
chromInfo=$(echo "$CONFIG" | jq -r '.chromInfo')

if [ ! -f "$mapping_index" ]; then
    echo "Error: Mapping fasta file not found: $mapping_index"
    exit 1
fi

if [ ! -f "$chromInfo" ]; then
    echo "Error: Chromosome info file not found: $chromInfo"
    exit 1
fi

IFS=',' read -ra RES_ARRAY <<< "$RESOLUTIONS"
for resolution in "${RES_ARRAY[@]}"; do
    # Validate resolution is a number
    if ! [[ "$resolution" =~ ^[0-9]+$ ]]; then
        echo "Invalid resolution value: $resolution. Must be a number." 1>&2
        exit 1
    fi
done

# 提取文件名
if [ "$LAYOUT" = "1" ]; then
    file_list=$(ls "$INPUT_DIRECTORY"/*.{fastq,fq}.gz 2>/dev/null | sed 's/.fastq.gz//;s/.fq.gz//' | sort -u)
elif [ "$LAYOUT" = "2" ]; then
    file_list=$(ls "$INPUT_DIRECTORY"/*.{fastq,fq}.gz 2>/dev/null | sed 's/_1.fastq.gz//;s/_2.fastq.gz//;s/_1.fq.gz//;s/_2.fq.gz//' | sort -u)
fi

if [ -z "$file_list" ]; then
    echo "Error: No valid FASTQ files found in directory '$INPUT_DIRECTORY'."
    echo "Please make sure the files end with .fastq.gz / .fq.gz if single-end; _1.fastq.gz,_2.fastq.gz, / _1.fq.gz,_2.fq.gz if paired-end "
    exit 1
fi

pairsam_path="${OUTPUT_DIRECTORY}/Pairsam"
QC_path="${OUTPUT_DIRECTORY}/QC"
fastqc_path="${OUTPUT_DIRECTORY}/QC/fastqc"
stats_path="${OUTPUT_DIRECTORY}/QC/stats"
pairs_path="${OUTPUT_DIRECTORY}/Pairs"
cool_path="${OUTPUT_DIRECTORY}/Cool"

if [ ! -d "${OUTPUT_DIRECTORY}" ]; then
    mkdir -p ${OUTPUT_DIRECTORY}
fi
if [ ! -d "${QC_path}" ]; then
    mkdir -p ${QC_path}
fi
if [ ! -d "${fastqc_path}" ]; then
    mkdir -p ${fastqc_path}
fi
if [ ! -d "${pairsam_path}" ]; then
    mkdir -p ${pairsam_path}
fi
if [ ! -d "${stats_path}" ]; then
    mkdir -p ${stats_path}
fi
if [ ! -d "${pairs_path}" ]; then
    mkdir -p ${pairs_path}
fi
if [ ! -d "${cool_path}" ]; then
    mkdir -p ${cool_path}
fi

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

####################################################################################################
#################################          ANALYSIS STEPS           ################################
####################################################################################################
counter=0

for FILE in $file_list; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename ${FILE})
    echo "-----------------------------------------------------------------------------------------------------"
    echo "----------------------------第 ${counter} 个样本，样本名为 ${SAMPLE_PREFIX}--------------------------"
    echo "-----------------------------------------------------------------------------------------------------"

    FILE_NAME=""
    R1_FILE_NAME=""
    R2_FILE_NAME=""

    if [ "${LAYOUT}" == "1" ]; then
        if [ -e "${FILE}.fastq.gz" ] || [ -e "${FILE}.fq.gz" ]; then
            if [ -e "${FILE}.fastq.gz" ]; then
                FILE_NAME="${FILE}.fastq.gz"
                FILE_TYPE="fastq"
            else
                FILE_NAME="${FILE}.fq.gz"
                FILE_TYPE="fq"
            fi
        else
            echo "Error: Can not find ${SAMPLE_PREFIX} fq file, skip"
            continue
        fi
    elif [ "${LAYOUT}" == "2" ]; then
        if [ -e "${FILE}_1.fastq.gz" ] && [ -e "${FILE}_2.fastq.gz" ]; then
            R1_FILE_NAME="${FILE}_1.fastq.gz"
            R2_FILE_NAME="${FILE}_2.fastq.gz"
            FILE_TYPE="fastq"
        elif [ -e "${FILE}_1.fq.gz" ] && [ -e "${FILE}_2.fq.gz" ]; then
            R1_FILE_NAME="${FILE}_1.fq.gz"
            R2_FILE_NAME="${FILE}_2.fq.gz"
            FILE_TYPE="fq"
        else
            echo "Error: Can not find the 2 ${SAMPLE_PREFIX} fq files,skip"
            continue
        fi
    fi

    echo " "
    echo "--------------------------------------------INITIALIZING----------------------------------------------" 
    echo "HiC data analysis pipeline is now running..."
    echo "Number of threads ---------- ${THREAD}" 
    echo "Directory of data ---------- ${INPUT_DIRECTORY}"
    echo "Directory of result ---------- ${OUTPUT_DIRECTORY}"
    echo "Name of sample ---------- ${SAMPLE_PREFIX}"
    echo "File of data ---------- ${FILE_NAME}"
    echo "File of R1 data ---------- ${R1_FILE_NAME}"
    echo "File of R2 data ---------- ${R2_FILE_NAME}"
    echo "Method of sequencing ---------- ${LAYOUT}"
    echo "Mapping index ---------- ${mapping_index}"
    echo "ChromInfo ---------- ${chromInfo}"
    echo " "

    echo " "
    echo "-----------------------------------GENERATING DIRECTORY-------------------------------------" 
    echo "Please wait, the directory for storing the results is being generated..."
    echo "Store pairsam files for all samples ---------- ${pairsam_path}"
    echo "Store stats result for all samples ---------- ${stats_path}"
    echo "Store pairs result for all samples ---------- ${pairs_path}"
    echo "Store cool result for all samples ---------- ${cool_path}"
    echo " "

    echo " "
    echo "---------------------------------------START ANALYZING-----------------------------------------" 
    echo "Preparation work completed, start analysis..."

    # fastqc ----------------------------------------------------------------------------------------
    echo ""
    echo "QC `get_time` fastqc ..."
    if [[ "${FILE_TYPE}" == "fastq" ]]; then
        fastqc -o ${fastqc_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fastq.gz
    elif [[ "${FILE_TYPE}" == "fq" ]]; then
        fastqc -o ${fastqc_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fq.gz
    fi

    # mapping ----------------------------------------------------------------------------------------
    echo " "
    echo "paired data `get_time` mapping ..."
    echo " "
    if [ "${LAYOUT}" == "1" ]; then
        if [ "$Dofilter" = False ]; then
            bwa mem -SP5M -t ${THREAD} \
                ${mapping_index} \
                ${FILE_NAME} > ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam 2>/dev/null
            
            cat ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam | \
                samtools view -bhS - | \
                samtools view -h - | \
                pairtools parse -c ${chromInfo} -o ${pairsam_path}/${SAMPLE_PREFIX}.pairsam;
            
            rm ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam
        elif [ "$Dofilter" = True ]; then
            bwa mem -SP5M -t ${THREAD} \
                ${mapping_index} \
                ${FILE_NAME} > ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam 2>/dev/null
            
            cat ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam | \
                samtools view -bhS - | \
                samtools view -h - | \
                pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 -c ${chromInfo} \
                -o ${pairsam_path}/${SAMPLE_PREFIX}.pairsam;
            
            rm ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam
        fi
    elif [ "${LAYOUT}" == "2" ]; then
        if [ "$Dofilter" = False ]; then
            bwa mem -SP5M -t ${THREAD} \
                ${mapping_index} \
                ${R1_FILE_NAME} \
                ${R2_FILE_NAME} > ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam 2>/dev/null
            
            cat ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam | \
                samtools view -bhS - | \
                samtools view -h - | \
                pairtools parse -c ${chromInfo} -o ${pairsam_path}/${SAMPLE_PREFIX}.pairsam;
            
            rm ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam
        elif [ "$Dofilter" = True ]; then
            bwa mem -SP5M -t ${THREAD} \
                ${mapping_index} \
                ${R1_FILE_NAME} \
                ${R2_FILE_NAME} > ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam 2>/dev/null
            
            cat ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam | \
                samtools view -bhS - | \
                samtools view -h - | \
                pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 -c ${chromInfo} \
                -o ${pairsam_path}/${SAMPLE_PREFIX}.pairsam;
            
            rm ${pairsam_path}/${SAMPLE_PREFIX}.tmp.sam
        fi
    fi

    # stats
    echo "paired data `get_time` stats ..."
    echo " "
    pairtools stats ${pairsam_path}/${SAMPLE_PREFIX}.pairsam > ${stats_path}/${SAMPLE_PREFIX}.pairsam.stats.txt &
    pairtools sort --nproc ${THREAD} --memory 500G -o ${pairsam_path}/${SAMPLE_PREFIX}.sort.pairsam ${pairsam_path}/${SAMPLE_PREFIX}.pairsam
    pairtools stats ${pairsam_path}/${SAMPLE_PREFIX}.sort.pairsam > ${stats_path}/${SAMPLE_PREFIX}.sort.pairsam.stats.txt &
    pairtools dedup --mark-dups -o ${pairsam_path}/${SAMPLE_PREFIX}.sort.dedup.pairsam ${pairsam_path}/${SAMPLE_PREFIX}.sort.pairsam;
    pairtools stats ${pairsam_path}/${SAMPLE_PREFIX}.sort.dedup.pairsam > ${stats_path}/${SAMPLE_PREFIX}.sort.dedup.pairsam.stats.txt
    rm -fr ${pairsam_path}/${SAMPLE_PREFIX}.sort.pairsam;

    pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' -o ${pairsam_path}/${SAMPLE_PREFIX}.flt.pairsam ${pairsam_path}/${SAMPLE_PREFIX}.sort.dedup.pairsam;
    rm -fr ${pairsam_path}/${SAMPLE_PREFIX}.sort.dedup.pairsam;

    pairtools split --output-pairs ${pairsam_path}/${SAMPLE_PREFIX}.pairs ${pairsam_path}/${SAMPLE_PREFIX}.flt.pairsam;
    rm -fr ${pairsam_path}/${SAMPLE_PREFIX}.flt.pairsam;

    grep -v " " ${pairsam_path}/${SAMPLE_PREFIX}.pairs > ${pairsam_path}/${SAMPLE_PREFIX}.pairs2;
    rm -fr ${pairsam_path}/${SAMPLE_PREFIX}.pairs;
    mv ${pairsam_path}/${SAMPLE_PREFIX}.pairs2 ${pairs_path}/${SAMPLE_PREFIX}.pairs;

    bgzip ${pairs_path}/${SAMPLE_PREFIX}.pairs
    pairix -f ${pairs_path}/${SAMPLE_PREFIX}.pairs.gz

    # bin
    echo "cooler cload pairix `get_time` binning ..."
    echo " "
    for resolution in "${RES_ARRAY[@]}"; do
        if (( resolution >= 1000000 )); then
            suffix="${resolution:0:$((${#resolution}-6))}M"
        elif (( resolution >= 1000 )); then
            suffix="${resolution:0:$((${#resolution}-3))}k"
        else
            suffix="${resolution}"
        fi
        
        cooler cload pairix -p ${THREAD} ${chromInfo}:${resolution} ${pairs_path}/${SAMPLE_PREFIX}.pairs.gz ${cool_path}/${SAMPLE_PREFIX}.${suffix}.cool >/dev/null 2>&1
    done
done

multiqc $fastqc_path --outdir $fastqc_path
cool_files=($cool_path/*.cool)
total_files=${#cool_files[@]}

for cool_file in "${cool_files[@]}"; do
    # 获取分辨率信息
    resolution=$(cooler info "$cool_file" | grep "bin-size" | awk '{gsub(",", "", $2); print $2}')
    # 统计信息
    cooler dump -t pixels "$cool_file" 2>/dev/null | awk -v threshold=1000 -v file="$(basename "$cool_file")" -v res="$resolution" '
    {
        bin_sums[$1] += $3
        bins[$1] = 1
        total_pixels++
        total_contacts_all += $3
    }
    END {
        total_bins = 0
        high_contact_bins = 0
        total_contacts = 0
        max_contacts = 0
        min_contacts = 999999999
        
        for (bin in bins) {
            total_bins++
            contacts = bin_sums[bin]
            total_contacts += contacts
            if (contacts > max_contacts) max_contacts = contacts
            if (contacts < min_contacts) min_contacts = contacts
            if (contacts > 1000) high_contact_bins++
        }
        
        avg_contacts = total_contacts / total_bins
        ratio = high_contact_bins / total_bins
        avg_pixel_contacts = total_contacts_all / total_pixels
        
        # 输出结果
        printf "========== %s ==========\n",file
        printf "Resolutions: %sbp\n", res
        printf "Total bins: %d\n", total_bins
        printf "Total pixels: %d\n", total_pixels
        printf "Total contacts: %d\n", total_contacts_all
        printf "Bins contact > 1000: %d\n", high_contact_bins
        printf "Proportion of bin contact > 1000: %.2f%%\n", ratio * 100
        printf "Average bin contact: %.2f\n", avg_contacts
        printf "Average pixel contact: %.2f\n", avg_pixel_contacts
        printf "Max bin contact: %d\n", max_contacts
        printf "Min bin contact: %d\n", (min_contacts == 999999999 ? 0 : min_contacts)
        printf "\n"
    }' >> "${OUTPUT_DIRECTORY}/QC/Resolution_evaluation.txt"
done

echo "HiC preprocessing has been finished!!!"

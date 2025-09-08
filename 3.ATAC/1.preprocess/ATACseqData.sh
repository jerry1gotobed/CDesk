#!/usr/bin/env bash
config_json=${!#}
SCRIPT_DIR=$(dirname $(realpath $0))

Dofastqc=false

TRIM_GALORE=$(jq -r '.software.trim_galore' $config_json)
BOWTIE2=$(jq -r '.software.bowtie2' $config_json)
SAMTOOLS=$(jq -r '.software.samtools' $config_json)
bamToBed=$(jq -r '.software.bamToBed' $config_json)
bedSort=$(jq -r '.software.bedSort' $config_json)
bedClip=$(jq -r '.software.bedClip' $config_json)
bedGraphToBigWig=$(jq -r '.software.bedGraphToBigWig' $config_json)
bedToBigBed=$(jq -r '.software.bedToBigBed' $config_json)
genomeCoverageBed=$(jq -r '.software.genomeCoverageBed' $config_json)
macs14=$(jq -r '.software.macs14' $config_json)
FASTQC=$(jq -r '.software.fastqc' $config_json)
MULTIQC=$(jq -r '.software.multiqc' $config_json)

show_help() {
    echo "Usage: ATACseqData.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -s SPECIES		Specify the species"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (fastq.gz)"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file"
    echo "  -p INT		Specify number of threads (default is 8)"
    echo "  -l 1 or 2		1:Single sequencing, 2:Pair sequencing (default is 2)"
    echo "  -q		Do fastqc"
}

THREAD=8
LAYOUT=2

while getopts ":hs:i:o:p:l:q:g:" opt; do
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
        q )
            Dofastqc=true
            ;;
	g )
            Group=${OPTARG}
            ;;
        l )
            LAYOUT=${OPTARG}
            if [[ "$LAYOUT" != "1" && "$LAYOUT" != "2" ]]; then
                echo "Invalid value for -l: $LAYOUT. It must be 1 or 2." 1>&2
                show_help
                exit 1
            fi
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
if [ -z "$INPUT_DIRECTORY" ] || [ -z "$SPECIES" ] || [ -z "$OUTPUT_DIRECTORY" ]; then
    echo "Error: INPUT_DIRECTORY, OUTPUT_DIRECTORY and SPECIES must be provided."
    show_help
    exit 1
fi


CONFIG=$(cat $config_json | jq -r --arg species "$SPECIES" '.data[$species]')
mapping_index=$(echo "$CONFIG" | jq -r '.mapping_index')
species=$(echo "$CONFIG" | jq -r '.effective_genome_size')
chromInfo=$(echo "$CONFIG" | jq -r '.chromInfo')
promoter_file=$(echo "$CONFIG" | jq -r '.promoter_file')

# 提取文件名
if [ "$LAYOUT" = "1" ]; then
    file_list=$(ls "$INPUT_DIRECTORY"/*.{fastq,fq}.gz 2>/dev/null | sed 's/.fastq.gz//;s/.fq.gz//' | sort -u)
elif [ "$LAYOUT" = "2" ]; then
    file_list=$(ls "$INPUT_DIRECTORY"/*.{fastq,fq}.gz 2>/dev/null | sed 's/_1.fastq.gz//;s/_2.fastq.gz//;s/_1.fq.gz//;s/_2.fq.gz//' | sort -u)
fi

log_path="${OUTPUT_DIRECTORY}/Log"
bam_path="${OUTPUT_DIRECTORY}/Bam"
bw_path="${OUTPUT_DIRECTORY}/Bw"
peak_path="${OUTPUT_DIRECTORY}/Peak"
qc_path="${OUTPUT_DIRECTORY}/QC"
bed_path="${OUTPUT_DIRECTORY}/Bed"
signal_path="${OUTPUT_DIRECTORY}/Signal"
fastqc_path="${OUTPUT_DIRECTORY}/QC/fastqc"
temp_path="${OUTPUT_DIRECTORY}/Temp"

if [ ! -d "${OUTPUT_DIRECTORY}" ]; then
    mkdir -p ${OUTPUT_DIRECTORY}
fi
if [ ! -d "${log_path}" ]; then
    mkdir -p ${log_path}
fi
if [ ! -d "${bam_path}" ]; then
    mkdir -p ${bam_path}
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
if [ ! -d "${bed_path}" ]; then
    mkdir -p ${bed_path}
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


get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}

if [ -f "$Group" ]; then
    prefixes=$(tail -n +2 "$Group" | cut -d ',' -f 1 | sort -u)

    for prefix in $prefixes; do
        # 清理前缀中的特殊字符
        clean_prefix=$(echo "$prefix" | tr -d '\r')
        
        # 检查是否有匹配的文件
        if ! ls "${INPUT_DIRECTORY}/${clean_prefix}"* >/dev/null 2>&1; then
            echo "错误：找不到以 '$clean_prefix' 开头的文件"
            exit 1
        fi
    done

    # 检查第二列内容是否包含空格
    #second_columns=$(tail -n +2 "$Group" | cut -d ',' -f 2)
    second_columns=$(tail -n +2 "$Group" | cut -d ',' -f 2 | tr -d '\r')
    IFS=$'\n'  # 将 IFS 设置为换行符
    for second_column in $second_columns; do
        if [[ "$second_column" =~ [[:space:]] ]]; then
            echo "错误：第二列的内容包含空格"
            exit 1
        fi
    done
    unset IFS

fi

####################################################################################################
#################################          ANALYSIS STEPS           ################################
####################################################################################################
echo "检查是否有存在的bam文件"
find "$INPUT_DIRECTORY" -maxdepth 1 -name "*.bam" -exec sh -c '
    bam_file="$1"
    bam_dir="$2"
    sample=$(basename "$bam_file" .bam)
    target="${bam_dir}/${sample}.bam"

    # 如果目标文件不存在，创建硬链接
    if [ ! -e "$target" ]; then
        if ln "$bam_file" "$target" 2>/dev/null; then
            echo "Linked $bam_file -> $target"
        else
            cp "$bam_file" "$target"
            echo "Copy $bam_file -> $target"
        fi

        # 快速检查
        if ! "$3" quickcheck -q "$target"; then
            echo "BAM文件校验失败，可能需要重新处理: $target"
            rm "$target"
        fi
    fi
' sh {} "$bam_path" "$SAMTOOLS" \;

echo "---------------------------------------处理fq文件-----------------------------------------------------"
counter=0

for FILE in $file_list; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename ${FILE})
    BAM_FILE="${bam_path}/${SAMPLE_PREFIX}.bam"
    echo "-----------------------------------------------------------------------------------------------------"
    echo "----------------------------第 ${counter} 个样本，样本名为 ${SAMPLE_PREFIX}--------------------------"
    echo "-----------------------------------------------------------------------------------------------------"

    if [ -f "$BAM_FILE" ]; then
        echo "检测到现有BAM文件，跳过比对"
        SKIP_MAPPING=true
    else
        echo "未检测到现有BAM文件，进行比对"
        SKIP_MAPPING=false
    fi

    FILE_NAME=""
    R1_FILE_NAME=""
    R2_FILE_NAME=""

    if [ "${LAYOUT}" == "1" ]; then
        # 检查是否存在对应的文件
        if [ -e "${FILE}.fastq.gz" ] || [ -e "${FILE}.fq.gz" ]; then
            if [ -e "${FILE}.fastq.gz" ]; then
                FILE_NAME="${FILE}.fastq.gz"
            else
                FILE_NAME="${FILE}.fq.gz"
            fi
        else
            echo "错误：未找到样本 ${SAMPLE_REFIX} 的fq文件。跳过此样本。"
            continue
        fi
    elif [ "${LAYOUT}" == "2" ]; then
        # 检查是否存在对应的文件
        if [ -e "${FILE}_1.fastq.gz" ] && [ -e "${FILE}_2.fastq.gz" ]; then
            R1_FILE_NAME="${FILE}_1.fastq.gz"
            R2_FILE_NAME="${FILE}_2.fastq.gz"
            FILE_TYPE="fastq"
        elif [ -e "${FILE}_1.fq.gz" ] && [ -e "${FILE}_2.fq.gz" ]; then
            R1_FILE_NAME="${FILE}_1.fq.gz"
            R2_FILE_NAME="${FILE}_2.fq.gz"
            FILE_TYPE="fq"
        else
            echo "错误：未找到样本 ${SAMPLE_REFIX} 的2个fq文件。跳过此样本。"
            continue
        fi
    fi

    if $Dofastqc; then
        echo "fastqc `get_time` fastqc ..."
        if [[ "${FILE_TYPE}" == "fastq" ]]; then
            $FASTQC -o ${qc_path}/fastqc ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fastq.gz
        elif [[ "${FILE_TYPE}" == "fq" ]]; then
            $FASTQC -o ${qc_path}/fastqc ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fq.gz
        fi
    fi
    
    if [ "$SKIP_MAPPING" = false ]; then
        echo " "
        echo "--------------------------------------------INITIALIZING----------------------------------------------" 
        echo "ATAC-seq data analysis pipeline is now running..."
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
        echo "PromoterInfo ---------- ${promoter_file}"
        echo " "

        echo " "
        echo "-----------------------------------GENERATING DIRECTORY-------------------------------------" 
        echo "Please wait, the directory for storing the results is being generated..."
        echo "Store LOG records for all samples ---------- ${log_path}" 
        echo "Store BAM files for all samples ---------- ${bam_path}"
        echo "Store BW results for all samples ---------- ${bw_path}"
        echo "Store PEAK data for all samples ---------- ${peak_path}"
        echo "Store QC result for all samples ---------- ${qc_path}"
        echo "Store BED result for all samples ---------- ${bed_path}"
        echo "Store PromoterSignal result for all samples ---------- ${signal_path}"
        echo "Store temporary result for all samples ---------- ${temp_path}"
        echo " "

        echo " "
        echo "---------------------------------------START ANALYZING-----------------------------------------" 
        echo "Preparation work completed, start analysis..."

        # 1.mapping ----------------------------------------------------------------------------------------
        echo " "
        echo "1.mapping..."
        echo " "

        # trim
        if [ "${LAYOUT}" == "1" ]; then
            echo "single data `get_time` trim_galore ..."
            echo " "
            $TRIM_GALORE --gzip -j ${THREAD} ${FILE_NAME} --trim-n -o ${temp_path}/ --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1
        elif [ "${LAYOUT}" == "2" ]; then
            echo "paired data `get_time` trim_galore ..."
            echo " "
            $TRIM_GALORE --gzip -j ${THREAD} --paired ${R1_FILE_NAME} ${R2_FILE_NAME} --trim-n -o ${temp_path} --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1
        fi

        # map
        if [ "${LAYOUT}" == "1" ]; then
            echo "single data `get_time` bowtie2 ..."
            echo " "
            $BOWTIE2 -q -p ${THREAD} -x ${mapping_index} -U ${temp_path}/${SAMPLE_PREFIX}_trimmed.fq.gz -S ${bam_path}/${SAMPLE_PREFIX}.align.sam > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1
            rm -rf ${temp_path}/${SAMPLE_PREFIX}_trimmed.fq.gz
        elif [ "${LAYOUT}" == "2" ]; then
            echo "paired data `get_time` bowtie2 ..."
            echo " "
            $BOWTIE2 -q -p ${THREAD} -x ${mapping_index} -1 ${temp_path}/${SAMPLE_PREFIX}_val_1.fq.gz -2 ${temp_path}/${SAMPLE_PREFIX}_val_2.fq.gz -S ${bam_path}/${SAMPLE_PREFIX}.align.sam > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1
            rm -rf ${temp_path}/${SAMPLE_PREFIX}_val_1.fq.gz ${temp_path}/${SAMPLE_PREFIX}_val_2.fq.gz
        fi

        # sam2bam
        echo "sam2bam `get_time` samtools view ..."
        echo " "
        $SAMTOOLS view -@ ${THREAD} -bS ${bam_path}/${SAMPLE_PREFIX}.align.sam > ${bam_path}/${SAMPLE_PREFIX}.bam
    else
        echo "使用现有BAM文件: $BAM_FILE"
    fi

done

if $Dofastqc; then
    $MULTIQC ${qc_path}/fastqc --outdir ${qc_path}/fastqc
fi

echo "---------------------------------------处理bam文件----------------------------------------------------"
file_list=$(ls "$bam_path"/*.bam 2>/dev/null | sort -u)
for FILE in $file_list; do
    SAMPLE_PREFIX=$(basename "$FILE" | cut -d '.' -f 1)
    echo "处理 ${FILE} 文件"
    echo " "
    # sort bam
    echo "bam2sortbam `get_time` samtools sort ..."
    echo " "
    if ! $SAMTOOLS view -H "${bam_path}/${SAMPLE_PREFIX}.bam" | grep -q "SO:coordinate"; then
        echo "重新排序BAM文件..."
        tmp_sorted=$(mktemp)
        $SAMTOOLS sort -@ $THREAD "${bam_path}/${SAMPLE_PREFIX}.bam" -o "$tmp_sorted"
        mv "$tmp_sorted" "${bam_path}/${SAMPLE_PREFIX}.bam"
    else
        echo "BAM文件已排序"
    fi
done

if [ -f "$Group" ]; then
    echo "Check replicates ..."
    python ${SCRIPT_DIR}/Merge_replicate.py $Group $bam_path $SAMTOOLS  
fi

file_list=$(ls "$bam_path"/*.bam 2>/dev/null | sort -u)
for FILE in $file_list; do
    SAMPLE_PREFIX=$(basename "$FILE" | cut -d '.' -f 1)
    echo "处理 ${FILE} 文件"
    # chromosome distribution
    if [ "${LAYOUT}" == "1" ]; then
        echo "single data `get_time` chromosome distribution ..."
        echo " "
        $SAMTOOLS view -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam | awk '/^[^#]/ {print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${qc_path}/${SAMPLE_PREFIX}.ChromosomeDistribution.txt
        # 去除重复的reads，去除质量低于30的reads
        $SAMTOOLS view -@ ${THREAD} -F 0x400 ${bam_path}/${SAMPLE_PREFIX}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3 !~ /_/) print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${qc_path}/${SAMPLE_PREFIX}.Filtered.ChromosomeDistribution.txt
    elif [ "${LAYOUT}" == "2" ]; then
        echo "paired data `get_time` chromosome distribution ..."
        echo " "
        $SAMTOOLS view -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam | awk '/^[^#]/ {print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${qc_path}/${SAMPLE_PREFIX}.ChromosomeDistribution.txt
        # 去除未成功匹配的reads，去除质量低于30的reads
        $SAMTOOLS view -@ ${THREAD} -f 0x2 ${bam_path}/${SAMPLE_PREFIX}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3 !~ /_/) print $3}' | sort | uniq -c | sort -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${qc_path}/${SAMPLE_PREFIX}.Filtered.ChromosomeDistribution.txt
    fi

    # filtering
    if [ "${LAYOUT}" == "1" ]; then
        echo "single data `get_time` filtering ..."
        echo " "
        $SAMTOOLS view -@ ${THREAD} -H ${bam_path}/${SAMPLE_PREFIX}.bam > ${bam_path}/${SAMPLE_PREFIX}.Filtered.sam
        $SAMTOOLS view -@ ${THREAD} -F 0x400 ${bam_path}/${SAMPLE_PREFIX}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3!="chrM" && $3!="chrEBV" && $3 !~ /_/) print forward"\n"$0}' >> ${bam_path}/${SAMPLE_PREFIX}.Filtered.sam
	$SAMTOOLS view -@ ${THREAD} -bS ${bam_path}/${SAMPLE_PREFIX}.Filtered.sam > ${bam_path}/${SAMPLE_PREFIX}.Filtered.bam
        rm ${bam_path}/${SAMPLE_PREFIX}.Filtered.sam
    elif [ "${LAYOUT}" == "2" ]; then
        echo "paired data `get_time` filtering ..."
        echo " "
        $SAMTOOLS view -@ ${THREAD} -H ${bam_path}/${SAMPLE_PREFIX}.bam > ${bam_path}/${SAMPLE_PREFIX}.Filtered.sam
        $SAMTOOLS view -@ ${THREAD} -f 0x2 ${bam_path}/${SAMPLE_PREFIX}.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $3!="chrM" && $3!="chrEBV" && $3 !~ /_/) print forward"\n"$0}' >> ${bam_path}/${SAMPLE_PREFIX}.Filtered.sam
        $SAMTOOLS view -@ ${THREAD} -bS ${bam_path}/${SAMPLE_PREFIX}.Filtered.sam > ${bam_path}/${SAMPLE_PREFIX}.Filtered.bam
	rm ${bam_path}/${SAMPLE_PREFIX}.Filtered.sam
    fi

    # 2.BW ----------------------------------------------------------------------------------------
    echo " "
    echo "2.generating BigWig..."
    echo " "

    if [ "${LAYOUT}" == "1" ]; then
        echo "single data `get_time` BigWig ..."
        echo " "
        $bamToBed -i ${bam_path}/${SAMPLE_PREFIX}.Filtered.bam > ${temp_path}/${SAMPLE_PREFIX}.tmp.bed
        $bedSort ${temp_path}/${SAMPLE_PREFIX}.tmp.bed ${bed_path}/${SAMPLE_PREFIX}.bed
        n=`wc -l ${bed_path}/${SAMPLE_PREFIX}.bed | cut -f 1 -d " "`
        c=`bc -l <<< "1000000 / $n"` 
        $genomeCoverageBed -bga -scale ${c} -i ${bed_path}/${SAMPLE_PREFIX}.bed -g ${chromInfo} > ${temp_path}/${SAMPLE_PREFIX}.bdg 
        $bedClip ${temp_path}/${SAMPLE_PREFIX}.bdg ${chromInfo} ${temp_path}/${SAMPLE_PREFIX}.Clip.bdg
        $bedSort ${temp_path}/${SAMPLE_PREFIX}.Clip.bdg ${temp_path}/${SAMPLE_PREFIX}.bdg.temp
        $bedGraphToBigWig ${temp_path}/${SAMPLE_PREFIX}.bdg.temp ${chromInfo} ${bw_path}/${SAMPLE_PREFIX}.bw
        rm ${temp_path}/${SAMPLE_PREFIX}.bdg
        rm ${temp_path}/${SAMPLE_PREFIX}.bdg.temp
        rm ${temp_path}/${SAMPLE_PREFIX}.Clip.bdg
        rm ${temp_path}/${SAMPLE_PREFIX}.tmp.bed
    elif [ "${LAYOUT}" == "2" ]; then
        echo "paired data `get_time` BigWig ..."
        echo " "
        $bamToBed -i ${bam_path}/${SAMPLE_PREFIX}.Filtered.bam > ${bed_path}/${SAMPLE_PREFIX}.RawReads.bed
        $SAMTOOLS sort -n -@ 20 ${bam_path}/${SAMPLE_PREFIX}.Filtered.bam > ${bam_path}/${SAMPLE_PREFIX}.FilteredSorted.bam
        rm ${bam_path}/${SAMPLE_PREFIX}.Filtered.bam
        $bamToBed -bedpe -i ${bam_path}/${SAMPLE_PREFIX}.FilteredSorted.bam 2>/dev/null | awk '$1==$4 && $2>0 && $5>0 {print $1"\t"$2"\t"$6}' | sort -k1,1 -k2,2n | uniq | grep -v WARNING > ${bed_path}/${SAMPLE_PREFIX}.Fragments.bed
        awk '{print $3-$2}' ${bed_path}/${SAMPLE_PREFIX}.Fragments.bed | sort | uniq -c | sort -k2,2g | awk 'BEGIN{print "fragment_length\tnumber"} {print $2"\t"$1}' > ${qc_path}/${SAMPLE_PREFIX}.FragmentsLength.txt
        Rscript ${SCRIPT_DIR}/FragmentLengthDistribution.R ${qc_path}/${SAMPLE_PREFIX}.FragmentsLength.txt ${qc_path}
        $bedSort ${bed_path}/${SAMPLE_PREFIX}.Fragments.bed ${bed_path}/${SAMPLE_PREFIX}.bed
        n=`wc -l ${bed_path}/${SAMPLE_PREFIX}.bed | cut -f 1 -d " "`
        c=`bc -l <<< "1000000 / $n"`
        $genomeCoverageBed -bga -scale $c -i ${bed_path}/${SAMPLE_PREFIX}.bed -g ${chromInfo} > ${temp_path}/${SAMPLE_PREFIX}.bdg
        $bedClip ${temp_path}/${SAMPLE_PREFIX}.bdg ${chromInfo} ${temp_path}/${SAMPLE_PREFIX}.Clip.bdg
        $bedSort ${temp_path}/${SAMPLE_PREFIX}.Clip.bdg ${temp_path}/${SAMPLE_PREFIX}.bdg.temp
        $bedGraphToBigWig ${temp_path}/${SAMPLE_PREFIX}.bdg.temp ${chromInfo} ${bw_path}/${SAMPLE_PREFIX}.bw
        rm ${temp_path}/${SAMPLE_PREFIX}.bdg
        rm ${temp_path}/${SAMPLE_PREFIX}.bdg.temp
        rm ${temp_path}/${SAMPLE_PREFIX}.Clip.bdg
    fi

    echo "generate bigbed `get_time` bedToBigBed ..."
    echo " "
    $bedToBigBed ${bed_path}/${SAMPLE_PREFIX}.bed ${chromInfo} ${bed_path}/${SAMPLE_PREFIX}.bb

    # 3.peak calling ----------------------------------------------------------------------------------------
    echo " "
    echo "3.peak calling..."
    echo " "

    $macs14 -t ${bed_path}/${SAMPLE_PREFIX}.bed -n ${peak_path}/${SAMPLE_PREFIX}_e3 -f BED -g ${species} --keep-dup all --nomodel --shiftsize 25 -p 1e-3 > ${log_path}/${SAMPLE_PREFIX}_e3.PeakCalling.log 2>&1
    $macs14 -t ${bed_path}/${SAMPLE_PREFIX}.bed -n ${peak_path}/${SAMPLE_PREFIX}_e5 -f BED -g ${species} --keep-dup all --nomodel --shiftsize 25 -p 1e-5 > ${log_path}/${SAMPLE_PREFIX}_e5.PeakCalling.log 2>&1
    $macs14 -t ${bed_path}/${SAMPLE_PREFIX}.bed -n ${peak_path}/${SAMPLE_PREFIX}_e7 -f BED -g ${species} --keep-dup all --nomodel --shiftsize 25 -p 1e-7 > ${log_path}/${SAMPLE_PREFIX}_e7.PeakCalling.log 2>&1

    # 4.Promoter signal calculation ----------------------------------------------------------------------------------------
    echo " "
    echo "4.Promoter signal calculation..."
    echo " "
    python ${SCRIPT_DIR}/GetBigwigSignal.py ${promoter_file} ${bw_path}/${SAMPLE_PREFIX}.bw ${signal_path}/${SAMPLE_PREFIX}.PromoterSignal.csv
    
done

rm -rf ${temp_path}
echo "ATACseq preprocessing has been finished!!!"

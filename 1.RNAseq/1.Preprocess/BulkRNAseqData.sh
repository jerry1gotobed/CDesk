#!/bin/bash

# Description : Preprocessing+Mapping+QC
# Author      : WEI SHI, Lin Zejie
# Version     : 1.0
# Time        : 2024-10-19 22:40:00

config_json=${!#}
SCRIPT_DIR=$(dirname $(realpath $0))

show_help() {
    echo "Usage: BulkRNAseqData.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help		Show this help message"
    echo "  -s SPECIES		Specify the species (mm9, mm10)"
    echo "  -i INPUT_DIRECTORY	Specify the absolute directory of input file (fastq.gz), the ending does not require a '/'"
    echo "  -o OUTPUT_DIRECTORY	Specify the absolute directory of output file, the ending does not require a '/'"
    echo "  -p INT		Specify number of threads (default is 8)"
    echo "  -l 1 or 2		1:Single sequencing, 2:Pair sequencing (default is 2)"
}

# 检查是否安装了 `jq`（用于解析 JSON）
if ! command -v jq &> /dev/null; then
    echo "Error: 'jq' is required but not installed. Please install it first."
    exit 1
fi

# 检查软件是否存在
HISAT2=$(jq -r '.software.hisat2' $config_json)
STRINGTIE=$(jq -r '.software.stringtie' $config_json)
SAMTOOLS=$(jq -r '.software.samtools' $config_json)
TRIM_GALORE=$(jq -r '.software.trim_galore' $config_json)
FASTQC=$(jq -r '.software.fastqc' $config_json)
MULTIQC=$(jq -r '.software.multiqc' $config_json)

for tool in "$HISAT2" "$STRINGTIE" "$SAMTOOLS" "$TRIM_GALORE" "$FASTQC" "$MULTIQC"; do
    if [ ! -x "$tool" ]; then
        echo "Error: Tool not found or not executable: $tool"
        exit 1
    fi
done

THREAD=8
LAYOUT=2
ALLOWED_SPECIES=("mm10" "hg38" "rn7" "susScr11" "galGal6")

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
	    if [[ ! " ${ALLOWED_SPECIES[@]} " =~ " ${SPECIES} " ]]; then
                echo "Error: Invalid species '${SPECIES}'. Allowed species are: ${ALLOWED_SPECIES[*]}"
                exit 1
            fi
            ;;
        o)
            OUTPUT_DIRECTORY=${OPTARG}
            ;;
        p )
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
    show_helpTRIM_GALORE=$
    exit 1
fi

CONFIG=$(cat $config_json | jq -r --arg species "$SPECIES" '.data[$species]')
if [ -z "$CONFIG" ]; then
    echo "Error: No configuration found for species '$SPECIES' in data.json"
    exit 1
fi

# 分配参考数据
mapping_index=$(echo "$CONFIG" | jq -r '.mapping_index')
mapping_gtf=$(echo "$CONFIG" | jq -r '.refseq_gtf')
gfold_gtf=$(echo "$CONFIG" | jq -r '.refseq_gtf') 
RSeQC_bed=$(echo "$CONFIG" | jq -r '.refseq_bed')
# 检查变量是否为空
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
# 检查文件是否存在
# 索引目录进行检查
if ! ls "${mapping_index}"*.ht2 1> /dev/null 2>&1; then
  echo "HISAT2 index files are missing. Please check the index path."
fi
# 检查剩余文件是否存在
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


# 提取文件名
if [ "$LAYOUT" = "1" ]; then
    file_elist=$(ls "$INPUT_DIRECTORY"/*.fastq.gz 2>/dev/null | sed 's/.fastq.gz//' | sort -u)
    fqlist=$(ls "$INPUT_DIRECTORY"/*.fq.gz 2>/dev/null | sed 's/.fq.gz//' | sort -u)
    file_list=$(echo "$file_elist $fqlist" | tr ' ' '\n' | sort -u)
elif [ "$LAYOUT" = "2" ]; then
    file_elist=$(ls "$INPUT_DIRECTORY"/*.fastq.gz 2>/dev/null | sed 's/_1.fastq.gz//;s/_2.fastq.gz//' | sort -u)
    fqlist=$(ls "$INPUT_DIRECTORY"/*.fq.gz 2>/dev/null | sed 's/_1.fq.gz//;s/_2.fq.gz//' | sort -u)
    file_list=$(echo "$file_elist $fqlist" | tr ' ' '\n' | sort -u)
fi


log_path="${OUTPUT_DIRECTORY}/Log"
bam_path="${OUTPUT_DIRECTORY}/Bam"
QC_path="${OUTPUT_DIRECTORY}/QC"
Expression_path="${OUTPUT_DIRECTORY}/Expression"


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

get_time(){
    printf "%-19s" "`date +\"%Y-%m-%d %H:%M:%S\"`"
}


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
        ln "$bam_file" "$target"
	echo "检查到Bam文件"
        echo "Linked input BAM: $bam_file -> $target"
        
        # 快速检查
        if ! samtools quickcheck -q "$target"; then
            echo "BAM文件校验失败，可能需要重新处理: $target"
            rm "$target"
        fi
    fi
' sh {} "$bam_path" \;

echo "---------------------------------------处理fq文件-----------------------------------------------------"
counter=0

for FILE in $file_list; do
    counter=$((counter + 1))
    SAMPLE_PREFIX=$(basename ${FILE})
    BAM_FILE="${bam_path}/${SAMPLE_PREFIX}.bam"
    echo "----------------------------第 ${counter} 个fq样本，样本名为 ${SAMPLE_PREFIX}--------------------------"
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
        # 检查是否存在对应的SE文件
        if [ -e "${FILE}.fastq.gz" ] || [ -e "${FILE}.fq.gz" ]; then
            if [ -e "${FILE}.fastq.gz" ]; then
                FILE_NAME="${FILE}.fastq.gz"
            else
                FILE_NAME="${FILE}.fq.gz"
            fi
        else
            echo "错误：未找到样本 ${SAMPLE_REFIX} 的SE文件。跳过此样本。"
            continue
        fi
    elif [ "${LAYOUT}" == "2" ]; then
        # 检查是否存在对应的PE文件
        if [ -e "${FILE}_1.fastq.gz" ] && [ -e "${FILE}_2.fastq.gz" ]; then
            R1_FILE_NAME="${FILE}_1.fastq.gz"
            R2_FILE_NAME="${FILE}_2.fastq.gz"
	    FILE_TYPE="fastq"
        elif [ -e "${FILE}_1.fq.gz" ] && [ -e "${FILE}_2.fq.gz" ]; then
            R1_FILE_NAME="${FILE}_1.fq.gz"
            R2_FILE_NAME="${FILE}_2.fq.gz"
	    FILE_TYPE="fq"
        else
            echo "错误：未找到样本 ${SAMPLE_REFIX} 的PE文件。跳过此样本。"
            continue
        fi
    fi

    if [ "$SKIP_MAPPING" = false ]; then
        echo " "
        echo "--------------------------------------------INITIALIZING----------------------------------------------" 
        echo "RNA-seq data analysis pipeline is now running..."
        echo "Number of threads ---------- ${THREAD}" 
        echo "Directory of data ---------- ${INPUT_DIRECTORY}"
        echo "Directory of result ---------- ${OUTPUT_DIRECTORY}"
        echo "Name of sample ---------- ${SAMPLE_PREFIX}"
        echo "File of data ---------- ${FILE_NAME}"
        echo "File of R1 data ---------- ${R1_FILE_NAME}"
        echo "File of R2 data ---------- ${R2_FILE_NAME}"
        echo "Method of sequencing ---------- ${LAYOUT}"
        echo "Mapping index ---------- ${mapping_index}"
        echo "Mapping gtf ---------- ${mapping_gtf}"
        echo "GFOLD gtf ---------- ${gfold_gtf}"
        echo "RSeQC bed ---------- ${RSeQC_bed}"
        echo "----------------------------------INITIALIZATION FINISHED-------------------------------------"
        echo " "

        echo " "
        echo "-----------------------------------GENERATING DIRECTORY-------------------------------------" 
        echo "Please wait, the directory for storing the results is being generated..."


        echo "Save LOG records for all samples ---------- ${log_path}" 
        echo "Store BAM files for all samples ---------- ${bam_path}"
        echo "Store quality control results, including fastqc results, mapping results, and RSeQC results ---------- ${QC_path}"
        echo "Store gene expression data, including count and fpkm versions ---------- ${Expression_path}"
        echo "---------------------------DIRECTORY GENERATION FINISHED-----------------------------"
        echo " "

        echo " "
        echo "---------------------------------------START ANALYZING-----------------------------------------" 
        echo "Preparation work completed, start analysis..."

        # 1.mapping ----------------------------------------------------------------------------------------
        echo " "
        echo "1.mapping..."
        echo " "
        if [ "${LAYOUT}" == "1" ]; then
            echo "67 `get_time` trim_galore ..."
            echo " "
        $TRIM_GALORE --gzip -j ${THREAD} ${FILE_NAME} --trim-n -o ${OUTPUT_DIRECTORY} --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1
        elif [ "${LAYOUT}" == "2" ]; then
            echo "86 `get_time` trim_galore ..."
            echo " "
        $TRIM_GALORE --gzip -j ${THREAD} --paired ${R1_FILE_NAME} ${R2_FILE_NAME} --trim-n -o ${OUTPUT_DIRECTORY} --no_report_file --basename ${SAMPLE_PREFIX} > ${log_path}/${SAMPLE_PREFIX}.TrimGalore.log 2>&1
        fi


        if [ "${LAYOUT}" == "1" ]; then
            echo "70 `get_time` hisat2 ..."
            echo " "
        $HISAT2 -p ${THREAD} --dta -x ${mapping_index} -U ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_trimmed.fq.gz -S ${bam_path}/${SAMPLE_PREFIX}.align.sam --summary-file ${QC_path}/${SAMPLE_PREFIX}_mapping_summary.txt > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1 
            rm -rf ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_trimmed.fq.gz
        elif [ "${LAYOUT}" == "2" ]; then
            echo "91 `get_time` hisat2 ..."
            echo " "
        $HISAT2 -p ${THREAD} --dta -x ${mapping_index} -1 ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_val_1.fq.gz -2 ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_val_2.fq.gz -S ${bam_path}/${SAMPLE_PREFIX}.align.sam --summary-file ${QC_path}/${SAMPLE_PREFIX}_mapping_summary.txt > ${log_path}/${SAMPLE_PREFIX}.Mapping.log 2>&1
            rm -rf ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_val_1.fq.gz ${OUTPUT_DIRECTORY}/${SAMPLE_PREFIX}_val_2.fq.gz
        fi

        echo "74 & 95 `get_time` samtools view ..."
        echo " "
        $SAMTOOLS view -@ ${THREAD} -bS ${bam_path}/${SAMPLE_PREFIX}.align.sam > ${bam_path}/${SAMPLE_PREFIX}.bam

        # echo "77 & 98 `get_time` samtools sort ..."
        # echo " "
        # $SAMTOOLS sort -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam > ${bam_path}/${SAMPLE_PREFIX}.Sorted.bam
        # rm -rf ${bam_path}/${SAMPLE_PREFIX}.align.sam ${bam_path}/${SAMPLE_PREFIX}.bam
        # mv ${bam_path}/${SAMPLE_PREFIX}.Sorted.bam ${bam_path}/${SAMPLE_PREFIX}.bam

        # echo "112 `get_time` samtools index ..."
        # echo " "
        # $SAMTOOLS index -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam
    else
        echo "使用现有BAM文件: $BAM_FILE"
    fi

    echo " "
    echo "120 `get_time` fastqc ..."
    echo " "
    if [[ "${FILE_TYPE}" == "fastq" ]]; then
        $FASTQC -o ${QC_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fastq.gz
    elif [[ "${FILE_TYPE}" == "fq" ]]; then
        $FASTQC -o ${QC_path} ${INPUT_DIRECTORY}/${SAMPLE_PREFIX}*.fq.gz
    fi
done
$MULTIQC $QC_path --outdir $QC_path

echo "---------------------------------------处理bam文件----------------------------------------------------"
file_list=$(ls "$bam_path"/*.bam 2>/dev/null | sort -u)
for FILE in $file_list; do
    SAMPLE_PREFIX=$(basename "$FILE" | cut -d '.' -f 1)
    echo "处理 ${FILE} 文件"
    
    # 检查排序状态
    if ! samtools view -H "${bam_path}/${SAMPLE_PREFIX}.bam" | grep -q "SO:coordinate"; then
        echo "重新排序BAM文件..."
        tmp_sorted=$(mktemp)
        $SAMTOOLS sort -@ $THREAD "${bam_path}/${SAMPLE_PREFIX}.bam" -o "$tmp_sorted"
        mv "$tmp_sorted" "${bam_path}/${SAMPLE_PREFIX}.bam"
    else
        echo "BAM文件已排序"
    fi

    # 检查是否需要建立索引
    if [ ! -f "${bam_path}/${SAMPLE_PREFIX}.bam.bai" ]; then
        echo "创建BAM索引..."
        $SAMTOOLS index -@ $THREAD "${bam_path}/${SAMPLE_PREFIX}.bam"
    fi

    # 2.Calculate gene expression levels --------------------------------------------------------------------
    echo " "
    echo "2.Calculating gene expression levels ..."
    echo " "
    echo "116 `get_time` gfold count ..."
    echo " "
    $SAMTOOLS view -@ ${THREAD} ${bam_path}/${SAMPLE_PREFIX}.bam | gfold count -ann ${gfold_gtf} -tag stdin -o ${Expression_path}/${SAMPLE_PREFIX}.read_cnt


    if [ ! -d "${Expression_path}/Stringtie" ]; then
        mkdir -p ${Expression_path}/Stringtie
    fi

    
    echo "160 `get_time` stringtie count ..."
    echo " "
    nohup $STRINGTIE ${bam_path}/${SAMPLE_PREFIX}.bam -p ${THREAD} -G ${mapping_gtf} -A ${Expression_path}/Stringtie/${SAMPLE_PREFIX}_gene_abund.tab > ${log_path}/${SAMPLE_PREFIX}.Stringtie.log 2>&1


    # 3.QC --------------------------------------------------------------------
    echo "124 `get_time` RSeQC ..."
    echo " "
    ${SCRIPT_DIR}/read_distribution.py -i ${bam_path}/${SAMPLE_PREFIX}.bam -r ${RSeQC_bed} > ${QC_path}/${SAMPLE_PREFIX}_dirstribution.txt

    echo "---------------------------ANALYZING FINISHED-----------------------------"
    echo " "

done

echo " "
echo "begin merge fpkm expression..."
echo " "
# -------------------------------------------合并fpkm----------------------------------------------
# 提取第一个.read_cnt文件的第一列并追加到CSV文件中
awk '{print $1}' $(ls ${Expression_path}/*.read_cnt | head -n 1) >> ${Expression_path}/output1.csv

# 遍历所有的.read_cnt文件
for file in ${Expression_path}/*.read_cnt; do
    paste -d ',' ${Expression_path}/output1.csv <(awk '{print $5}' "$file") > ${Expression_path}/temp1.csv
    mv ${Expression_path}/temp1.csv ${Expression_path}/output1.csv
done

# 加表头
echo "gene_name,$(ls ${Expression_path}/*.read_cnt | xargs -n 1 basename | sed 's/.read_cnt//;s/^/,/' | tr -d '\n' | sed 's/^,//')" > ${Expression_path}/temp1.csv

# 将原始数据追加到临时文件中
cat ${Expression_path}/output1.csv >> ${Expression_path}/temp1.csv
# 用临时文件替换原始文件
mv ${Expression_path}/temp1.csv ${Expression_path}/merged_fpkm.csv

rm -rf ${Expression_path}/output1.csv
echo " "
echo "fpkm expression has been merged..."
echo " "

echo " "
echo "begin merge counts expression..."
echo " "
# -------------------------------------------合并count----------------------------------------------
# 提取第一个.read_cnt文件的第一列并追加到CSV文件中
awk '{print $1}' $(ls ${Expression_path}/*.read_cnt | head -n 1) >> ${Expression_path}/output1.csv

# 遍历所有的.read_cnt文件
for file in ${Expression_path}/*.read_cnt; do
    paste -d ',' ${Expression_path}/output1.csv <(awk '{print $3}' "$file") > ${Expression_path}/temp1.csv
    mv ${Expression_path}/temp1.csv ${Expression_path}/output1.csv
done

echo "gene_name,$(ls ${Expression_path}/*.read_cnt | xargs -n 1 basename | sed 's/.read_cnt//;s/^/,/' | tr -d '\n' | sed 's/^,//')" > ${Expression_path}/temp1.csv

# 将原始数据追加到临时文件中
cat ${Expression_path}/output1.csv >> ${Expression_path}/temp1.csv
# 用临时文件替换原始文件
mv ${Expression_path}/temp1.csv ${Expression_path}/merged_count.csv

rm -rf ${Expression_path}/output1.csv
echo " "
echo "counts expression has been merged..."
echo "运行完成"


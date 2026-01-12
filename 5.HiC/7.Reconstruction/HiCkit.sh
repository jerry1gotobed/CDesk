#!/usr/bin/env bash

ulimit -c unlimited

# Parameters initialize
VCF=""
MAT=""
PAT=""
R1=""
R2=""
SPECIES=""
OUT_DIR=""
OUT=""
THREAD=""

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -vcf)
      VCF="$2"; shift; shift ;;
    -mat)
      MAT="$2"; shift; shift ;;
    -pat)
      PAT="$2"; shift; shift ;;
    -r1)
      R1="$2"; shift; shift ;;
    -r2)
      R2="$2"; shift; shift ;;
    -species)
      SPECIES="$2"; shift; shift ;;
    -out)
      OUT="$2"; shift; shift ;;
    -t)
      THREAD="$2"; shift; shift ;;
    -o)
      OUT_DIR="$2"; shift; shift ;;
    *)
      echo "Unknown option $1"
      echo "Usage: $0 (-vcf <SNP.vcf.gz> -mat <MaternalSample> -pat <PaternalSample>) -r1 <HiC_R1.fastq.gz> -r2 <HiC_R2.fastq.gz> -species <Species> -out <OutputPrefix> -o <OutputDir> -t <Thread>"
      exit 1 ;;
  esac
done

if [ -z "$R1" ] || [ -z "$R2" ] || [ -z "$SPECIES" ] || [ -z "$OUT" ] || [ -z "$OUT_DIR" ] || [ -z "$THREAD" ]; then
  echo "Missing one or more required parameters!"
  echo "Usage: $0 (-vcf <SNP.vcf.gz> -mat <MaternalSample> -pat <PaternalSample>) -r1 <HiC_R1.fastq.gz> -r2 <HiC_R2.fastq.gz> -species <Species> -out <OutputPrefix> -o <OutputDir> -t <Thread>"
  exit 1
fi

if [ -n "$VCF" ] && { [ -z "$MAT" ] || [ -z "$PAT" ]; }; then
  echo "Error: When VCF is provided, both MAT and PAT must be specified!"
  exit 1
fi

tools=("jq" "seqkit" "bcftools" "bwa" "python2.7" "bgzip" "bash")
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

if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p "${OUT_DIR}"
fi

CONFIG=$(cat $CDesk_config | jq -r --arg species "$SPECIES" '.data[$species]')
fasta=$(echo "$CONFIG" | jq -r '.fasta')

DIP_C_DIR=$(jq -r '.software.dipc' $CDesk_config)
SCRIPTS_DIR="${DIP_C_DIR}/scripts"
OUT="${OUT_DIR}/${OUT}"
hickit=$(jq -r '.software.hickit' $CDesk_config)
hickit_js=$(jq -r '.software.hickit_js' $CDesk_config)
if [ ! -d "$DIP_C_DIR" ]; then
  echo "Error: DIP-C directory $DIP_C_DIR does not exist!"
  exit 1
fi
if [ ! -f "$hickit" ]; then
  echo "Error: hickit executable $hickit does not exist!"
  exit 1
fi
if [ ! -f "$hickit_js" ]; then
  echo "Error: hickit.js executable $hickit_js does not exist!"
  exit 1
fi

echo "[Step 1] Extracting and filtering informative SNPs..."
if [ -n "$VCF" ]; then
  bcftools view -s "$MAT","$PAT" "$VCF" | bgzip > "${OUT}.vcf.gz"
  bcftools query -f "%CHROM\t%POS\t[%TGT\t]\n" "${OUT}.vcf.gz" \
      | sed -e "s/\//\t/g" \
      | awk '{if($3==$4 && $5==$6 && $3!=$5){print"chr"$1"\t"$2"\t"$3"\t"$5}}' \
      | grep -v "\." > "${OUT}.SNP.tsv"
  awk '{if($1!="chrX" && $1!="chrY"){print$0}}' "${OUT}.SNP.tsv" > "${OUT}.SNP.rmsex.tsv"
else
  echo 'No VCF detected, skip this step'
fi

echo "[Step 2] Mapping Hi-C data to reference genome..."
seqkit seq -n "$fasta"  | grep -v "_" | grep -v 'chrM' > "${OUT_DIR}/keep_ids.txt"
seqkit grep -f "${OUT_DIR}/keep_ids.txt" "$fasta" > "${OUT_DIR}/clean.fa"
bwa index "${OUT_DIR}/clean.fa" > /dev/null 2>&1
bwa mem -5SP -t "$THREAD" "${OUT_DIR}/clean.fa" "$R1" "$R2" 2>"${OUT}.bwa.err" | gzip > "${OUT}.sam.gz"

echo "[Step 3] Phasing and generating contact pairs..."
if [ -n "$VCF" ]; then
  $hickit.js sam2seg -v "${OUT}.SNP.rmsex.tsv" "${OUT}.sam.gz" 2>"${OUT}.phase.err" | \
  $hickit.js chronly -y - 2>>"${OUT}.phase.err" | gzip > "${OUT}.contacts.seg.gz"
else
  $hickit.js sam2seg "${OUT}.sam.gz" 2>"${OUT}.phase.err" | \
  $hickit.js chronly -y - 2>>"${OUT}.phase.err" | gzip > "${OUT}.contacts.seg.gz"
fi
$hickit -i "${OUT}.contacts.seg.gz" -o - | bgzip > "${OUT}.contacts.pairs.gz"
$hickit -i "${OUT}.contacts.pairs.gz" -u -o - 2>"${OUT}.impute.err" | bgzip > "${OUT}.impute.pairs.gz"

echo "[Step 4] 3D genome reconstruction (1M resolution)..."
if [ -n "$VCF" ]; then
  $hickit -i "${OUT}.impute.pairs.gz" -Sr1m -c1 -r10m -c2 -b4m -b1m -O "${OUT}.mb.1m_temp.3dg" 2>"${OUT}.3dg.err"
  grep -v 'chrX' "${OUT}.mb.1m_temp.3dg" | grep -v 'chrY' > "${OUT}.mb.1m.3dg" 

  echo "[Step 5] Post-processing and visualization..."
  bash "${SCRIPTS_DIR}/hickit_3dg_to_3dg_rescale_unit.sh" "${OUT}.mb.1m.3dg"
  bash "${SCRIPTS_DIR}/hickit_pairs_to_con.sh" "${OUT}.contacts.pairs.gz"
  bash "${SCRIPTS_DIR}/hickit_impute_pairs_to_con.sh" "${OUT}.impute.pairs.gz"

  python2.7 "${DIP_C_DIR}/dip-c" clean3 -c "${OUT}.impute.con.gz" "${OUT}.mb.1m.dip-c.3dg" > "${OUT}.mb.1m.dip-c.clean.3dg"
  python2.7 "${DIP_C_DIR}/dip-c" exp "${OUT}.mb.1m.dip-c.clean.3dg" > "${OUT}.mb.1m.dip-c.clean.exp.3dg"

  CHR_FILE="${OUT}.chr.txt"
  grep "^>" "${OUT_DIR}/clean.fa" | sed 's/^>//' | sort -V > "$CHR_FILE"

  python2.7 "${DIP_C_DIR}/dip-c" color -n "$CHR_FILE" "${OUT}.mb.1m.dip-c.clean.exp.3dg" | \
  python2.7 "${DIP_C_DIR}/dip-c" vis -c /dev/stdin "${OUT}.mb.1m.dip-c.clean.exp.3dg" | \
  sed 's/(mat)/♀/g; s/(pat)/♂/g' > "${OUT}.mb.1m.dip-c.clean.exp.cif"
else
  $hickit -i "${OUT}.impute.pairs.gz" -r1m -c1 -r10m -c2 -b4m -b1m -O "${OUT}.mb.1m_temp.3dg" 2>"${OUT}.3dg.err"
  
  echo "[Step 5] Post-processing and visualization..."
  bash "${SCRIPTS_DIR}/hickit_3dg_to_3dg_rescale_unit.sh" "${OUT}.mb.1m_temp.3dg"
  grep -v "^$" "${OUT}.mb.1m_temp.3dg" | grep -v "^#" | awk 'NF==5'  > "${OUT}.mb.1m.3dg"
  python2.7 "${DIP_C_DIR}/dip-c" exp "${OUT}.mb.1m.3dg" > "${OUT}.mb.1m.dip-c.clean.exp.3dg"
  
  CHR_FILE="${OUT}.chr.txt"
  grep "^>" "${OUT_DIR}/clean.fa" | sed 's/^>//' | sort -V > "$CHR_FILE"
  python2.7 "${DIP_C_DIR}/dip-c" vis "${OUT}.mb.1m.dip-c.clean.exp.3dg" | \
  sed 's/(mat)/♀/g; s/(pat)/♂/g' > "${OUT}.mb.1m.dip-c.clean.exp.cif"
fi

echo "Pipeline finished! Main result file: ${OUT}.mb.1m.dip-c.clean.exp.cif"

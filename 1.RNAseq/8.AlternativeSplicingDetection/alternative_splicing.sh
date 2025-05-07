#! /usr/bin/bash



# Description : Detect the differentail splicing events by rMATs based on the bam files. 
# Author      : Zheting Zhang
# Time        : 2024-11-17


# Function to display help

display_help() {
  echo "Usage: $0 [options]"
  echo ""
  echo "Options:"
  echo "  --b1 <b1 file path>        Path to the b1 file containing BAM file paths seprated by comma"
  echo "  --b2 <b2 file path>        Path to the b2 file containing BAM file paths seprated by comma"
  echo "  --gtf <GTF file path>      Path to the GTF file"
  echo "  --type <paired/single>     Type of read used in the analysis"
  echo "  --length <read length>     The length of each read"
  echo "  --thread <thread number>   Number of threads used"
  echo "  --od <output directory>    Output directory for the results"
  echo "  --config <config_file>     The configuration file"
  echo "  --help                     Display this help message"
  exit 0

}

# Parse command line arguments

while [[ $# -gt 0 ]]; do

  key="$1"
  case $key in

    --b1)
      b1_file="$2"
      shift # past argument

      shift # past value

      ;;
    --b2)
      b2_file="$2"
      shift # past argument

      shift # past value

      ;;

    --gtf)
      gtf_file="$2"
      shift # past argument

      shift # past value

      ;;

    --type)
      Type="$2"
      shift # past argument

      shift # past value

      ;;
    --length)
      Length="$2"
      shift # past argument

      shift # past value

      ;;
    --thread)
      n_thread="$2"
      shift # past argument

      shift # past value

      ;;
    --od)
      output_dir="$2"
      shift # past argument

      shift # past value

      ;;
    --config)
      config_json="$2"
      shift # past argument

      shift # past value

      ;;
    --help)
      display_help

      ;;
    *)
      echo "Error: Unknown option $1"
      display_help

      ;;
  esac

done

# Check if required arguments are provided

if [ -z "$b1_file" ] || [ -z "$b2_file" ] || [ -z "$gtf_file" ] || [ -z "$output_dir" ]; then

  echo "Error: All required arguments must be provided"
  display_help

  exit 1

fi

if [ ! -f "$config_json" ]; then

  echo "Error: GTF file $config_json does not exist"
  exit 1

fi

if [ ! -f "$gtf_file" ]; then

  echo "Error: GTF file $gtf_file does not exist"
  exit 1

fi

# Run rmats.py command
if [ ! -d "${output_dir}" ]; then
    mkdir -p ${output_dir}
fi

if [ ! -d "${output_dir}/tmp" ]; then
    mkdir -p ${output_dir}/tmp
fi

rmats_conda=$(jq -r '.conda_env.rmats' $config_json)
rmats=$(jq -r '.software.rmats' $config_json)

$rmats_conda/bin/python $rmats/rmats.py \
  --b1 "$b1_file" \
  --b2 "$b2_file" \
  --gtf "$gtf_file" \
  -t "$Type" \
  --readLength "$Length" \
  --nthread "$n_thread" \
  --od "$output_dir" \
  --tmp "$output_dir"/tmp


echo "rmats.py command has finished, output is saved in $output_dir, log files are rmats.log and rmats.err"

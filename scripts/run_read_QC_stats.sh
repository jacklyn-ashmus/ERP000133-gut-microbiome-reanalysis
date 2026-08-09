#!/bin/bash
set -euo pipefail


usage() {
    cat <<EOF
Usage:
  $0 \\
    --project_name NAME \\
    --raw_fastq_dir PATH \\
    --trimmed_fastq_dir PATH \\
    --bad_fastq_dir PATH \\
    --output_dir PATH 

Required arguments:

  --project_name NAME
      Arbitrary project name used as a prefix for output files.

  --raw_fastq_dir PATH
      Full path to the directory with raw FASTQ files (gzipped).

  --trimmed_fastq_dir PATH
      Full path to the directory with primer-trimmed FASTQ files (gzipped).

  --bad_fastq_dir PATH
      Full path to the directory with the reads that were removed from the primer-trimmed FASTQ files (gzipped).

  --output_dir SEQUENCE
      Full path to the output directory for this script.


Example:
  $0 \\ 
   --project_name ERP000133 \\ 
   --raw_fastq_dir /home/jashmus/gut_meta/data/fastqs/raw_fastqs/ \\ 
   --output_dir /home/jashmus/gut_meta/initial_qc/ \\ 
    --trimmed_fastq_dir /home/jashmus/gut_meta/data/fastqs/trimmed_reads/ \\
    --bad_fastq_dir /home/jashmus/gut_meta/data/fastqs/trimmed_out_reads/
EOF
}



PROJECT_NAME=""
RAW_FASTQ_DIR=""
TRIMMED_FASTQ_DIR=""
BAD_FASTQ_DIR=""
OUTPUT_DIR=""


while [[ $# -gt 0 ]]; do
    case "$1" in
        --project_name)
            PROJECT_NAME="$2"
            shift 2
            ;;
        --raw_fastq_dir)
            RAW_FASTQ_DIR="$2"
            shift 2
            ;;
        --trimmed_fastq_dir)
            TRIMMED_FASTQ_DIR="$2"
            shift 2
            ;;
        --bad_fastq_dir)
            BAD_FASTQ_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

# Check required arguments.
required_vars=(
    PROJECT_NAME
    RAW_FASTQ_DIR
    TRIMMED_FASTQ_DIR
    BAD_FASTQ_DIR
    OUTPUT_DIR
)

for var_name in "${required_vars[@]}"; do
    if [[ -z "${!var_name}" ]]; then
        echo "ERROR: Missing required argument for ${var_name}" >&2
        usage >&2
        exit 1
    fi
done


# make sure project dir string ends in forward slash 
if [[ "${RAW_FASTQ_DIR}" != */ ]]; then
    RAW_FASTQ_DIR="${RAW_FASTQ_DIR}/"
fi

# make sure qc dir string ends in forward slash 
if [[ "${TRIMMED_FASTQ_DIR}" != */ ]]; then
    TRIMMED_FASTQ_DIR="${TRIMMED_FASTQ_DIR}/"
fi

# make sure fastq dir string ends in forward slash 
if [[ "${OUTPUT_DIR}" != */ ]]; then
    OUTPUT_DIR="${OUTPUT_DIR}/"
fi


# make sure fastq dir string ends in forward slash 
if [[ "${BAD_FASTQ_DIR}" != */ ]]; then
    BAD_FASTQ_DIR="${BAD_FASTQ_DIR}/"
fi



####################### Directories ######################

read_qc_dir="${OUTPUT_DIR}read_qc/" # the main directory for outputs 
mkdir -p "${read_qc_dir}"


# fastqc directories
fastqc_base_dir="${read_qc_dir}fastqc/" 
mkdir -p "${fastqc_base_dir}"
fastqc_raw_dir="${fastqc_base_dir}raw/" 
mkdir -p "${fastqc_raw_dir}"
fastqc_trimmed_dir="${fastqc_base_dir}trimmed/" 
mkdir -p "${fastqc_trimmed_dir}"
fastqc_bad_dir="${fastqc_base_dir}trimmed_out_bad/" 
mkdir -p "${fastqc_bad_dir}"

# multiqc directories
multiqc_base_dir="${read_qc_dir}multiqc/" 
mkdir -p "${multiqc_base_dir}"
multiqc_raw_dir="${multiqc_base_dir}raw/" 
mkdir -p "${multiqc_raw_dir}"
multiqc_trimmed_dir="${multiqc_base_dir}trimmed/" 
mkdir -p "${multiqc_trimmed_dir}"
multiqc_bad_dir="${multiqc_base_dir}trimmed_out_bad/" 
mkdir -p "${multiqc_bad_dir}"

################### File names #############################
seqkit_stats_raw="${read_qc_dir}seqkit_stats_raw_fastq.tsv"
seqkit_stats_trimmed="${read_qc_dir}seqkit_stats_trimmed_fastq.tsv"
seqkit_stats_bad="${read_qc_dir}seqkit_stats_bad_fastq.tsv"






###################### 
# 1. run seqkit on FASTQ files
######################


# Raw 
echo "Running seqkit stats on raw FASTQ files"
shopt -s nullglob
fastqs_raw=( "${RAW_FASTQ_DIR}"*.fastq.gz )
if (( ${#fastqs_raw[@]} == 0 )); then
    echo "ERROR: No FASTQ files found in ${RAW_FASTQ_DIR}" >&2
    exit 1
fi
seqkit stats --all --tabular "${fastqs_raw[@]}" > "${seqkit_stats_raw}"


# Trimmed 
echo "Running seqkit stats on primer-trimmed FASTQ files"
shopt -s nullglob
fastqs_trimmed=( "${TRIMMED_FASTQ_DIR}"*.fastq.gz )
if (( ${#fastqs_trimmed[@]} == 0 )); then
    echo "ERROR: No FASTQ files found in ${TRIMMED_FASTQ_DIR}" >&2
    exit 1
fi
seqkit stats --all --tabular "${fastqs_trimmed[@]}" > "${seqkit_stats_trimmed}"


# Bad
echo "Running seqkit stats on trimmed out read (bad) FASTQ files"
shopt -s nullglob
fastqs_bad=( "${BAD_FASTQ_DIR}"*.fastq.gz )
if (( ${#fastqs_bad[@]} == 0 )); then
    echo "ERROR: No FASTQ files found in ${BAD_FASTQ_DIR}" >&2
    exit 1
fi
seqkit stats --all --tabular "${fastqs_bad[@]}" > "${seqkit_stats_bad}"





###################### 
# 2. run FastQC on FASTQ files
######################

# Raw 
echo "Running FastQC on raw FASTQ files"
fastqc "${fastqs_raw[@]}" -o "${fastqc_raw_dir}"

# Trimmed 
echo "Running FastQC on trimmed FASTQ files"
fastqc "${fastqs_trimmed[@]}" -o "${fastqc_trimmed_dir}"

# Bad
echo "Running FastQC on trimmed out read (bad) FASTQ files"
fastqc "${fastqs_bad[@]}" -o "${fastqc_bad_dir}"



###################### 
# 3. run MultiQC on FastQC files
######################

# Raw 
echo "Running MultiQC on raw FASTQC files"
multiqc "${fastqc_raw_dir}" -o "${multiqc_raw_dir}"

# Trimmed 
echo "Running MultiQC on trimmed FASTQC files"
multiqc "${fastqc_trimmed_dir}" -o "${multiqc_trimmed_dir}"

# Bad
echo "Running MultiQC on trimmed out read (bad) FASTQC files"
multiqc "${fastqc_bad_dir}" -o "${multiqc_bad_dir}"





exit 0
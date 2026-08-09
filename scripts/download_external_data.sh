#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  download_external_data.sh BASE_PROJECT_DIR

Downloads:
  1. ERP000133 raw FASTQ files to:
       BASE_PROJECT_DIR/data/fastqs/raw_fastqs/
  2. SILVA 138.2 DADA2 taxonomy training files to:
       BASE_PROJECT_DIR/taxonomy_training_classifiers/

Example:
  download_external_data.sh "$HOME/gut_meta"
EOF
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
fi

if ! command -v wget >/dev/null 2>&1; then
    echo "ERROR: wget is required but was not found in PATH." >&2
    exit 1
fi

PROJECT_DIR="${1/#\~/$HOME}"
RAW_FASTQ_DIR="${PROJECT_DIR}/data/fastqs/raw_fastqs"
TAXONOMY_DIR="${PROJECT_DIR}/taxonomy_training_classifiers"

mkdir -p "${RAW_FASTQ_DIR}" "${TAXONOMY_DIR}"

RUN_ACCESSIONS=(
    ERR011073 ERR011086 ERR011076 ERR011070 ERR011067
    ERR011060 ERR011063 ERR011077 ERR011080 ERR011083
    ERR011072 ERR011082 ERR011058 ERR011062 ERR011068
    ERR011065 ERR011075 ERR011085 ERR011081 ERR011078
    ERR011084 ERR011079 ERR011069 ERR011066 ERR011074
    ERR011064 ERR011059 ERR011061 ERR011071
)

echo "Downloading ${#RUN_ACCESSIONS[@]} raw FASTQ files..."
for acc in "${RUN_ACCESSIONS[@]}"; do
    url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/${acc}/${acc}.fastq.gz"
    echo "  ${acc}"
    wget -nc -P "${RAW_FASTQ_DIR}" "${url}"
done

echo
echo "Downloading SILVA 138.2 taxonomy training sets..."
wget -nc -P "${TAXONOMY_DIR}" \
    "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz"

wget -nc -P "${TAXONOMY_DIR}" \
    "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz"

echo
echo "Downloads complete."
echo "Raw FASTQs:"
echo "  ${RAW_FASTQ_DIR}"
echo "Taxonomy training files:"
echo "  ${TAXONOMY_DIR}"

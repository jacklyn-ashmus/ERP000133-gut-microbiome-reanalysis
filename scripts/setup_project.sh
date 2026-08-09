#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  setup_project.sh BASE_PROJECT_DIR

Example:
  setup_project.sh "$HOME/gut_meta"

Creates the required project directory hierarchy.
EOF
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
fi

PROJECT_DIR="${1/#\~/$HOME}"

mkdir -p \
    "${PROJECT_DIR}/notebooks" \
    "${PROJECT_DIR}/scripts" \
    "${PROJECT_DIR}/data/fastqs/raw_fastqs" \
    "${PROJECT_DIR}/initial_qc" \
    "${PROJECT_DIR}/taxonomy_training_classifiers"

PROJECT_DIR="$(cd "${PROJECT_DIR}" && pwd)"

echo "Created project structure under:"
echo "  ${PROJECT_DIR}"
echo
echo "Required input locations:"
echo "  ${PROJECT_DIR}/ERP000133_clean_metadata.tsv"
echo "  ${PROJECT_DIR}/ERP000133_fastq_to_paper_id.tsv"
echo "  ${PROJECT_DIR}/scripts/"
echo "  ${PROJECT_DIR}/notebooks/"

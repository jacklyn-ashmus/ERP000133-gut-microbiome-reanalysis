#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

usage() {
    cat <<'EOF'
Usage:
  run_R_analysis.sh BASE_PROJECT_DIR [PROJECT_NAME]

Arguments:
  BASE_PROJECT_DIR   Absolute or ~/... path to the project root.
  PROJECT_NAME       Project prefix. Default: ERP000133

All R scripts receive the same two positional arguments:
  1. PROJECT_NAME
  2. BASE_PROJECT_DIR

The scripts are run sequentially from n1 through n9.

Example:
  run_R_analysis.sh "$HOME/gut_meta" ERP000133
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage >&2
    exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript was not found in PATH." >&2
    exit 1
fi

PROJECT_DIR="${1/#\~/$HOME}"
PROJECT_NAME="${2:-ERP000133}"

if [[ ! -d "${PROJECT_DIR}" ]]; then
    echo "ERROR: Project directory does not exist: ${PROJECT_DIR}" >&2
    exit 1
fi

PROJECT_DIR="$(cd "${PROJECT_DIR}" && pwd)"
SCRIPT_DIR="${PROJECT_DIR}/scripts"
LOG_DIR="${PROJECT_DIR}/logs"
mkdir -p "${LOG_DIR}"

find_one_script() {
    local pattern="$1"
    local matches=( "${SCRIPT_DIR}"/${pattern} )

    if (( ${#matches[@]} == 0 )); then
        echo "ERROR: No script matching '${pattern}' was found in ${SCRIPT_DIR}" >&2
        exit 1
    fi

    if (( ${#matches[@]} > 1 )); then
        echo "ERROR: More than one script matches '${pattern}':" >&2
        printf '  %s\n' "${matches[@]}" >&2
        echo "Keep only one copy of each analysis script in scripts/." >&2
        exit 1
    fi

    printf '%s\n' "${matches[0]}"
}

SCRIPT_PATTERNS=(
    "n1_Create_qc_summary_tables*.r"
    "n2_Assign_Taxonomy*.r"
    "n3_Find_Taxa_Abundance_Diversity*.r"
    "n4_Make_Phylogeny_Trees*.r"
    "n5_analyze_alpha_diversity*.r"
    "n6_analyze_beta_diversity*.r"
    "n7_Differential_Abundance_Analysis*.r"
    "n8_Age_x_Pop_Analysis*.r"
    "n9_CoreMicrobiomeComparison*.r"
)

echo "Project name: ${PROJECT_NAME}"
echo "Project dir:  ${PROJECT_DIR}"
echo

for pattern in "${SCRIPT_PATTERNS[@]}"; do
    script="$(find_one_script "${pattern}")"
    base="$(basename "${script}")"
    log="${LOG_DIR}/${base%.r}.log"

    echo "============================================================"
    echo "Running: ${base}"
    echo "Log:     ${log}"
    echo "============================================================"

    Rscript "${script}" "${PROJECT_NAME}" "${PROJECT_DIR}" 2>&1 | tee "${log}"
    echo
done

echo "All R analysis scripts completed successfully."

#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  generate_manifest.sh BASE_PROJECT_DIR

Creates:
  BASE_PROJECT_DIR/manifest.tsv

The sample-to-run mapping is hardcoded for the ERP000133 subset used in
this analysis, while the absolute FASTQ paths are generated from the
local BASE_PROJECT_DIR.

Example:
  generate_manifest.sh "$HOME/gut_meta"
EOF
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
fi

PROJECT_DIR="${1/#\~/$HOME}"

if [[ ! -d "${PROJECT_DIR}" ]]; then
    echo "ERROR: Project directory does not exist: ${PROJECT_DIR}" >&2
    exit 1
fi

PROJECT_DIR="$(cd "${PROJECT_DIR}" && pwd)"
RAW_FASTQ_DIR="${PROJECT_DIR}/data/fastqs/raw_fastqs"
MANIFEST="${PROJECT_DIR}/manifest.tsv"

mkdir -p "${RAW_FASTQ_DIR}"

printf 'sample-id\tabsolute-filepath\tdirection\n' > "${MANIFEST}"

while IFS=$'\t' read -r sample_id run_id; do
    [[ -z "${sample_id}" ]] && continue
    printf '%s\t%s/%s.fastq.gz\tforward\n' \
        "${sample_id}" "${RAW_FASTQ_DIR}" "${run_id}" >> "${MANIFEST}"
done <<'EOF'
2BF	ERR011075
3BF	ERR011077
4BF	ERR011079
6BF	ERR011081
7BF	ERR011083
8BF	ERR011084
9BF	ERR011086
10BF	ERR011058
11BF	ERR011060
12BF	ERR011062
13BF	ERR011064
15BF	ERR011066
16BF	ERR011067
17BF	ERR011068
1EU	ERR011072
2EU	ERR011076
3EU	ERR011078
5EU	ERR011080
6EU	ERR011082
8EU	ERR011085
10EU	ERR011059
11EU	ERR011061
12EU	ERR011063
13EU	ERR011065
17EU	ERR011069
18EU	ERR011070
19EU	ERR011071
20EU	ERR011073
21EU	ERR011074
EOF

echo "Manifest written to:"
echo "  ${MANIFEST}"

missing=0
while IFS=$'\t' read -r sample_id fastq direction; do
    [[ "${sample_id}" == "sample-id" ]] && continue
    if [[ ! -f "${fastq}" ]]; then
        echo "WARNING: FASTQ not found yet: ${fastq}" >&2
        missing=$((missing + 1))
    fi
done < "${MANIFEST}"

if (( missing > 0 )); then
    echo "Manifest is valid, but ${missing} FASTQ file(s) are currently missing." >&2
    echo "Run download_external_data.sh first if needed." >&2
else
    echo "All 29 FASTQ files referenced by the manifest are present."
fi

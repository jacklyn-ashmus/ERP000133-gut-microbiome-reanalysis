#!/usr/bin/env bash
set -euo pipefail


# usage: summarize_primer_check.sh PROJECT_NAME  MANIFEST_FILE CUTADAPT_REPORT_DIR OUTPUT_TSV
# example: ./summarize_primer_check.sh 
#                     ERP000133 \
#                    /home/jashmus/gut_meta/manifest_outside_docker.tsv   \
#                    /home/jashmus/gut_meta/initial_qc/cutadapt/ \
#                    /home/jashmus/gut_meta/initial_qc/ERP000133.cutadapt_check_summary.tsv

# how it's actually called in the run_QIIME2.sh script: 
# ./"${SUMMARIZE_PRIMER_CHECK_SCRIPT}" \
#                     "${PROJECT_NAME}" \
#                     "${MANIFEST_FILE}" \
#                     "${cutadapt_out_dir}" \
#                     "${cutadapt_summary_file}" 

# Required variables are inherited from the run_QIIME2.sh script that calls this one.
PROJECT_NAME="$1"
MANIFEST_FILE="$2"
CUTADAPT_REPORT_DIR="$3"
OUTPUT_TSV="$4"


# make sure CUTADAPT_REPORT_DIR string ends in forward slash 
if [[ "${CUTADAPT_REPORT_DIR}" != */ ]]; then
    CUTADAPT_REPORT_DIR="${CUTADAPT_REPORT_DIR}/"
fi





# Extract a comma-formatted integer from a Cutadapt summary line.
# Example:
#   Total reads processed: 21,736
# becomes:
#   21736
extract_count() {
    local label="$1"
    local report="$2"

    awk -v label="${label}" '
        index($0, label) == 1 {
            line = $0
            sub(/^[^:]*:[[:space:]]*/, "", line)
            split(line, fields, /[[:space:]]+/)
            gsub(/,/, "", fields[1])
            print fields[1]
            exit
        }
    ' "${report}"
}


# Extract the percentage enclosed in parentheses.
# Example:
#   Reads with adapters: 20,826 (95.8%)
# becomes:
#   95.8
extract_percent() {
    local label="$1"
    local report="$2"

    awk -v label="${label}" '
        index($0, label) == 1 {
            if (match($0, /\([0-9.]+%\)/)) {
                value = substr($0, RSTART + 1, RLENGTH - 3)
                print value
            }
            exit
        }
    ' "${report}"
}


# Extract all relevant metrics from one Cutadapt report.
#
# Results are placed into global variables prefixed with:
#   forward_
# or:
#   reverse_
parse_cutadapt_report() {
    local report="$1"
    local prefix="$2"

    local total_reads="NA"
    local reads_with_adapters="NA"
    local reads_with_adapters_pct="NA"
    local reads_passing_filters="NA"
    local reads_passing_filters_pct="NA"
    local total_bp_processed="NA"
    local total_bp_written="NA"
    local total_bp_written_pct="NA"
    local mean_read_length="NA"
    local mean_written_length="NA"

    if [[ -f "${report}" ]]; then
        total_reads="$(extract_count "Total reads processed:" "${report}")"
        reads_with_adapters="$(extract_count "Reads with adapters:" "${report}")"
        reads_with_adapters_pct="$(extract_percent "Reads with adapters:" "${report}")"
        reads_passing_filters="$(extract_count "Reads written (passing filters):" "${report}")"
        reads_passing_filters_pct="$(extract_percent "Reads written (passing filters):" "${report}")"
        total_bp_processed="$(extract_count "Total basepairs processed:" "${report}")"
        total_bp_written="$(extract_count "Total written (filtered):" "${report}")"
        total_bp_written_pct="$(extract_percent "Total written (filtered):" "${report}")"

        # Avoid division by zero and only calculate when values are numeric.
        if [[ "${total_reads}" =~ ^[0-9]+$ ]] &&
           [[ "${total_bp_processed}" =~ ^[0-9]+$ ]] &&
           (( total_reads > 0 )); then
            mean_read_length="$(
                awk -v bp="${total_bp_processed}" -v reads="${total_reads}" \
                    'BEGIN { printf "%.2f", bp / reads }'
            )"
        fi

        if [[ "${reads_passing_filters}" =~ ^[0-9]+$ ]] &&
           [[ "${total_bp_written}" =~ ^[0-9]+$ ]] &&
           (( reads_passing_filters > 0 )); then
            mean_written_length="$(
                awk -v bp="${total_bp_written}" -v reads="${reads_passing_filters}" \
                    'BEGIN { printf "%.2f", bp / reads }'
            )"
        fi
    fi

    printf -v "${prefix}_total_reads"                 '%s' "${total_reads}"
    printf -v "${prefix}_reads_with_adapters"         '%s' "${reads_with_adapters}"
    printf -v "${prefix}_reads_with_adapters_pct"     '%s' "${reads_with_adapters_pct}"
    printf -v "${prefix}_reads_passing_filters"       '%s' "${reads_passing_filters}"
    printf -v "${prefix}_reads_passing_filters_pct"   '%s' "${reads_passing_filters_pct}"
    printf -v "${prefix}_total_bp_processed"          '%s' "${total_bp_processed}"
    printf -v "${prefix}_total_bp_written"            '%s' "${total_bp_written}"
    printf -v "${prefix}_total_bp_written_pct"        '%s' "${total_bp_written_pct}"
    printf -v "${prefix}_mean_read_length"            '%s' "${mean_read_length}"
    printf -v "${prefix}_mean_written_length"         '%s' "${mean_written_length}"
}


# Write the TSV header.
printf '%s\n' \
"sample_id	run_id	fastq	direction	status	\
forward_total_reads	forward_reads_with_adapters	forward_reads_with_adapters_pct	\
forward_reads_passing_filters	forward_reads_passing_filters_pct	\
forward_total_basepairs_processed	forward_total_basepairs_written	\
forward_total_basepairs_written_pct	forward_mean_read_length	\
forward_mean_written_length	\
reverse_total_reads	reverse_reads_with_adapters	reverse_reads_with_adapters_pct	\
reverse_reads_passing_filters	reverse_reads_passing_filters_pct	\
reverse_total_basepairs_processed	reverse_total_basepairs_written	\
reverse_total_basepairs_written_pct	reverse_mean_read_length	\
reverse_mean_written_length	\
forward_report	reverse_report" \
> "${OUTPUT_TSV}"



while IFS=$'\t' read -r sample_id fastq direction; do
    [[ -z "${sample_id}" ]] && continue     # skip blank lines
    [[ "${sample_id}" == "sample-id" ]] && continue # skip header line
    [[ -z "${fastq}" ]] && continue

    # Remove common FASTQ suffixes.
    run_id="$(basename "${fastq}")"
    run_id="${run_id%.fastq.gz}"
    run_id="${run_id%.fq.gz}"
    run_id="${run_id%.fastq}"
    run_id="${run_id%.fq}"

    forward_report="${CUTADAPT_REPORT_DIR}${PROJECT_NAME}.${run_id}.${sample_id}.forward.txt"
    reverse_report="${CUTADAPT_REPORT_DIR}${PROJECT_NAME}.${run_id}.${sample_id}.reverse_comp.txt"

    status="OK"

    if [[ ! -f "${fastq}" ]]; then
        echo "WARNING: FASTQ not found for ${sample_id}: ${fastq}" >&2
        status="FASTQ_MISSING"
    fi

    if [[ ! -f "${forward_report}" ]]; then
        echo "WARNING: Forward report not found: ${forward_report}" >&2
        status="REPORT_MISSING"
    fi

    if [[ ! -f "${reverse_report}" ]]; then
        echo "WARNING: Reverse report not found: ${reverse_report}" >&2
        status="REPORT_MISSING"
    fi

    parse_cutadapt_report "${forward_report}" "forward"
    parse_cutadapt_report "${reverse_report}" "reverse"

    printf '%s\t%s\t%s\t%s\t%s\t' \
        "${sample_id}" \
        "${run_id}" \
        "${fastq}" \
        "${direction}" \
        "${status}" \
        >> "${OUTPUT_TSV}"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' \
        "${forward_total_reads}" \
        "${forward_reads_with_adapters}" \
        "${forward_reads_with_adapters_pct}" \
        "${forward_reads_passing_filters}" \
        "${forward_reads_passing_filters_pct}" \
        "${forward_total_bp_processed}" \
        "${forward_total_bp_written}" \
        "${forward_total_bp_written_pct}" \
        "${forward_mean_read_length}" \
        "${forward_mean_written_length}" \
        >> "${OUTPUT_TSV}"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' \
        "${reverse_total_reads}" \
        "${reverse_reads_with_adapters}" \
        "${reverse_reads_with_adapters_pct}" \
        "${reverse_reads_passing_filters}" \
        "${reverse_reads_passing_filters_pct}" \
        "${reverse_total_bp_processed}" \
        "${reverse_total_bp_written}" \
        "${reverse_total_bp_written_pct}" \
        "${reverse_mean_read_length}" \
        "${reverse_mean_written_length}" \
        >> "${OUTPUT_TSV}"

    printf '%s\t%s\n' \
        "${forward_report}" \
        "${reverse_report}" \
        >> "${OUTPUT_TSV}"

done < <(tail -n +2 "${MANIFEST_FILE}") # read from manifest file, skip first line



echo "     cutadapt summary written to:"
echo "          ${OUTPUT_TSV}"
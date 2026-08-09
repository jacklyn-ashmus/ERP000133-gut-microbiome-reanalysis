#!/bin/bash
set -euo pipefail


usage() {
    cat <<EOF
Usage:
  $0 \\
    --project_name NAME \\
    --manifest_file PATH \\
    --project_dir PATH \\
    --fastq_dir PATH \\
    --qc_dir PATH \\
    --forward_primer SEQUENCE \\
    --reverse_complement_primer SEQUENCE \\
    --summarize_primer_check_script PATH \\
    --forward_allowed_truncation INTEGER \\
    --reverse_allowed_truncation INTEGER \\
    --n_threads_dada2 INTEGER \\
    --n_threads_phylo INTEGER \\
    --min_read_length INTEGER \\
    --sample_metadata_tsv PATH

Required arguments:

  --project_name NAME
      Arbitrary project name used as a prefix for output files.

  --manifest_file PATH
      Full path to the manifest TSV.

      Required columns:
        sample-id
        absolute-filepath
        direction

      Each row represents one sample. The absolute-filepath column must
      contain the full path to a gzipped FASTQ file.

    --project_dir PATH
    Full path to the main project directory. The script will create a
    QIIME 2 export subdirectory at this location.
    
    --fastq_dir PATH
    Full path to the directory containing the FASTQ files. The script will
    create subdirectories here for the trimmed FASTQ files.
    
    --qc_dir PATH
    Full path to the project QC directory. The script will create the
    required Cutadapt and QIIME 2 output subdirectories here.

  --forward_primer SEQUENCE
      Forward 16S primer sequence, such as the 784F primer.

      Primer sequences should be obtained from the source publication.

      Example:
        AGGATTAGATACCCTGGTA

  --reverse_complement_primer SEQUENCE
      Reverse complement of the reverse 16S primer, such as 1061R.

      Primer sequences should be obtained from the source publication.

      Example:
        Original 1061R:
          CRRCACGAGCTGACGAC

        Reverse complement supplied to this script:
          GTCGTCAGCTCGTGYYG

  --summarize_primer_check_script PATH
      Full path to the Cutadapt all-sample summary script.

      Example:
        summarize_primer_check.sh

  --forward_allowed_truncation INTEGER
      Maximum number of bases by which the forward primer may be
      truncated while still being recognized.

      Run this workflow through the Cutadapt primer-detection and
      match-length summary steps first, then choose this value based
      on the observed forward-primer truncation distribution.

      The minimum overlap is calculated as:

        forward primer length - forward allowed truncation

  --reverse_allowed_truncation INTEGER
      Maximum number of bases by which the reverse-complement primer
      may be truncated while still being recognized.

      Run this workflow through the Cutadapt primer-detection and
      match-length summary steps first, then choose this value based
      on the observed reverse-primer truncation distribution.

      The minimum overlap is calculated as:

        reverse-complement primer length - reverse allowed truncation

  --sample_metadata_tsv PATH
      Sample metadata table containing, at the very least, columns = sample-id, 
          run_accession, sample_alias, population_group
  
  --min_read_length INTEGER
      Minimum acceptable length of a read



Optional arguments:

  --n_threads_dada2 
      The number of threads to use with QIIME DADA2. [default: 0]
      From QIIME DADA2 help manual: "If 0 is provided, all available cores will be used."

    --n_threads_phylo
      The number of threads to use with QIIME phylogeny align-to-tree-mafft-fasttree. [default: 'auto']
      From QIIME phylogeny help manual: The number of threads. (Use 'auto' to automatically
      use all available cores) This value is used when aligning the sequences and creating the tree 
      with fasttree.
    
  -h, --help
      Display this help message and exit.

Example:
  $0 \\ 
   --project_name ERP000133 \\
   --manifest_file /home/jashmus/gut_meta/manifest.tsv \\
   --project_dir /home/jashmus/gut_meta/ \\
   --fastq_dir /home/jashmus/gut_meta/data/fastqs/ \\
   --qc_dir /home/jashmus/gut_meta/initial_qc/ \\
   --forward_primer AGGATTAGATACCCTGGTA \\
   --reverse_complement_primer GTCGTCAGCTCGTGYYG \\
   --summarize_primer_check_script /home/jashmus/gut_meta/scripts/summarize_primer_check.sh \\
   --forward_allowed_truncation 3 \\
   --reverse_allowed_truncation 5 \\
   --n_threads_dada2 2 \\
   --n_threads_phylo 2 \\
   --min_read_length 150 \\
   --sample_metadata_tsv /home/jashmus/gut_meta/ERP000133_clean_metadata.tsv
EOF
}


# example: ./run_QIIME2.sh ERP000133 /home/jashmus/gut_meta/manifest_outside_docker.tsv  /home/jashmus/gut_meta/initial_qc/ AGGATTAGATACCCTGGTA GTCGTCAGCTCGTGYYG /home/jashmus/gut_meta/scripts/summarize_primer_check.sh


# Inputs definitions: 
# PROJECT_NAME - Arbitrary name for project to prefix output files with
# MANIFEST_FILE - Full path for the manifest tsv (columns = sample-id, (gzipped FASTQ) absolute-filepath, direction, rows = samples)
# QC_DIR - Main project qc results directory. This script will create the qiime2 dir.
# FORWARD_PRIMER -  784F primer. Primer sequences can be found in the source publication
# REVERSE_COMPLEMENT_PRIMER -  1061R. Primer sequences can be found in the source publication. rev= CRRCACGAGCTGACGAC -> rev_comp = GTCGTCAGCTCGTGYYG
# SUMMARIZE_PRIMER_CHECK_SCRIPT -  path for cutadapt all sample summary script, summarize_primer_check.sh
# FORWARD_ALLOWED_TRUNCATION - truncation length cutoff for trimming forward primer. Run this scripts' steps up through cutadapt rate check, then set based on truncation frequency
# REVERSE_ALLOWED_TRUNCATION- truncation length cutoff for trimming reverse primer. Run this scripts' steps up through cutadapt rate check, then set based on truncation frequency

# ERP000133 / the De Filippo Burkina Faso–Europe dataset primers: 
# FORWARD_PRIMER="AGGATTAGATACCCTGGTA"   # 784F
# REVERSE_PRIMER="CRRCACGAGCTGACGAC"    # 1061R
# REVERSE_COMPLEMENT_PRIMER="GTCGTCAGCTCGTGYYG"    # 1061R rev comp


PROJECT_NAME=""
MANIFEST_FILE=""
PROJECT_DIR=""
QC_DIR=""
FASTQ_DIR=""
FORWARD_PRIMER=""
REVERSE_COMPLEMENT_PRIMER=""
SUMMARIZE_PRIMER_CHECK_SCRIPT=""
FORWARD_ALLOWED_TRUNCATION=0
REVERSE_ALLOWED_TRUNCATION=0
NTHREADS_DADA2=0
NTHREADS_PHYLO="auto"
SAMPLE_METADATA_TABLE=""
MIN_READ_LENGTH=1


while [[ $# -gt 0 ]]; do
    case "$1" in
        --project_name)
            PROJECT_NAME="$2"
            shift 2
            ;;
        --manifest_file)
            MANIFEST_FILE="$2"
            shift 2
            ;;
        --project_dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --qc_dir)
            QC_DIR="$2"
            shift 2
            ;;
        --fastq_dir)
            FASTQ_DIR="$2"
            shift 2
            ;;
        --forward_primer)
            FORWARD_PRIMER="$2"
            shift 2
            ;;
        --reverse_complement_primer)
            REVERSE_COMPLEMENT_PRIMER="$2"
            shift 2
            ;;
        --summarize_primer_check_script)
            SUMMARIZE_PRIMER_CHECK_SCRIPT="$2"
            shift 2
            ;;
        --forward_allowed_truncation)
            FORWARD_ALLOWED_TRUNCATION="$2"
            shift 2
            ;;
        --reverse_allowed_truncation)
            REVERSE_ALLOWED_TRUNCATION="$2"
            shift 2
            ;;
        --n_threads_dada2)
            NTHREADS_DADA2="$2"
            shift 2
            ;;
        --n_threads_phylo)
            NTHREADS_PHYLO="$2"
            shift 2
            ;;
        --sample_metadata_tsv)
            SAMPLE_METADATA_TABLE="$2"
            shift 2
            ;;
        --min_read_length)
            MIN_READ_LENGTH="$2"
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
    MANIFEST_FILE
    PROJECT_DIR
    QC_DIR
    FASTQ_DIR
    FORWARD_PRIMER
    REVERSE_COMPLEMENT_PRIMER
    SUMMARIZE_PRIMER_CHECK_SCRIPT
    FORWARD_ALLOWED_TRUNCATION
    REVERSE_ALLOWED_TRUNCATION
    SAMPLE_METADATA_TABLE
    MIN_READ_LENGTH
    
)

for var_name in "${required_vars[@]}"; do
    if [[ -z "${!var_name}" ]]; then
        echo "ERROR: Missing required argument for ${var_name}" >&2
        usage >&2
        exit 1
    fi
done

# check input manifest file
if [[ ! -f "${MANIFEST_FILE}" ]]; then
    echo "ERROR: Manifest file not found: ${MANIFEST_FILE}" >&2
    exit 1
fi

if [[ ! -f "${SUMMARIZE_PRIMER_CHECK_SCRIPT}" ]]; then
    echo "ERROR: Summary script not found: ${SUMMARIZE_PRIMER_CHECK_SCRIPT}" >&2
    exit 1
fi



# calculate overlap minimums 
FORWARD_MIN_OVERLAP=$(( ${#FORWARD_PRIMER} - FORWARD_ALLOWED_TRUNCATION ))
REVERSE_MIN_OVERLAP=$(( ${#REVERSE_COMPLEMENT_PRIMER} - REVERSE_ALLOWED_TRUNCATION ))


# make sure project dir string ends in forward slash 
if [[ "${PROJECT_DIR}" != */ ]]; then
    PROJECT_DIR="${PROJECT_DIR}/"
fi

# make sure qc dir string ends in forward slash 
if [[ "${QC_DIR}" != */ ]]; then
    QC_DIR="${QC_DIR}/"
fi

# make sure fastq dir string ends in forward slash 
if [[ "${FASTQ_DIR}" != */ ]]; then
    FASTQ_DIR="${FASTQ_DIR}/"
fi



############################################## Echo inputs ##########################################################################



echo ""
echo "Inputs:"
echo "    PROJECT_NAME = project_name = "
echo "             ${PROJECT_NAME}"
echo "    MANIFEST_FILE = manifest_file = "
echo "             ${MANIFEST_FILE}"
echo "    PROJECT_DIR = project_dir = "
echo "             ${PROJECT_DIR}"
echo "    QC_DIR = qc_dir = "
echo "             ${QC_DIR}"
echo "    FASTQ_DIR = fastq_dir = "
echo "             ${FASTQ_DIR}"
echo "    FORWARD_PRIMER = forward_primer = "
echo "             ${FORWARD_PRIMER}"
echo "    REVERSE_COMPLEMENT_PRIMER = reverse_complement_primer = "
echo "             ${REVERSE_COMPLEMENT_PRIMER}"
echo "    SUMMARIZE_PRIMER_CHECK_SCRIPT = summarize_primer_check_script = "
echo "             ${SUMMARIZE_PRIMER_CHECK_SCRIPT}"
echo "    FORWARD_ALLOWED_TRUNCATION = forward_allowed_truncation = "
echo "             ${FORWARD_ALLOWED_TRUNCATION}"
echo "    REVERSE_ALLOWED_TRUNCATION = reverse_allowed_truncation = "
echo "             ${REVERSE_ALLOWED_TRUNCATION}"
echo "    FORWARD_MIN_OVERLAP = (from forward_allowed_truncation) = "
echo "             ${FORWARD_MIN_OVERLAP}"
echo "    REVERSE_MIN_OVERLAP = (from reverse_allowed_truncation) = "
echo "             ${REVERSE_MIN_OVERLAP}"
echo "    NTHREADS_DADA2 = n_threads_dada2 = "
echo "             ${NTHREADS_DADA2}"
echo "    NTHREADS_PHYLO = n_threads_phylo = "
echo "             ${NTHREADS_PHYLO}"
echo "    MIN_READ_LENGTH = min_read_length = "
echo "             ${MIN_READ_LENGTH}"
echo "    SAMPLE_METADATA_TABLE = sample_metadata_tsv = "
echo "             ${SAMPLE_METADATA_TABLE}"
echo ""


############################################## Directories ##########################################################################

# Steps 1-2
qiime_out_dir="${QC_DIR}qiime2/" # main results dir (qiime results)
demux_dir="${qiime_out_dir}demux/" # main results dir (qiime results)
# Step 3
cutadapt_base_dir="${QC_DIR}cutadapt/" # cutadapt dir
cutadapt_check_dir="${cutadapt_base_dir}rate_check_logs/" # Cutadapt check primer rate results dir
# Step 4
cutadapt_trim_out_dir="${cutadapt_base_dir}actual_trim_logs/" # Cutadapt results dir, logs of read trimming
trimmed_reads_dir="${FASTQ_DIR}trimmed_reads/" # Step 4 dir - trimmed reads
untrimmed_reads_dir="${FASTQ_DIR}trimmed_out_reads/" # Step 4 dir - reads that are removed from trimmed fastq files
# Step 5
cutadapt_doublecheck_out_dir="${cutadapt_base_dir}double_check_logs/" # Cutadapt results dir, double check after trim
# Step 7
qiime_dada2_dir="${qiime_out_dir}dada2/" # dada2 results, including denoising and ASV table files
# Step 8
qiime_phylo_dir="${qiime_out_dir}phylo/" # phylogeny file dir
# Step 9
qiime_out_export_dir="${PROJECT_DIR}qiime2_results_exported/" # Directory for R/Python readable QIIME results
dada2_asv_table_exported_dir="${qiime_out_export_dir}dada2_asv_table_export/" # Step 9 out - DADA2 ASV table export prefix
qiime_out_rep_seq_fasta_dir="${qiime_out_export_dir}rep_seq_FASTAs/" # Step 9 out - dir for DADA2 representative sequences, converted to FASTA
dada2_denoising_stats_exported_dir="${qiime_out_export_dir}dada2_denoising_stats_export/" # Step 9 out - dir for DADA2 denoising stats, R/Python readable
phylo_export_rooted_tree_dir="${qiime_out_export_dir}phylo_rooted_tree_export/" # Step 9 out - dir for Phylogeny - rooted fasttree tree, converted to Newick

############################################### Output File Paths ###################################################################


# Steps 1-2
qza_outpath="${demux_dir}${PROJECT_NAME}.qiime.1.demux.raw.qza" # Step 1 out - qza file
summary1_out="${demux_dir}${PROJECT_NAME}.qiime.1.demux.raw.summary_1.qzv" # Step 2 out - 1st summary qzv file
# Step 3.5
cutadapt_summary_file="${cutadapt_base_dir}${PROJECT_NAME}.cutadapt.check_summary.tsv"  # Step 3.5 summary table
# Step 4.5
cutadapt_trim_summary_file="${cutadapt_base_dir}${PROJECT_NAME}.cutadapt.trim_summary.tsv"  # Step 4.5 summary table
# Step 5.5
cutadapt_double_check_summary_file="${cutadapt_base_dir}${PROJECT_NAME}.cutadapt.double_check_summary.tsv"  # Step 5.5 summary table

# Step 6

postCA_qza_outpath="${demux_dir}${PROJECT_NAME}.qiime.2.cutadapt.pre_trim.qza" # Step 6 out - qza file, re-summarize after cutadapt - Untrimmed, Cutadapt analyzed reads
postCA_qza_outpath_stats="${demux_dir}${PROJECT_NAME}.qiime.2.cutadapt.pre_trim.stats.qza" # Step 6 out - stats qza file, re-summarize after cutadapt - Untrimmed, Cutadapt analyzed reads
postCA_summary2_out="${demux_dir}${PROJECT_NAME}.qiime.2.cutadapt.pre_trim.summary_2.qzv" # Step 6 out - 2nd summary qzv file, re-summarize after cutadapt  - Untrimmed, Cutadapt analyzed reads

postTrim_qza_outpath_badTrim="${demux_dir}${PROJECT_NAME}.qiime.4.cutadapt.bad.trimmed_out.qza" # Step 6 out - qza file, re-summarize after read trim - Just the bad reads
postTrim_qza_outpath_badTrim_stats="${demux_dir}${PROJECT_NAME}.qiime.4.cutadapt.bad.trimmed_out.stats.qza" # Step 6 out - stats qza file, re-summarize after cutadapt - Untrimmed, Cutadapt analyzed reads
postTrim_summary4_out_badTrim="${demux_dir}${PROJECT_NAME}.qiime.4.cutadapt.bad.trimmed_out.summary_4.qzv" # Step 6 out - 2nd summary qzv file, re-summarize after read trim - Just the bad reads

postTrim_qza_outpath="${demux_dir}${PROJECT_NAME}.qiime.3.cutadapt.trimmed.qza" # Step 6 out - qza file, re-summarize after read trim
postTrim_qza_outpath_stats="${demux_dir}${PROJECT_NAME}.qiime.3.cutadapt.trimmed.stats.qza" # Step 6 out - stats qza file, re-summarize after cutadapt - Untrimmed, Cutadapt analyzed reads
postTrim_summary3_out="${demux_dir}${PROJECT_NAME}.qiime.3.cutadapt.trimmed.summary_3.qzv" # Step 6 out - 2nd summary qzv file, re-summarize after read trim


# Step 7
dada2_asv_table="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.asv_table.qza" # Step 7 out - DADA2 ASV table
dada2_rep_seqs="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.rep_seqs.qza" # Step 7 out - DADA2 representative sequences
dada2_denoise_stats="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.denoising_stats.qza" # Step 7 out - DADA2 denoising statistics
dada2_base_transition_stats="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.base_transition_stats.qza" # Step 7 out - DADA2 base transitions statistics
# Step 7.5
dada2_denoise_stats_vis="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.denoising_stats.vis.qzv" # Step 7.5 out - DADA2 denoising statistics visualizations
# Step 7.6
dada2_asv_table_summary="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.asv_table.summary.qzv" # Step 7.6 out - DADA2 ASV table summary
dada2_asv_table_feat_freq="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.asv_table.feature_freqs" # Step 7.6 out - DADA2 ASV feature frequency table
dada2_asv_table_sample_freq="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.asv_table.sample_freqs" # Step 7.6 out - DADA2 ASV sample frequency table
# Step 7.7
dada2_rep_seqs_vis="${qiime_dada2_dir}${PROJECT_NAME}.qiime.dada2.rep_seqs.vis.qzv" # Step 7.7 out - DADA2 representative sequence visualizations
# Step 8
phylo_aligned_seqs="${qiime_phylo_dir}${PROJECT_NAME}.qiime.phylo.aligned_rep_seqs.qza" # Step 8 out - Phylogeny - aligned representative sequences
phylo_masked_aligned_seqs="${qiime_phylo_dir}${PROJECT_NAME}.qiime.phylo.aligned_rep_seqs.masked.qza" # Step 8 out - Phylogeny - aligned representative sequences, 
phylo_unrooted_tree="${qiime_phylo_dir}${PROJECT_NAME}.qiime.phylo.unrooted_tree.qza" # Step 8 out - Phylogeny - unrooted fasttree tree
phylo_rooted_tree="${qiime_phylo_dir}${PROJECT_NAME}.qiime.phylo.rooted_tree.qza" # Step 8 out - Phylogeny - rooted fasttree tree
# Step 9 
dada2_asv_table_biom="${dada2_asv_table_exported_dir}feature-table.biom" # Step 9.2 in - DADA2 ASV table, BIOM (generated in Step 9.1)
dada2_asv_table_tsv="${dada2_asv_table_exported_dir}${PROJECT_NAME}.qiime.dada2.asv_table.tsv" # Step 9.2 out - DADA2 ASV table, converted from BIOM to TSV




############################################################# Run Things ##################################################################


#######################################################
# 1. Import FASTQ files, generate qza file
#######################################################

# Make QIIME 2 dir
mkdir -p "${qiime_out_dir}"
mkdir -p "${demux_dir}"




echo ""
echo "1. Importing FASTQ files with QIIME, generating qza file"
# echo "     input file:   ${MANIFEST_FILE}    output file:  ${qza_outpath}"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "${MANIFEST_FILE}" \
  --output-path "${qza_outpath}" \
  --input-format SingleEndFastqManifestPhred33V2



# stop if that didn't work
if [[ ! -e "${qza_outpath}" ]]; then  
    echo "Error: ${qza_outpath} does not exist." >&2
    exit 1
fi
#######################################################
# 2. Generate the quality summary
#######################################################


echo "2. Summarizing Quality with QIIME Demux"
# echo "                input file:  ${qza_outpath}    output file: ${summary1_out}"

qiime demux summarize \
  --i-data "${qza_outpath}" \
  --o-visualization "${summary1_out}"





# stop if that didn't work
if [[ ! -e "${summary1_out}" ]]; then
    echo "Error: ${summary1_out} does not exist." >&2
    exit 1
fi


#######################################################
# 3. Check primer match rate with cutadapt
#######################################################


echo "3. Checking primer match rates with cutadapt."

mkdir -p "${cutadapt_base_dir}"
mkdir -p "${cutadapt_check_dir}"

while IFS=$'\t' read -r sample_id fastq direction; do
    [[ "${sample_id}" == "sample-id" ]] && continue # skip header line
    [[ -z "${fastq}" ]] && continue # skip empty fastq name

    run_id=$(basename "${fastq}" .fastq.gz) 

    if [[ ! -f "${fastq}" ]]; then
        echo "     ERROR: FASTQ not found for ${sample_id}: ${fastq}" >&2
        continue
    fi

    echo "     Checking: ${sample_id}"

    cutadapt \
        -g "forward=${FORWARD_PRIMER}" \
        --action=none \
        -o /dev/null \
        "${fastq}" \
        > "${cutadapt_check_dir}${PROJECT_NAME}.${run_id}.${sample_id}.forward.txt"

    cutadapt \
        -a "reverse_rc=${REVERSE_COMPLEMENT_PRIMER}" \
        --overlap=12 \
        --action=none \
        -o /dev/null \
        "$fastq" \
        > "${cutadapt_check_dir}${PROJECT_NAME}.${run_id}.${sample_id}.reverse_comp.txt"

done < "${MANIFEST_FILE}"


# #######################################################
# # 3.5. Get table summarizing cutadapt check
# #######################################################
# SUMMARIZE_PRIMER_CHECK_SCRIPT = gut_meta/scripts/summarize_primer_check.sh
# example: ./summarize_primer_check.sh 
#                     ERP000133 \
#                    /home/jashmus/gut_meta/manifest_outside_docker.tsv   \
#                    /home/jashmus/gut_meta/initial_qc/cutadapt_check_dir/ \
#                    /home/jashmus/gut_meta/initial_qc/ERP000133.cutadapt_check_summary.tsv




echo "3.5. Making cutadapt primer rate check summary table: ${cutadapt_summary_file}"

"${SUMMARIZE_PRIMER_CHECK_SCRIPT}" \
                    "${PROJECT_NAME}" \
                    "${MANIFEST_FILE}" \
                    "${cutadapt_check_dir}" \
                    "${cutadapt_summary_file}" 



#######################################################
# 4. Trim primers with cutadapt
#######################################################

echo "4. Trimming primers wth cutadapt."

mkdir -p "${trimmed_reads_dir}"  "${untrimmed_reads_dir}"  "${cutadapt_trim_out_dir}" 


while IFS=$'\t' read -r sample_id fastq direction; do
    [[ -z "${sample_id}" || -z "${fastq}" ]] && continue

    run_id=$(basename "${fastq}" .fastq.gz)

    if [[ ! -f "${fastq}" ]]; then
        echo "      ERROR: FASTQ not found for ${sample_id} ${run_id} : $fastq" >&2
        continue
    fi

    echo "      Trimming: ${sample_id}"

    # Pass 1: Require and remove the forward primer.
    cutadapt \
        -g "forward=${FORWARD_PRIMER}" \
        -O "${FORWARD_MIN_OVERLAP}" \
        --untrimmed-output \
            "${untrimmed_reads_dir}/${PROJECT_NAME}.${run_id}.${sample_id}.no-forward.fastq.gz" \
        -o - \
        "${fastq}" \
        2> "${cutadapt_trim_out_dir}/${PROJECT_NAME}.${run_id}.${sample_id}.forward.txt" |
    # Pass 2: Remove reverse primer when present, but retain reads without it.
    # filter out reads len < 50 after primer trim of both primers

    cutadapt \
        -a "reverse_rc=${REVERSE_COMPLEMENT_PRIMER}" \
        -O "${REVERSE_MIN_OVERLAP}" \
        -m ${MIN_READ_LENGTH} \
        -o "${trimmed_reads_dir}/${PROJECT_NAME}.${run_id}.${sample_id}.trimmed.fastq.gz" \
        - \
        > "${cutadapt_trim_out_dir}/${PROJECT_NAME}.${run_id}.${sample_id}.reverse_comp.txt"

done < <(tail -n +2 "$MANIFEST_FILE")




#######################################################
# 4.5.  Get table summarizing cutadapt trim
#######################################################
echo "4.5.  Making cutadapt trim log summary table: ${cutadapt_trim_summary_file}"

"${SUMMARIZE_PRIMER_CHECK_SCRIPT}" \
                    "${PROJECT_NAME}" \
                    "${MANIFEST_FILE}" \
                    "${cutadapt_trim_out_dir}" \
                    "${cutadapt_trim_summary_file}" 



#######################################################
# 5. Double check primer match rate after trim with cutadapt
#######################################################


echo "5. Double checking primer match rate after trim with cutadapt."

mkdir -p "${cutadapt_doublecheck_out_dir}"
shopt -s nullglob # if there aren't any fastq files, expand to an empty array instead of a literal "*".

trimmed_fastq_files=(
    "${trimmed_reads_dir%/}/${PROJECT_NAME}."*.trimmed.fastq.gz
)

if (( ${#trimmed_fastq_files[@]} == 0 )); then
    echo "ERROR: No trimmed FASTQ files found in ${trimmed_reads_dir}" >&2
    exit 1
fi

for fastq in "${trimmed_fastq_files[@]}"; do
    filename="$(basename "${fastq}")"

    identifiers="${filename#${PROJECT_NAME}.}"
    identifiers="${identifiers%.trimmed.fastq.gz}"

    run_id="${identifiers%%.*}"
    sample_id="${identifiers#*.}"

    echo "     Checking: ${sample_id}"

    cutadapt \
        -g "forward=${FORWARD_PRIMER}" \
        -O "${FORWARD_MIN_OVERLAP}" \
        --action=none \
        -o /dev/null \
        "${fastq}" \
        > "${cutadapt_doublecheck_out_dir}${PROJECT_NAME}.${run_id}.${sample_id}.forward.txt"

    cutadapt \
        -a "reverse_rc=${REVERSE_COMPLEMENT_PRIMER}" \
        -O "${REVERSE_MIN_OVERLAP}" \
        --action=none \
        -o /dev/null \
        "${fastq}" \
        > "${cutadapt_doublecheck_out_dir}${PROJECT_NAME}.${run_id}.${sample_id}.reverse_comp.txt"

done

    #   left as default: 
            # -e 0.1 


#######################################################
# 5.5. Get table summarizing cutadapt double check
#######################################################
echo "5.5. Making cutadapt trim log summary table: ${cutadapt_double_check_summary_file}"

"${SUMMARIZE_PRIMER_CHECK_SCRIPT}" \
                    "${PROJECT_NAME}" \
                    "${MANIFEST_FILE}" \
                    "${cutadapt_doublecheck_out_dir}" \
                    "${cutadapt_double_check_summary_file}" 







#######################################################
# 6.  Trim primers with QIIME cutadapt just so there's a trimmed qza file to summarize and denoise
#######################################################

echo "6.0. Trim primers with QIIME cutadapt just so there's a trimmed qza file to summarize and denoise"
echo "6.0.a. Reads with adapters, untrimmed"

qiime cutadapt trim-single \
  --i-demultiplexed-sequences "${qza_outpath}" \
  --p-adapter \
    "${FORWARD_PRIMER};min_overlap=${FORWARD_MIN_OVERLAP}...${REVERSE_COMPLEMENT_PRIMER};min_overlap=${REVERSE_MIN_OVERLAP}" \
  --p-minimum-length ${MIN_READ_LENGTH} \
  --p-cores ${NTHREADS_DADA2} \
  --o-trimmed-sequences "${postCA_qza_outpath}" \
  --o-stats "${postCA_qza_outpath_stats}" 

# left as default
#   --p-error-rate 0.1 

echo "6.0.b. Reads with adapters, trimmed"
qiime cutadapt trim-single \
  --i-demultiplexed-sequences "${qza_outpath}" \
  --p-adapter \
    "${FORWARD_PRIMER};min_overlap=${FORWARD_MIN_OVERLAP}...${REVERSE_COMPLEMENT_PRIMER};min_overlap=${REVERSE_MIN_OVERLAP}" \
  --p-minimum-length ${MIN_READ_LENGTH} \
  --p-cores ${NTHREADS_DADA2} \
  --p-discard-untrimmed \
  --o-trimmed-sequences "${postTrim_qza_outpath}" \
  --o-stats "${postTrim_qza_outpath_stats}" 


  


echo "6.5.c. Summarize - Reads without adapters"
qiime cutadapt trim-single \
  --i-demultiplexed-sequences "${qza_outpath}" \
  --p-adapter \
    "${FORWARD_PRIMER};min_overlap=${FORWARD_MIN_OVERLAP}...${REVERSE_COMPLEMENT_PRIMER};min_overlap=${REVERSE_MIN_OVERLAP}" \
  --p-minimum-length ${MIN_READ_LENGTH} \
  --p-cores ${NTHREADS_DADA2} \
  --p-discard-trimmed \
  --o-trimmed-sequences "${postTrim_qza_outpath_badTrim}" \
  --o-stats "${postTrim_qza_outpath_badTrim_stats}" 





# stop if that didn't work
if [[ ! -e "${postTrim_qza_outpath}" ]]; then
    echo "Error: ${postTrim_qza_outpath} does not exist." >&2
    exit 1
fi



#######################################################
# 6.5  Run QIIME summarize again
#######################################################
echo "6.5. Running QIIME summarize again, after read trim"


echo "6.5.a. Summarize - Reads with adapters, untrimmed"
qiime demux summarize \
  --i-data "${postCA_qza_outpath}" \
  --o-visualization "${postCA_summary2_out}"



echo "6.5.b. Summarize - Reads with adapters, trimmed"
qiime demux summarize \
  --i-data "${postTrim_qza_outpath}" \
  --o-visualization "${postTrim_summary3_out}"

echo "6.5.c. Summarize - Reads without adapters"
qiime demux summarize \
  --i-data "${postTrim_qza_outpath_badTrim}" \
  --o-visualization "${postTrim_summary4_out_badTrim}"




# stop if that didn't work
if [[ ! -e "${postTrim_summary3_out}" ]]; then
    echo "Error: ${postTrim_summary3_out} does not exist." >&2
    exit 1
fi


#######################################################
# 7.  Denoise with DADA2
#######################################################
echo "7. Denoising with DADA2"


mkdir -p "${qiime_dada2_dir}"

qiime dada2 denoise-pyro \
  --i-demultiplexed-seqs "${postTrim_qza_outpath}" \
  --p-trunc-len 0 \
  --p-n-threads ${NTHREADS_DADA2} \
  --o-table "${dada2_asv_table}" \
  --o-representative-sequences "${dada2_rep_seqs}" \
  --o-denoising-stats "${dada2_denoise_stats}" \
  --o-base-transition-stats "${dada2_base_transition_stats}" \
  --verbose




# # --p-trunc-len INTEGER   Position at which sequences should be truncated due  to decrease in quality. This truncates the 3' end of
# #                           the of the input sequences, which will be the bases that were sequenced in the last cycles. Reads that
# #                           are shorter than this value will be discarded. If 0 is provided, no truncation or length filtering will be performed  [required]
# #   --p-n-threads NTHREADS  If 0 is provided, all available cores will be used.                           [default: 1]


# # left as default: 
# #    --p-trim-left 0
# #    --p-max-ee 2.0
# #    --p-trunc-q 2
# #    --p-max-len 0
# #    --p-pooling-method independent
# #    --p-chimera-method consensus
# #   --p-min-fold-parent-over-abundance 1.0
# #   --p-allow-one-off / --p-no-allow-one-off False
# #   --p-n-reads-learn 250000
# #   --p-hashed-feature-ids / --p-no-hashed-feature-ids True
# #   --p-retain-all-samples / --p-no-retain-all-samples True





#######################################################
# 7.5.   Tabulate DADA2 Denoising Statistics with QIIME
#######################################################
echo "7.5. Tabulating DADA2 Denoising Statistics with QIIME"

# 1. Denoising statistics
qiime metadata tabulate \
  --m-input-file ${dada2_denoise_stats} \
  --o-visualization ${dada2_denoise_stats_vis}

# 2. view visualizations
# view visualizations
# qiime tools view ${dada2_denoise_stats_vis}

#######################################################
# 7.6.  Generate ASV summary tables with QIIME
#######################################################
echo "7.6. Generating ASV summary tables with QIIME"

qiime feature-table summarize \
  --i-table "${dada2_asv_table}" \
  --m-metadata-file "${SAMPLE_METADATA_TABLE}" \
  --o-feature-frequencies "${dada2_asv_table_feat_freq}" \
  --o-sample-frequencies "${dada2_asv_table_sample_freq}" \
  --o-summary "${dada2_asv_table_summary}"


# Outputs:
#   --o-feature-frequencies ARTIFACT
#     ImmutableMetadata     Per-sample and total frequencies per feature.
#                                                                     [required]
#   --o-sample-frequencies ARTIFACT
#     ImmutableMetadata     Observed feature count and total frequencies per
#                           sample.                                   [required]
#   --o-summary VISUALIZATION
#                           Visual summary of feature table           [required]
# view visualizations
# qiime "${dada2_asv_table_summary}"

#######################################################
# 7.7.  Extracting representative sequences from DADA2 results with QIIME 
#######################################################
echo "7.7. Extracting representative sequences from DADA2 results with QIIME"

qiime feature-table tabulate-seqs \
  --i-data ${dada2_rep_seqs} \
  --o-visualization ${dada2_rep_seqs_vis}


# view visualizations
# qiime tools view ${dada2_rep_seqs_vis}



#######################################################
# 8.  Build the phylogenetic tree with QIIME phyologeny
#######################################################
echo "8. Building the phylogenetic tree with QIIME phyologeny"

mkdir -p "${qiime_phylo_dir}"


qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${dada2_rep_seqs}" \
  --p-n-threads "${NTHREADS_PHYLO}" \
  --o-alignment "${phylo_aligned_seqs}" \
  --o-masked-alignment "${phylo_masked_aligned_seqs}" \
  --o-tree "${phylo_unrooted_tree}" \
  --o-rooted-tree "${phylo_rooted_tree}"

#left as default: 
#max-gap-frequency = 1.0 = no columns are removed merely for containing many gaps.
 #min-conservation = 0.4 = column is retained only when at least one nucleotide is present in at least 40% of the sequences.

######################################################################
# 9.1.  Export data to R / Python readable formats 1: Export ASV Table
######################################################################

mkdir -p "${qiime_phylo_dir}"
mkdir -p "${qiime_dada2_dir}"

echo "9. Exporting data to R / Python readable formats:"
echo "9.1.  Exporting data to R / Python readable formats 1: Export ASV table"

mkdir -p "${qiime_out_dir}"
mkdir -p "${dada2_asv_table_exported_dir}"

qiime tools export \
  --input-path "${dada2_asv_table}" \
  --output-path "${dada2_asv_table_exported_dir}"


###########################################################################
# 9.2.  Export data to R / Python readable formats 2: Convert BIOM to TSV
###########################################################################

echo "9.2.  Exporting data to R / Python readable formats 2: Convert BIOM ASV table to TSV"


biom convert \
  -i "${dada2_asv_table_biom}" \
  -o "${dada2_asv_table_tsv}" \
  --to-tsv




#################################################################################################
# 9.3.  Export data to R / Python readable formats 3: Export ASV sequence FASTA files from qza 
#################################################################################################

echo "9.3.  Exporting data to R / Python readable formats 3: Export ASV sequence FASTA files from qza "

# ASV sequences
qiime tools export \
  --input-path ${dada2_rep_seqs} \
  --output-path "${qiime_out_rep_seq_fasta_dir}"

#################################################################################################
# 9.4.  Export data to R / Python readable formats 1: Export Denoising Statistics to Newick tree 
#################################################################################################

echo "9.4.  Exporting data to R / Python readable formats 3: Export Denoising Statistics "

mkdir -p "${dada2_asv_table_exported_dir}"

# Denoising statistics
qiime tools export \
  --input-path "${dada2_denoise_stats}" \
  --output-path "${dada2_asv_table_exported_dir}"

#################################################################################################
# 9.5.  Export data to R / Python readable formats 5: Make rooted Newick tree 
#################################################################################################

echo "9.5.  Exporting data to R / Python readable formats 5: Make rooted Newick tree "

mkdir -p "${phylo_export_rooted_tree_dir}"

qiime tools export \
  --input-path "${phylo_rooted_tree}" \
  --output-path "${phylo_export_rooted_tree_dir}"



echo "Done!!! All steps complete."

wait  # Wait for any remaining background processes to complete
exit 0
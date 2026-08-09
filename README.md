# ERP000133 Gut Microbiome Reanalysis: Execution README

This repository contains the scripts and notebooks used to reproduce the ERP000133 gut microbiome reanalysis. The workflow is intentionally described as a **manual, sequential workflow**; Snakemake and Docker are not required.

The notebooks contain the same analysis code as the corresponding R scripts and are included for convenience and interactive inspection. **For a clean command-line rerun, execute the R scripts.**

---

## Final Report

The completed analysis report is available here:

- [HTML report](report/Gut_Microbiome_Analysis_Report.html)
- [PDF report](report/Gut_Microbiome_Analysis_Report.pdf)

The corresponding notebook is available in `reports/` as well.



## 1. Required project structure

Everything should be placed under one base project directory.

For the commands below, define:

```bash
export PROJECT_NAME="ERP000133"
export BASE_PROJECT_DIR="$HOME/gut_meta"
```

The required starting hierarchy is:

```text
~/[baseProjectDir Name]/
├── ERP000133_clean_metadata.tsv
├── ERP000133_fastq_to_paper_id.tsv
├── conda_env_microbiome.yml             # main R/QC environment
├── conda_env_microbiome_qiime2.yml      # QIIME 2 preprocessing environment
├── manifest.tsv                         # generated locally
├── notebooks/                           # all notebooks
├── scripts/                             # all shell and R scripts
├── taxonomy_training_classifiers/       # downloaded locally
│   ├── silva_nr99_v138.2_toGenus_trainset.fa.gz
│   └── silva_nr99_v138.2_toSpecies_trainset.fa.gz
├── data/
│   └── fastqs/
│       └── raw_fastqs/                  # downloaded raw FASTQs
└── initial_qc/
```

Additional output directories are created automatically by the supplied scripts during execution, including `qiime2_results_exported/`, `taxonomy_assignment/`, `analysis/`, `figures/`, `report_tables/`, `summary_tables/`, and `quarto_plot_tables/`.

Create the starting directory structure with:

```bash
bash setup_project.sh "$BASE_PROJECT_DIR"
```

If `setup_project.sh` is already inside the repository's `scripts/` directory, use:

```bash
bash "$BASE_PROJECT_DIR/scripts/setup_project.sh" "$BASE_PROJECT_DIR"
```

---

## 2. Project-specific files that must be supplied

Aside from the scripts/notebooks themselves, the project-specific input tables that should be included are:

```text
ERP000133_clean_metadata.tsv
ERP000133_fastq_to_paper_id.tsv
```

Place both directly in the project root:

```text
$BASE_PROJECT_DIR/ERP000133_clean_metadata.tsv
$BASE_PROJECT_DIR/ERP000133_fastq_to_paper_id.tsv
```

`ERP000133_clean_metadata.tsv` is used directly by the QIIME 2 and R analysis steps.

`ERP000133_fastq_to_paper_id.tsv` is retained as the run-accession-to-paper-sample mapping/provenance table. The executable `manifest.tsv` is generated locally so that its FASTQ paths are correct for the user's own machine.

---

## 3. Put scripts and notebooks in place

Put all `.R`/`.r` and `.sh` files in:

```text
$BASE_PROJECT_DIR/scripts/
```

Put all notebooks in:

```text
$BASE_PROJECT_DIR/notebooks/
```

For clarity, the shell scripts can be named:

```text
run_QIIME2.sh
run_read_QC_stats.sh
summarize_primer_check.sh
setup_project.sh
download_external_data.sh
generate_manifest.sh
run_R_analysis.sh
```

The R analysis scripts should retain their `n1` through `n9` prefixes so the execution order remains obvious.

Make shell scripts executable:

```bash
chmod +x "$BASE_PROJECT_DIR"/scripts/*.sh
```

---

## 4. Create the Conda environments

Two Conda environments are supplied with the project and should be kept in the **project root**:

```text
$BASE_PROJECT_DIR/conda_env_microbiome.yml
$BASE_PROJECT_DIR/conda_env_microbiome_qiime2.yml
```

They define two separate environments:

```text
conda_env_microbiome.yml          -> environment name: microbiome
conda_env_microbiome_qiime2.yml   -> environment name: microbiome_qiime2
```

The `microbiome_qiime2` environment is used for the QIIME 2 preprocessing workflow, including QIIME 2, Cutadapt, DADA2 through QIIME 2, phylogeny construction, and related command-line tools.

The `microbiome` environment is used for read-level QC and the downstream R analysis, including SeqKit/FastQC/MultiQC and the R/Bioconductor packages used in scripts `n1` through `n9`.


Alternatively, the Conda environment can be installed using these commands: 
```

###################### Main Environment ################################################

# create the environment
conda create -n microbiome -y
conda activate microbiome
# install mamba
conda install -c conda-forge mamba	-y


# install packages 
mamba install \
  --override-channels --strict-channel-priority \
  -c conda-forge -c bioconda \
  quarto jupyterlab jupyter-cache r-irkernel \
  fastqc multiqc seqkit \
  pandas numpy scipy scikit-learn matplotlib seaborn \
  r-base r-essentials r-biocmanager r-devtools r-tictoc \
  r-tidyverse r-dplyr r-readr r-stringr r-data.table \
  r-matrix r-matrixstats r-compositions r-gridextra \
  r-ggplot2 r-ggpubr r-ggthemes r-ggfortify r-ggrepel \
  r-cowplot r-patchwork r-scales r-rcolorbrewer \
  r-pheatmap r-gplots r-upsetr r-beeswarm \
  r-future r-knitr r-rmarkdown r-complexupset \
  r-vegan r-ape r-phangorn r-rstatix r-coin r-gt \
  bioconductor-biostrings bioconductor-biomformat \
  bioconductor-decipher bioconductor-dada2 \
  bioconductor-treeSummarizedexperiment \
  bioconductor-mia bioconductor-miaviz \
  bioconductor-phyloseq bioconductor-microbiome \
  bioconductor-dirichletmultinomial \
  bioconductor-ancombc bioconductor-aldex2 \
  bioconductor-deseq2 r-qiime2r \
  bioconductor-scater bioconductor-scran \
  bioconductor-delayedarray \
  bioconductor-tidysinglecellexperiment \ 
  -y   
    
mamba install \
  --override-channels --strict-channel-priority \
  -c conda-forge -c bioconda \
  r-rbiom=2.2.1 -y
  
  
  # In R
  R 
  BiocManager::install("sechm")
  # when prompted: 
  #   Update all/some/none? [a/s/n]: n
  # Save workspace image? [y/n/c]: n
  # (Don't risk breaking the environment)
  
  
  
###################### QIIME 2 Environment ################################################

 # make new env for QIIME 2 
  conda deactivate
  conda env create --name microbiome_qiime2 --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-linux-64-conda.yml \ 
          -y
  conda install -c conda-forge mamba  -y 
 
	# note: qiime2 env file installs cutadapt

 
 
###################### RDP Environment  ################################################

 conda deactivate
 conda create -n rdp_classifier -y
 conda activate rdp_classifier
 conda install -c conda-forge mamba -y
 mamba  install -c conda-forge  openjdk=17 unzip -y
 
 cd /home/jashmus/gut_meta
 mkdir -p software
 cd software

 wget https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.1.zip

unzip rdp_classifier_2.1.zip
 
 
```





### Create the environments

From the project root:

```bash
cd "$BASE_PROJECT_DIR"

conda env create -f conda_env_microbiome.yml
conda env create -f conda_env_microbiome_qiime2.yml
```

If the environments have already been created, **do not recreate them**. Check with:

```bash
conda env list
```

You should see environments named:

```text
microbiome
microbiome_qiime2
```

To update an existing environment from the supplied YAML instead of recreating it:

```bash
conda env update -n microbiome \
  -f conda_env_microbiome.yml \
  --prune

conda env update -n microbiome_qiime2 \
  -f conda_env_microbiome_qiime2.yml \
  --prune
```

### Which environment is used for each step?

```text
Step                                      Conda environment
----------------------------------------  -------------------
Create project directories                none required
Download FASTQs/SILVA                     microbiome_qiime2
Generate manifest.tsv                     none required
run_QIIME2.sh                             microbiome_qiime2
run_read_QC_stats.sh                      microbiome
n1-n9 R analysis scripts                  microbiome
Notebooks                                 microbiome
```

The download helper only requires `wget`; the supplied `microbiome_qiime2` environment contains it, so the simplest reproducible approach is to activate that environment before downloading.

Activate an environment with:

```bash
conda activate microbiome_qiime2
```

or:

```bash
conda activate microbiome
```

You can return to the base environment with:

```bash
conda deactivate
```

---

## 5. Download the raw FASTQ files and SILVA taxonomy training sets

Run:

```bash
bash "$BASE_PROJECT_DIR/scripts/download_external_data.sh" "$BASE_PROJECT_DIR"
```

This downloads the 29 ERP000133 FASTQ files into:

```text
$BASE_PROJECT_DIR/data/fastqs/raw_fastqs/
```

and downloads:

```text
silva_nr99_v138.2_toGenus_trainset.fa.gz
silva_nr99_v138.2_toSpecies_trainset.fa.gz
```

into:

```text
$BASE_PROJECT_DIR/taxonomy_training_classifiers/
```

The download script uses `wget -nc`, so files that are already present are not downloaded again.

---

## 6. Generate `manifest.tsv`

QIIME 2 requires absolute FASTQ paths in the manifest. Because those paths depend on the local project directory, **do not distribute a manifest containing paths from another computer**.

Generate the manifest locally:

```bash
bash "$BASE_PROJECT_DIR/scripts/generate_manifest.sh" "$BASE_PROJECT_DIR"
```

This creates:

```text
$BASE_PROJECT_DIR/manifest.tsv
```

with the columns:

```text
sample-id    absolute-filepath    direction
```

and hardcodes the ERP000133 sample-to-run mapping while generating the correct absolute path for each local FASTQ.

For example, if:

```text
BASE_PROJECT_DIR=/home/user/gut_meta
```

then the first manifest path will be:

```text
/home/user/gut_meta/data/fastqs/raw_fastqs/ERR011075.fastq.gz
```

rather than a path from the original analysis computer.

---

## 7. Run the QIIME 2 preprocessing workflow

### Required command-line inputs

The supplied `run_QIIME2.sh` script requires:

```text
--project_name
--manifest_file
--project_dir
--fastq_dir
--qc_dir
--forward_primer
--reverse_complement_primer
--summarize_primer_check_script
--forward_allowed_truncation
--reverse_allowed_truncation
--sample_metadata_tsv
--min_read_length
```

It also accepts:

```text
--n_threads_dada2
--n_threads_phylo
```

For this analysis, use:

```text
Project name:                  ERP000133
Forward primer (784F):         AGGATTAGATACCCTGGTA
Reverse-complement 1061R:      GTCGTCAGCTCGTGYYG
Forward allowed truncation:    3
Reverse allowed truncation:    5
Minimum read length:           150 bp
```

Activate the QIIME 2 environment:

```bash
conda activate microbiome_qiime2
```

Then run:

```bash
bash "$BASE_PROJECT_DIR/scripts/run_QIIME2.sh" \
  --project_name "$PROJECT_NAME" \
  --manifest_file "$BASE_PROJECT_DIR/manifest.tsv" \
  --project_dir "$BASE_PROJECT_DIR" \
  --fastq_dir "$BASE_PROJECT_DIR/data/fastqs" \
  --qc_dir "$BASE_PROJECT_DIR/initial_qc" \
  --forward_primer AGGATTAGATACCCTGGTA \
  --reverse_complement_primer GTCGTCAGCTCGTGYYG \
  --summarize_primer_check_script "$BASE_PROJECT_DIR/scripts/summarize_primer_check.sh" \
  --forward_allowed_truncation 3 \
  --reverse_allowed_truncation 5 \
  --n_threads_dada2 2 \
  --n_threads_phylo 2 \
  --min_read_length 150 \
  --sample_metadata_tsv "$BASE_PROJECT_DIR/${PROJECT_NAME}_clean_metadata.tsv"
```

This script performs the FASTQ import, primer checks/trimming, DADA2 denoising, phylogenetic tree construction, and export of QIIME 2 outputs used by the downstream R scripts.

The DADA2 call uses `denoise-pyro` with no fixed 3' truncation (`--p-trunc-len 0`); other DADA2 filtering/chimera settings left at their QIIME 2 defaults are documented directly in `run_QIIME2.sh`.

---

## 8. Run read-level QC

This step should be run **after `run_QIIME2.sh`**, because the QIIME/Cutadapt workflow creates the primer-trimmed and trimmed-out FASTQ directories used here.

Activate the main analysis environment:

```bash
conda activate microbiome
```

Then run:

```bash
bash "$BASE_PROJECT_DIR/scripts/run_read_QC_stats.sh" \
  --project_name "$PROJECT_NAME" \
  --raw_fastq_dir "$BASE_PROJECT_DIR/data/fastqs/raw_fastqs" \
  --trimmed_fastq_dir "$BASE_PROJECT_DIR/data/fastqs/trimmed_reads" \
  --bad_fastq_dir "$BASE_PROJECT_DIR/data/fastqs/trimmed_out_reads" \
  --output_dir "$BASE_PROJECT_DIR/initial_qc"
```

This generates SeqKit statistics, FastQC results, and MultiQC summaries under:

```text
$BASE_PROJECT_DIR/initial_qc/read_qc/
```

---

## 9. Run the R analysis scripts sequentially

Every R script uses the same two positional command-line inputs:

```text
1. PROJECT_NAME
2. BASE_PROJECT_DIR
```

Equivalent direct syntax for any individual script is:

```bash
Rscript "$BASE_PROJECT_DIR/scripts/SCRIPT_NAME.r" \
  "$PROJECT_NAME" \
  "$BASE_PROJECT_DIR"
```

The required order is:

```text
n1  Create QC summary tables
n2  Assign taxonomy
n3  Taxonomic abundance/diversity analysis
n4  Make phylogeny/taxonomy tree figures
n5  Alpha-diversity analysis
n6  Beta-diversity analysis
n7  Differential-abundance analysis
n8  Age × population analysis
n9  Core-microbiome comparison
```

To execute all nine in order, first activate the main analysis environment:

```bash
conda activate microbiome
```

Then run:

```bash
bash "$BASE_PROJECT_DIR/scripts/run_R_analysis.sh" \
  "$BASE_PROJECT_DIR" \
  "$PROJECT_NAME"
```

`run_R_analysis.sh` accepts either the clean filenames or filenames with copy suffixes such as `(1)` or `(2)`, as long as there is only one matching copy of each `n1`–`n9` script in `scripts/`.

The runner stops immediately if any script fails and writes one log per R script under:

```text
$BASE_PROJECT_DIR/logs/
```

---

## 10. R packages used by the analysis

Run the R scripts in an environment containing the packages used by the supplied code. These include, among others:

```text
ANCOMBC
ape
BiocSingular
colorspace
ComplexUpset
cowplot
dada2
dplyr
ggplot2
ggpubr
ggrepel
ggthemes
gridExtra
knitr
mia
miaViz
patchwork
pheatmap
qiime2R
RColorBrewer
readr
rstatix
scales
scater
sechm
SummarizedExperiment
tibble
tictoc
tidySingleCellExperiment
tidyverse
TreeSummarizedExperiment
vegan
viridis
```

The existing environment/package specification distributed with the project should be preferred if available.

---

## 11. Minimal complete execution order

Assuming the project-specific TSV files, scripts, notebooks, and both Conda YAML files have already been copied into the correct locations:

```bash
export PROJECT_NAME="ERP000133"
export BASE_PROJECT_DIR="$HOME/gut_meta"

# 1. Create folders
bash "$BASE_PROJECT_DIR/scripts/setup_project.sh" "$BASE_PROJECT_DIR"

# 2. Create the two Conda environments (one-time setup)
cd "$BASE_PROJECT_DIR"
conda env create -f conda_env_microbiome.yml
conda env create -f conda_env_microbiome_qiime2.yml

# If they already exist, skip the two commands above.

# 3. Download FASTQs + SILVA training data
conda activate microbiome_qiime2
bash "$BASE_PROJECT_DIR/scripts/download_external_data.sh" "$BASE_PROJECT_DIR"

# 4. Generate local absolute-path QIIME manifest
bash "$BASE_PROJECT_DIR/scripts/generate_manifest.sh" "$BASE_PROJECT_DIR"

# 5. Run QIIME 2 preprocessing
bash "$BASE_PROJECT_DIR/scripts/run_QIIME2.sh" \
  --project_name "$PROJECT_NAME" \
  --manifest_file "$BASE_PROJECT_DIR/manifest.tsv" \
  --project_dir "$BASE_PROJECT_DIR" \
  --fastq_dir "$BASE_PROJECT_DIR/data/fastqs" \
  --qc_dir "$BASE_PROJECT_DIR/initial_qc" \
  --forward_primer AGGATTAGATACCCTGGTA \
  --reverse_complement_primer GTCGTCAGCTCGTGYYG \
  --summarize_primer_check_script "$BASE_PROJECT_DIR/scripts/summarize_primer_check.sh" \
  --forward_allowed_truncation 3 \
  --reverse_allowed_truncation 5 \
  --n_threads_dada2 2 \
  --n_threads_phylo 2 \
  --min_read_length 150 \
  --sample_metadata_tsv "$BASE_PROJECT_DIR/${PROJECT_NAME}_clean_metadata.tsv"

# 6. Switch to the main R/QC environment
conda deactivate
conda activate microbiome

# 7. Run SeqKit/FastQC/MultiQC
bash "$BASE_PROJECT_DIR/scripts/run_read_QC_stats.sh" \
  --project_name "$PROJECT_NAME" \
  --raw_fastq_dir "$BASE_PROJECT_DIR/data/fastqs/raw_fastqs" \
  --trimmed_fastq_dir "$BASE_PROJECT_DIR/data/fastqs/trimmed_reads" \
  --bad_fastq_dir "$BASE_PROJECT_DIR/data/fastqs/trimmed_out_reads" \
  --output_dir "$BASE_PROJECT_DIR/initial_qc"

# 8. Run R scripts n1-n9 in order
bash "$BASE_PROJECT_DIR/scripts/run_R_analysis.sh" \
  "$BASE_PROJECT_DIR" \
  "$PROJECT_NAME"
```

If the Conda environments already exist, the normal rerun starts at the download step (or later if the external data are already present).

---

## 12. Notebook use

The notebooks duplicate the analysis contained in the R scripts and are included for convenience. They are **not required** for command-line reproduction of the analysis.

For a complete reproducible rerun, the command-line R scripts should be treated as the executable analysis workflow.

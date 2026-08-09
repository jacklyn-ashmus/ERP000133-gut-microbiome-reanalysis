# ERP000133 Snakemake workflow

Place these files at the root of `gut_meta`:

- `Snakefile`
- `config.yaml`
- your existing `conda_env_microbiome.yml`
- your existing `conda_env_microbiome_qiime2.yml`

The Snakefile assumes the existing `scripts/` directory names shown in your project tree.

## Important

Apply the edits in `required_script_fixes.txt` first. Those are clean-command-line blockers in the current scripts.

## Path behavior

Nothing above `gut_meta` is hard-coded. The Snakefile resolves the directory containing itself and sets that as the Snakemake working directory. The QIIME manifest is generated at run time with absolute FASTQ paths for the environment in which Snakemake is actually executing.

For example, inside a Docker image with the project copied to `/gut_meta`, generated manifest paths will naturally look like:

    /gut_meta/data/fastqs/raw_fastqs/ERR011058.fastq.gz

No `/home/jashmus/...` path is embedded in the workflow.

## Downloads

The workflow downloads:

- ERR011058 through ERR011086 into `data/fastqs/raw_fastqs/`
- SILVA 138.2 genus training set into `taxonomy_training_classifiers/`
- SILVA 138.2 species training set into `taxonomy_training_classifiers/`

Existing files are treated as already-complete Snakemake outputs and are not downloaded again unless removed or explicitly forced.

## Run

First inspect the DAG without executing:

    snakemake --dry-run --printshellcmds

If your Docker image already contains all tools/packages in its active environment:

    snakemake --cores 8 --printshellcmds

If instead you want Snakemake to build/use the two YAML environments:

    snakemake --cores 8 --software-deployment-method conda --printshellcmds

Logs are written under `logs/`.

## Main dependency chain

    FASTQ downloads
        -> portable manifest
        -> run_QIIME2.sh
        -> run_read_QC_stats.sh
        -> n1 QC summaries
        -> n2 taxonomy assignment
        -> n3 abundance/composition
        -> n5 alpha diversity
             -> n6 beta diversity
             -> n7 differential abundance
             -> n8 age x population
             -> n9 core microbiome

`n4_Make_Phylogeny_Trees.r` branches from n2 because it only needs the taxonomy-assigned TSE.

## Notes

The download rules use retry behavior and non-empty output validation. Analysis rules write one log per stage and use actual files as dependencies, rather than relying only on rule ordering.

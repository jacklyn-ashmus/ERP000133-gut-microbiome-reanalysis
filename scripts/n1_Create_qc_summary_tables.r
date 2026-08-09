suppressPackageStartupMessages({

########### data import ###########
library(readr)

########### data wrangling ###########
library(tidyverse)
library(tibble)

########### visualization libraries ###########   

library(ggplot2)
library(scales)
library(patchwork)


    

})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
    PROJECT_NAME <- args[1]
    PROJECT_DIR  <- normalizePath(args[2])
} else {
    # defaults for notebook 
    PROJECT_NAME <- "ERP000133"
    PROJECT_DIR  <- normalizePath("..")
}

INITIAL_QC_DIR <- file.path(PROJECT_DIR, "initial_qc")
DADA2_EXPORT_DIR <- file.path(PROJECT_DIR, "qiime2_results_exported")
READ_QC_DIR <- file.path(INITIAL_QC_DIR, "read_qc")
MULTI_QC_DIR <- file.path(READ_QC_DIR, "multiqc")

# seqkit Tables
#   raw
seqkit_raw_file <- file.path(READ_QC_DIR, "seqkit_stats_raw_fastq.tsv")
#   trimmed
seqkit_trimmed_file <- file.path(READ_QC_DIR, "seqkit_stats_trimmed_fastq.tsv")


multiqc_fastqc_raw_file <- file.path(MULTI_QC_DIR, "raw", "multiqc_data", "multiqc_fastqc.txt")
multiqc_general_raw_file <- file.path(MULTI_QC_DIR, "raw", "multiqc_data", "multiqc_general_stats.txt")
multiqc_fastqc_trimmed_file <- file.path(MULTI_QC_DIR, "trimmed", "multiqc_data", "multiqc_fastqc.txt")
multiqc_general_trimmed_file <- file.path(MULTI_QC_DIR, "trimmed", "multiqc_data", "multiqc_general_stats.txt")
dada2_stats_file <- file.path(DADA2_EXPORT_DIR, "dada2_asv_table_export", "stats.tsv")

# sequencing manifest file - provides the mapping between biological 
#            sample IDs (2BF, 10EU, etc.) and the ENA run IDs (ERR011075, etc.)
manifest_file <- file.path(PROJECT_DIR, "manifest.tsv")

# Main metaata file
metadata_file <- file.path( PROJECT_DIR, paste0(PROJECT_NAME, "_clean_metadata.tsv") )


# make report output directories if they don't already exist
REPORT_TABLE_DIR <- file.path(PROJECT_DIR, "report_tables")
SUMMARY_TABLE_DIR <- file.path(PROJECT_DIR, "summary_tables") # more extensive tables not for the report go here
QUARTO_PLOT_TABLE_DIR <- file.path(PROJECT_DIR, "quarto_plot_tables")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
FIG_QC_DIR <- file.path(FIG_DIR, "qc")
dir.create(REPORT_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QUARTO_PLOT_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_QC_DIR, recursive = TRUE, showWarnings = FALSE)



manifest <- read_tsv(manifest_file, show_col_types = FALSE)
seqkit_raw <- read_tsv(seqkit_raw_file, show_col_types = FALSE)
seqkit_trimmed <- read_tsv(seqkit_trimmed_file, show_col_types = FALSE)
fastqc_raw <- read_tsv(multiqc_fastqc_raw_file, show_col_types = FALSE)
fastqc_trimmed <- read_tsv(multiqc_fastqc_trimmed_file, show_col_types = FALSE)

dada2_stats <- read_tsv(
    dada2_stats_file,
    col_types = cols(.default = col_character()), show_col_types = FALSE
    
) |>
    filter(`sample-id` != "#q2:types")  # qiime metadata has #q2:types row

# sample_id <-> ENA run_id
manifest_ids <- manifest |>
    transmute(
        sample_id = `sample-id`,
        run_id = str_extract(basename(`absolute-filepath`), "ERR[0-9]+")
    )


# raw seqkit stats
seqkit_raw_qc <- seqkit_raw |>
    mutate(run_id = str_extract(basename(file), "ERR[0-9]+")) |>
    transmute(
        run_id,
        raw_reads = num_seqs,
        raw_mean_length = avg_len,
        raw_min_length = min_len,
        raw_max_length = max_len,
        raw_Q20_percent = `Q20(%)`,
        raw_Q30_percent = `Q30(%)`,
        raw_mean_quality = AvgQual,
        raw_GC_percent = `GC(%)`,
        raw_N_bases = sum_n
    )


# trimmed seqkit stats. map run_id back to sample_id
seqkit_trimmed_qc <- seqkit_trimmed |>
    mutate(run_id = str_extract(basename(file), "ERR[0-9]+")) |>
    left_join(manifest_ids, by = "run_id") |>
    transmute(
        sample_id, run_id,
        trimmed_reads = num_seqs,
        trimmed_mean_length = avg_len,
        trimmed_min_length = min_len,
        trimmed_max_length = max_len,
        trimmed_Q20_percent = `Q20(%)`,
        trimmed_Q30_percent = `Q30(%)`,
        trimmed_mean_quality = AvgQual,
        trimmed_GC_percent = `GC(%)`,
        trimmed_N_bases = sum_n
    )


# raw fastqc columns to retain
fastqc_raw_qc <- fastqc_raw |>
    transmute(
        run_id = Sample,
        raw_duplication_percent = 100 - total_deduplicated_percentage,
        raw_fastqc_per_base_quality = per_base_sequence_quality,
        raw_fastqc_per_sequence_quality = per_sequence_quality_scores,
        raw_fastqc_length_distribution = sequence_length_distribution,
        raw_fastqc_overrepresented_sequences = overrepresented_sequences,
        raw_fastqc_adapter_content = adapter_content
    )


# trimmed fastqc columns to retain
fastqc_trimmed_qc <- fastqc_trimmed |>
    mutate(run_id = str_extract(Sample, "ERR[0-9]+")) |>
    left_join(manifest_ids, by = "run_id") |>
    transmute(
        sample_id, run_id,
        trimmed_duplication_percent = 100 - total_deduplicated_percentage,
        trimmed_fastqc_per_base_quality = per_base_sequence_quality,
        trimmed_fastqc_per_sequence_quality = per_sequence_quality_scores,
        trimmed_fastqc_length_distribution = sequence_length_distribution,
        trimmed_fastqc_overrepresented_sequences = overrepresented_sequences,
        trimmed_fastqc_adapter_content = adapter_content
    )


# dada2 stats
dada2_qc <- dada2_stats |>
    transmute(
        sample_id = `sample-id`,
        dada2_input_reads = as.numeric(input),
        dada2_filtered_reads = as.numeric(filtered),
        dada2_denoised_reads = as.numeric(denoised),
        dada2_final_reads = as.numeric(`non-chimeric`),
        dada2_filter_retention_percent = as.numeric(`percentage of input passed filter`),
        dada2_final_retention_percent = as.numeric(`percentage of input non-chimeric`)
    )


# sample metadata

metadata <- readr::read_tsv(
    file.path(
        metadata_file
    ), show_col_types = FALSE
) |>
    rename(
        sample_id = `sample-id`,
        run_id = run_accession
    )

metadata


population_colors <- c( # Non-dynamic coloring
  "BF" = "#E16A86",
  "EU" = "#00AD9A"
)

# Dynamically assign colors to values in population_group
population_levels <- unique(metadata$population_group)

population_colors <- setNames(
    grDevices::hcl.colors(
        length(population_levels),
        palette = "Dark 3"
    ),
    population_levels
)

population_colors


# join everything together
sample_qc_metadata <- manifest_ids |>
    left_join(seqkit_raw_qc, by = "run_id") |>
    left_join(seqkit_trimmed_qc, by = c("sample_id", "run_id")) |>
    left_join(fastqc_raw_qc, by = "run_id") |>
    left_join(fastqc_trimmed_qc, by = c("sample_id", "run_id")) |>
    left_join(dada2_qc, by = "sample_id") |>
    mutate(
        primer_trim_retention_percent = 100 * trimmed_reads / raw_reads,
        total_pipeline_retention_percent = 100 * dada2_final_reads / raw_reads,
        reads_lost_primer_trimming = raw_reads - trimmed_reads,
        reads_lost_dada2 = dada2_input_reads - dada2_final_reads,
        mean_length_change = trimmed_mean_length - raw_mean_length,
        GC_change = trimmed_GC_percent - raw_GC_percent
    ) |>
    select(
        sample_id, run_id,
        raw_reads, trimmed_reads,
        dada2_input_reads, dada2_filtered_reads, dada2_denoised_reads, dada2_final_reads,
        primer_trim_retention_percent, dada2_filter_retention_percent,
        dada2_final_retention_percent, total_pipeline_retention_percent,
        raw_mean_length, trimmed_mean_length, mean_length_change,
        raw_min_length, raw_max_length, trimmed_min_length, trimmed_max_length,
        raw_Q20_percent, trimmed_Q20_percent,
        raw_Q30_percent, trimmed_Q30_percent,
        raw_mean_quality, trimmed_mean_quality,
        raw_GC_percent, trimmed_GC_percent, GC_change,
        raw_duplication_percent, trimmed_duplication_percent,
        raw_N_bases, trimmed_N_bases,
        raw_fastqc_per_base_quality, trimmed_fastqc_per_base_quality,
        raw_fastqc_per_sequence_quality, trimmed_fastqc_per_sequence_quality,
        raw_fastqc_length_distribution, trimmed_fastqc_length_distribution,
        raw_fastqc_overrepresented_sequences, trimmed_fastqc_overrepresented_sequences,
        raw_fastqc_adapter_content, trimmed_fastqc_adapter_content,
        reads_lost_primer_trimming, reads_lost_dada2
    ) |>
    arrange(
        str_extract(sample_id, "[A-Za-z]+"),
        as.numeric(str_extract(sample_id, "[0-9]+"))
    )




# now do metadata
metadata_join <- metadata |>
    select(
        sample_id, run_id,

        # sample metadata
        population_group,
        age_years,
        sex,
        reported_ethnicity,
        antibiotic_history,
        delivery_mode,
        clinical_status,
        location,
        coordinates,
        diet_description,

        # accession / sequencing metadata
        sample_accession,
        experiment_accession,
        study_accession,
        instrument_platform,
        instrument_model,
        library_layout,

        # QC values reported by ENA / publication
        read_count,
        base_count,
        mean_fastq_read_length,
        fastq_bytes,
        sra_avg_spot_length,
        sra_mbases,
        sra_mbytes
    ) |>
    rename(
        ENA_read_count = read_count,
        ENA_base_count = base_count,
        ENA_mean_fastq_read_length = mean_fastq_read_length,
        ENA_fastq_bytes = fastq_bytes,
        ENA_sra_avg_spot_length = sra_avg_spot_length,
        ENA_sra_mbases = sra_mbases,
        ENA_sra_mbytes = sra_mbytes
    )


# Join everything together again
sample_qc_metadata <- manifest_ids |>
    left_join(metadata_join, by = c("sample_id", "run_id")) |>
    left_join(seqkit_raw_qc, by = "run_id") |>
    left_join(seqkit_trimmed_qc, by = c("sample_id", "run_id")) |>
    left_join(fastqc_raw_qc, by = "run_id") |>
    left_join(fastqc_trimmed_qc, by = c("sample_id", "run_id")) |>
    left_join(dada2_qc, by = "sample_id") |>
    mutate(
        primer_trim_retention_percent = 100 * trimmed_reads / raw_reads,
        total_pipeline_retention_percent = 100 * dada2_final_reads / raw_reads,
        reads_lost_primer_trimming = raw_reads - trimmed_reads,
        reads_lost_dada2 = dada2_input_reads - dada2_final_reads,
        mean_length_change = trimmed_mean_length - raw_mean_length,
        GC_change = trimmed_GC_percent - raw_GC_percent
    ) |>
    select(
        # ids
        sample_id, run_id,

        # biological / sample metadata
        population_group,
        age_years,
        sex,
        reported_ethnicity,
        antibiotic_history,
        delivery_mode,
        clinical_status,
        location,
        coordinates,
        diet_description,

        # sequencing metadata
        sample_accession,
        experiment_accession,
        study_accession,
        instrument_platform,
        instrument_model,
        library_layout,

        # my QC
        raw_reads, trimmed_reads,
        dada2_input_reads, dada2_filtered_reads,
        dada2_denoised_reads, dada2_final_reads,

        primer_trim_retention_percent,
        dada2_filter_retention_percent,
        dada2_final_retention_percent,
        total_pipeline_retention_percent,

        raw_mean_length, trimmed_mean_length, mean_length_change,
        raw_min_length, raw_max_length,
        trimmed_min_length, trimmed_max_length,

        raw_Q20_percent, trimmed_Q20_percent,
        raw_Q30_percent, trimmed_Q30_percent,
        raw_mean_quality, trimmed_mean_quality,

        raw_GC_percent, trimmed_GC_percent, GC_change,
        raw_duplication_percent, trimmed_duplication_percent,
        raw_N_bases, trimmed_N_bases,

        raw_fastqc_per_base_quality,
        trimmed_fastqc_per_base_quality,
        raw_fastqc_per_sequence_quality,
        trimmed_fastqc_per_sequence_quality,
        raw_fastqc_length_distribution,
        trimmed_fastqc_length_distribution,
        raw_fastqc_overrepresented_sequences,
        trimmed_fastqc_overrepresented_sequences,
        raw_fastqc_adapter_content,
        trimmed_fastqc_adapter_content,

        reads_lost_primer_trimming,
        reads_lost_dada2,

        # QC / sequence statistics reported by ENA
        ENA_read_count,
        ENA_base_count,
        ENA_mean_fastq_read_length,
        ENA_fastq_bytes,
        ENA_sra_avg_spot_length,
        ENA_sra_mbases,
        ENA_sra_mbytes
    ) |>
    arrange(
        str_extract(sample_id, "[A-Za-z]+"),
        as.numeric(str_extract(sample_id, "[0-9]+"))
    )

sample_qc_metadata


# sanity check for missing data
cat("samples:             ", nrow(sample_qc_metadata), "\n")
cat("unique sample ids:   ", n_distinct(sample_qc_metadata$sample_id), "\n")
cat("unique run ids:      ", n_distinct(sample_qc_metadata$run_id), "\n")
cat("missing raw counts:  ", sum(is.na(sample_qc_metadata$raw_reads)), "\n")
cat("missing trim counts: ", sum(is.na(sample_qc_metadata$trimmed_reads)), "\n")
cat("missing dada2:       ", sum(is.na(sample_qc_metadata$dada2_final_reads)), "\n")



# Smaller table - readaable within notebook
sample_qc_metadata |>
    select(
        sample_id, run_id,
        raw_reads, trimmed_reads, dada2_final_reads,
        primer_trim_retention_percent, dada2_final_retention_percent,
        total_pipeline_retention_percent,
        raw_mean_length, trimmed_mean_length,
        raw_Q30_percent, trimmed_Q30_percent,
        raw_GC_percent, trimmed_GC_percent
    ) |>
    mutate(across(where(is.numeric), ~ round(.x, 2))) 

write_tsv( 
    sample_qc_metadata,
    file = file.path(SUMMARY_TABLE_DIR, paste0(PROJECT_NAME, ".sample_qc_table_big.tsv")),
    quote = "needed"
)


sequencing_depth_table <- tibble(
    Stage = c(
        "Raw FASTQ",
        "Primer trimmed",
        "DADA2 denoised",
        "DADA2 non-chimeric"
    ),

    `Total reads` = c(
        format(sum(sample_qc_metadata$raw_reads), big.mark = ","),
        format(sum(sample_qc_metadata$trimmed_reads), big.mark = ","),
        format(sum(sample_qc_metadata$dada2_denoised_reads), big.mark = ","),
        format(sum(sample_qc_metadata$dada2_final_reads), big.mark = ",")
    ),

    `Median reads/sample (range)` = c(
        paste0(
            format(round(median(sample_qc_metadata$raw_reads)), big.mark = ","),
            " (",
            format(min(sample_qc_metadata$raw_reads), big.mark = ","), "–",
            format(max(sample_qc_metadata$raw_reads), big.mark = ","), ")"
        ),

        paste0(
            format(round(median(sample_qc_metadata$trimmed_reads)), big.mark = ","),
            " (",
            format(min(sample_qc_metadata$trimmed_reads), big.mark = ","), "–",
            format(max(sample_qc_metadata$trimmed_reads), big.mark = ","), ")"
        ),

        paste0(
            format(round(median(sample_qc_metadata$dada2_denoised_reads)), big.mark = ","),
            " (",
            format(min(sample_qc_metadata$dada2_denoised_reads), big.mark = ","), "–",
            format(max(sample_qc_metadata$dada2_denoised_reads), big.mark = ","), ")"
        ),

        paste0(
            format(round(median(sample_qc_metadata$dada2_final_reads)), big.mark = ","),
            " (",
            format(min(sample_qc_metadata$dada2_final_reads), big.mark = ","), "–",
            format(max(sample_qc_metadata$dada2_final_reads), big.mark = ","), ")"
        )
    )
)

sequencing_depth_table

sequencing_quality_table <- tibble(
    Metric = c(
        "Mean read length (bp)",
        "Read length range (bp)",
        "Mean Q30 (%)",
        "Mean GC (%)"
    ),

    `Raw FASTQ` = c(
        round(mean(sample_qc_metadata$raw_mean_length), 1),
        paste0(min(sample_qc_metadata$raw_min_length), "–", max(sample_qc_metadata$raw_max_length)),
        round(mean(sample_qc_metadata$raw_Q30_percent), 1),
        round(mean(sample_qc_metadata$raw_GC_percent), 1)
    ),

    `Primer trimmed` = c(
        round(mean(sample_qc_metadata$trimmed_mean_length), 1),
        paste0(min(sample_qc_metadata$trimmed_min_length), "–", max(sample_qc_metadata$trimmed_max_length)),
        round(mean(sample_qc_metadata$trimmed_Q30_percent), 1),
        round(mean(sample_qc_metadata$trimmed_GC_percent), 1)
    )
)

sequencing_quality_table


retention_table <- tibble(
    Metric = c(
        # "Number of samples",
        "Median primer-trim retention (%)",
        "Median DADA2 retention (%)",
        "Median total pipeline retention (%)"
    ),

    Value = c(
        # nrow(sample_qc_metadata),
        round(median(sample_qc_metadata$primer_trim_retention_percent), 0),
        round(median(sample_qc_metadata$dada2_final_retention_percent), 1),
        round(median(sample_qc_metadata$total_pipeline_retention_percent), 1)
    )
)
retention_table


sequencing_depth_table

write_tsv( 
    sequencing_depth_table,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".Summary1.sequencing_depth_table.tsv")),
    quote = "needed"
)

write_tsv( 
    sequencing_quality_table,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".Summary2.sequencing_quality_table.tsv")),
    quote = "needed"
)

write_tsv( 
    retention_table,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".Summary3.read_retention_table.tsv")),
    quote = "needed"
)

knitr::kable(
    sequencing_depth_table,
    caption = "Sequencing depth throughout preprocessing.",
    align = "ccc"
)


knitr::kable(
    sequencing_quality_table,
    caption = "Summary of sequencing quality metrics before and after primer trimming.",
    align = "lcc"
)

knitr::kable(
    retention_table,
    caption = "Read retention % throughout preprocessing.",
    align = "lc"
)



sample_qc_report_table <- sample_qc_metadata |>
    mutate(
        group = ifelse(str_detect(sample_id, "BF$"), "BF", "EU")
    ) |>
    transmute(
        Sample_ID = sample_id,
        Group = group,
        Run_ID = run_id,
        `Raw reads` = raw_reads,
        `Trimmed reads` = trimmed_reads,
        `DADA2 denoised reads` = dada2_denoised_reads,
        `DADA2 non-chimeric reads` = dada2_final_reads,
        `Overall retention (%)` = round(total_pipeline_retention_percent, 1),
        # `Raw Q30 (%)` = round(raw_Q30_percent, 1),
        # `Trimmed Q30 (%)` = round(trimmed_Q30_percent, 1)
    )

sample_qc_report_table

knitr::kable(
    sample_qc_report_table,
    caption = "Per-sample sequencing and preprocessing quality-control statistics.",
    align = c("l", "c", rep("r", 7))
)

write_tsv(
    sample_qc_report_table,
    file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".sample_qc_report.tsv"))
)





sample_metadata_tse <- sample_qc_metadata |>
    mutate(
        group = ifelse(str_detect(sample_id, "BF$"), "BF", "EU")
    ) |>
    select(
        # ids
        sample_id,
        group,
        run_id,
        
        # biological / sample metadata
        population_group,
        age_years,
        sex,
        reported_ethnicity,
        antibiotic_history,
        delivery_mode,
        clinical_status,
        location,
        coordinates,
        diet_description,

        # qc
        raw_reads,
        trimmed_reads,
        dada2_filtered_reads,
        dada2_denoised_reads,
        dada2_final_reads,

        primer_trim_retention_percent,
        dada2_final_retention_percent,
        total_pipeline_retention_percent,

        raw_Q30_percent,
        trimmed_Q30_percent,

        raw_mean_length,
        trimmed_mean_length
    )
head(sample_metadata_tse)

write_tsv( 
    sample_metadata_tse,
    file = file.path(SUMMARY_TABLE_DIR, paste0(PROJECT_NAME, ".qc_summary_table_for_TSE.tsv")),
    quote = "needed"
)



# function: does each value of this variable occur in only one population group?
meta_df <- sample_qc_metadata

check_match_to_group <- function(var) {

    tmp <- meta_df |>
        select(
            population_group,
            value = all_of(var)
        ) |>
        filter(!is.na(value))

    # skip continuous / mostly unique variables
    if (n_distinct(tmp$value) > 5) {
        return(FALSE)
    }

    # every category of this variable must belong to only one population group
    value_group_map <- tmp |>
        distinct(value, population_group)

    all(
        value_group_map |>
            count(value) |>
            pull(n) == 1
    )
}



# get table of metadata fields which are not useful for analysis because they match pop groups exactly
uninformative_metadata <- tibble(
    variable = setdiff(
        names(meta_df),
        c(
            "sample_id",
            "group",
            "run_id",
            "population_group"
        )
    )
) |>
    mutate(
        redundant_with_population_group =
            sapply(variable, check_match_to_group)
    ) |>
    filter(redundant_with_population_group)

uninformative_metadata


write_tsv( 
    uninformative_metadata,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".uninformative_metadata_variables.tsv")),
    quote = "needed"
)



sample_qc_report_table

# reshape table to long format for use in the read retention figure 


qc_long <- sample_qc_report_table |>
  select(
    Sample_ID,
    Group,
    `Raw reads`,
    `Trimmed reads`,
    `DADA2 denoised reads`,
    `DADA2 non-chimeric reads`
  ) |>
  pivot_longer(
    cols = c(
      `Raw reads`,
      `Trimmed reads`,
      `DADA2 denoised reads`,
      `DADA2 non-chimeric reads`
    ),
    names_to = "Stage",
    values_to = "Reads"
  ) |>
  mutate(
    Stage = factor(
      Stage,
      levels = c(
        "Raw reads",
        "Trimmed reads",
        "DADA2 denoised reads",
        "DADA2 non-chimeric reads"
      ),
      labels = c(
        "Raw",
        "Primer-Trimmed",
        "Denoised",
        "Non-chimeric"
      )
    ),
    Group = factor(Group, levels = c("BF", "EU"))
  )


head(qc_long)

options(repr.plot.width = 8, repr.plot.height = 8)

qc_fig1a <- ggplot(
  qc_long,
  aes(
    x = Stage,
    y = Reads,
    group = Sample_ID,
    color = Group
  )
) +
 theme_minimal(base_size = 18)  +    
  geom_line(alpha = 0.55, linewidth = 0.9) +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_color_manual(values = population_colors) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Reads Retained During Quality Control",
    x = NULL,
    y = "Number of Reads",
    color = "Group"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold")
  ) 



qc_retention <- sample_qc_report_table |>
  mutate(
    Final_retention_pct = (`DADA2 non-chimeric reads` / `Raw reads`) * 100,
    Group = factor(Group, levels = c("BF", "EU"))
  )

qc_fig1b <- ggplot(
  qc_retention,
  aes(
    x = Group,
    y = Final_retention_pct,
    color = Group,
    fill = Group
  )
) +
  geom_violin(
    alpha = 0.25,
    width = 0.8,
    trim = FALSE,
    linewidth = 0.7
  ) +
  geom_boxplot(
    width = 0.22,
    alpha = 0.5,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.07,
    size = 2.4,
    alpha = 0.9,
      color = "black"
  ) +
  scale_color_manual(values = population_colors) +
  scale_fill_manual(values = population_colors) +
  labs(
    title = "Final Read Retention,\n Per Population Group",
    x = NULL,
    y = "Non-chimeric / raw reads (%)"
  ) +
 theme_minimal(base_size = 18)  +    
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold",  hjust = 0.4),
      # axis.text.x = element_text(size = 20, face = "bold"), 
      axis.text.x = element_text( face = "bold")
  )


options(repr.plot.width = 16, repr.plot.height = 8)
pA <- qc_fig1a + qc_fig1b
pA

# save figure
ggsave(
    filename = paste0(FIG_QC_DIR, '/', PROJECT_NAME, '.qc_read_retention.png'),
    plot = pA,
    width = 12,
    height = 14,
    units = "in",
    dpi = 300,
    bg = "white"
)

ggsave(
    filename = paste0(FIG_QC_DIR, '/', PROJECT_NAME, '.qc_read_retention.pdf'),
    plot = pA,
    width = 11,
    height = 14,
    units = "in",
    bg = "white"
)
# save table for making figure in quarto report

write_tsv( 
    qc_long,
    file = file.path(QUARTO_PLOT_TABLE_DIR, paste0(PROJECT_NAME, ".Read_Retention_fig1a.qc_long.tsv")),
    quote = "needed"
)

write_tsv( 
    qc_retention,
    file = file.path(QUARTO_PLOT_TABLE_DIR, paste0(PROJECT_NAME, ".Read_Retention_fig1b.qc_retention.tsv")),
    quote = "needed"
)









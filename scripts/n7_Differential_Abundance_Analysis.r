suppressPackageStartupMessages({

########### main libraries ###########
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(mia)
library(miaViz) # mia visualizatioon
# library(ALDEx2) # Main library for this diff abund
library(ANCOMBC) # Main library for this diff abund
    
########### data manipulation tools ###########
library(tidyverse)
library(tibble)
library(scales)

########### visualization libraries ###########   
library(ggrepel)
library(patchwork)
library(scater)

library(knitr)
library(ggplot2)
library(patchwork)
library(viridis)
    
########### resource management tools ###########
# library(tictoc)

################ notes ####################
### loaded by tidyverse
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(readr)
# library(stringr)

})

warnLevel <- getOption('warn')
options(warn = -1)







args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
    PROJECT_NAME <- args[1]
    PROJECT_DIR  <- normalizePath(args[2])
} else {
    # defaults for notebook 
    PROJECT_NAME <- "ERP000133"
    PROJECT_DIR  <- normalizePath("..")
}

ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
ALPHA_DIV_DIR <-  file.path(ANALYSIS_DIR, "alpha_div")



tse <- readRDS( file = file.path(ALPHA_DIV_DIR, paste0(PROJECT_NAME, ".tse.withAlphaDiv.rds")))
# tse_genus <- readRDS( file = file.path(ALPHA_DIV_DIR, paste0(PROJECT_NAME, ".tse_genus.withAlphaDiv.rds")))
# tse_phylum <- readRDS( file = file.path(ALPHA_DIV_DIR, paste0(PROJECT_NAME, ".tse_phylum.withAlphaDiv.rds")))
metadata <- as.data.frame(colData(tse))



# make file output directory if it doesn't already exist
ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
DIFF_ABUND_DIR <-  file.path(ANALYSIS_DIR, "diff_abund")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
FIG_DIFF_ABUND_DIR <- file.path(FIG_DIR, "diff_abund")
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIFF_ABUND_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIFF_ABUND_DIR, recursive = TRUE, showWarnings = FALSE)

REPORT_TABLE_DIR <- file.path(PROJECT_DIR, "report_tables")
SUMMARY_TABLE_DIR <- file.path(PROJECT_DIR, "summary_tables") # more extensive tables not for the report go here
dir.create(REPORT_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)


QUARTO_PLOT_TABLE_DIR <- file.path(PROJECT_DIR, "quarto_plot_tables")
dir.create(QUARTO_PLOT_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

population_levels <- unique(colData(tse)$population_group)
# Dynamically assign colors to values in population_group
population_levels <- unique(colData(tse)$population_group)

population_colors <- setNames(
    grDevices::hcl.colors(
        length(population_levels),
        palette = "Dark 3"
    ),
    population_levels
)

# population_colors

# Confirm required data are present
assayNames(tse)
colnames(colData(tse))
colnames(rowData(tse))

# Raw counts must be present
assay(tse, "counts")[1:5, 1:5]

# Confirm groups
table(colData(tse)$population_group)



colData(tse)$population_group <- factor(
    colData(tse)$population_group,
    levels = c("BF", "EU")
)

table(colData(tse)$population_group)

tse_phylum <- agglomerateByRank(
    tse,
    rank = "Phylum"
)



tse_genus <- agglomerateByRank(
    tse,
    rank = "Genus"
)



# get rid of any ASVs without an associated phylum
keep_phylum <- (
    !is.na(rowData(tse_phylum)$Phylum) &
    rowData(tse_phylum)$Phylum != ""
)

tse_phylum <- tse_phylum[keep_phylum, ]



# get rid of any ASVs without an associated genus
keep_genus <- (
    !is.na(rowData(tse_genus)$Genus) &
    rowData(tse_genus)$Genus != ""
)

tse_genus <- tse_genus[keep_genus, ]





# help(ancombc)

set.seed(123)

ancom_phylum <- ancombc2(
    data = tse_phylum,
    assay_name = "counts",
    tax_level = NULL,

    fix_formula = "population_group",
    rand_formula = NULL,

    p_adj_method = "BH",
    prv_cut = 0.10, # prevalence 10% of samples
    lib_cut = 0,

    s0_perc = 0.05,
    group = "population_group", # default

    struc_zero = TRUE,
    neg_lb = FALSE,

    alpha = 0.05,

    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,

    n_cl = 1,
    verbose = TRUE
)

ancom_genus <- ancombc2(
    data = tse_genus,
    assay_name = "counts",
    tax_level = NULL,

    fix_formula = "population_group",
    rand_formula = NULL,

    p_adj_method = "BH",
    prv_cut = 0.10, # prevalence 10% of samples
    lib_cut = 0,

    s0_perc = 0.05, # default
    group = "population_group",

    struc_zero = TRUE,
    neg_lb = FALSE,

    alpha = 0.05,

    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,

    n_cl = 1,
    verbose = TRUE
)

# ancom_phylum

# ancom_genus

genus_results <- ancom_genus$res |>
    select(
        taxon,
        lfc_population_groupEU,
        se_population_groupEU,
        W_population_groupEU,
        p_population_groupEU,
        q_population_groupEU,
        diff_population_groupEU,
        passed_ss_population_groupEU,
        diff_robust_population_groupEU
    ) |>
    arrange(q_population_groupEU)

head(genus_results, 20)

sum(genus_results$diff_population_groupEU)
sum(genus_results$diff_robust_population_groupEU)



genus_results |>
    filter(diff_robust_population_groupEU) |>
    arrange(lfc_population_groupEU)

head(ancom_genus$zero_ind)

genus_robust <- genus_results |>
    filter(diff_robust_population_groupEU) |>
    mutate(
        enriched_in = if_else(
            lfc_population_groupEU > 0,
            "EU",
            "BF"
        )
    ) |>
    arrange(lfc_population_groupEU)

genus_robust

genus_structural_zeros <- ancom_genus$zero_ind |>
    filter(
        `structural_zero (population_group = BF)` |
        `structural_zero (population_group = EU)`
    ) |>
    mutate(
        enriched_in = case_when(
            `structural_zero (population_group = BF)` ~ "EU",
            `structural_zero (population_group = EU)` ~ "BF"
        )
    )


# Robust genus-level results
genus_robust <- genus_results |>
    filter(diff_robust_population_groupEU) |>
    mutate(
        enriched_in = if_else(
            lfc_population_groupEU < 0,
            "BF",
            "EU"
        )
    )

genus_robust

p_da_genus <- genus_robust |>
    mutate(
        taxon = reorder(taxon, lfc_population_groupEU)
    ) |>
    ggplot(
        aes(
            x = lfc_population_groupEU,
            y = taxon,
            color = enriched_in
        )
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "grey50"
    ) +
    geom_errorbar(
        aes(
            xmin = lfc_population_groupEU - 1.96 * se_population_groupEU,
            xmax = lfc_population_groupEU + 1.96 * se_population_groupEU
        ),
        width = 0.15
    ) +
    geom_point(size = 4) +
    scale_color_manual(values = population_colors) +
    labs(
        title = "Differentially Abundant Genera",
        subtitle = "ANCOM-BC2 robust results",
        x = "Log fold change (EU vs BF)",
        y = NULL,
        color = "Enriched in"
    ) +
    theme_minimal(base_size = 16) +
    theme(
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.position = "bottom"
    )

p_da_genus





# save figure 
ggsave(
    filename = paste0(FIG_DIFF_ABUND_DIR, '/', PROJECT_NAME, '.diffAbund.robustGenera.png'),
    plot = p_da_genus,
    width = 16,
    height = 16,
    units = "in",
    dpi = 300,
    bg = "white"
)

ggsave(
    filename = paste0(FIG_DIFF_ABUND_DIR, '/', PROJECT_NAME, '.diffAbund.robustGenera.pdf'),
    plot = p_da_genus,
    width = 16,
    height = 16,
    units = "in",
    bg = "white"
)






write_tsv( 
    genus_results,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".diffAbund.genus.tsv")),
    quote = "needed"
)


write_tsv( 
    genus_structural_zeros,
    file = file.path(SUMMARY_TABLE_DIR, paste0(PROJECT_NAME, ".diffAbund.genus.StructuralZeros.tsv")),
    quote = "needed"
)







# extract abundances table
relabund_genus <- assay(tse_genus, "relabundance")

groups <- colData(tse_genus)$population_group

# Calculate mean relative abundance within each population
genus_mean_abundance <- data.frame(
    taxon = rownames(relabund_genus),
    mean_abundance_BF = rowMeans(
        relabund_genus[, groups == "BF", drop = FALSE]
    ),
    mean_abundance_EU = rowMeans(
        relabund_genus[, groups == "EU", drop = FALSE]
    )
)

# add abundance information to structural-zero table
genus_structural_zeros_abundant <- genus_structural_zeros |>
    left_join(
        genus_mean_abundance,
        by = "taxon"
    ) |>
    mutate(
        mean_abundance = if_else(
            enriched_in == "BF",
            mean_abundance_BF,
            mean_abundance_EU
        )
    ) |>
    # Keep taxa averaging >= 1% in the group where present
    filter(mean_abundance >= 0.01) |>
    arrange(desc(mean_abundance))

genus_structural_zeros_abundant

abund_structural_zeros <- genus_structural_zeros_abundant |>
    transmute(
        Genus = taxon,
        `Associated group` = enriched_in,
        `Mean relative abundance (%)` = mean_abundance * 100
    )

abund_structural_zeros




write_tsv( 
    abund_structural_zeros,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".diffAbund.genus.abund_and_zero.tsv")),
    quote = "needed"
)




phylum_results <- ancom_phylum$res |>
    select(
        taxon,
        lfc_population_groupEU,
        se_population_groupEU,
        W_population_groupEU,
        p_population_groupEU,
        q_population_groupEU,
        diff_population_groupEU,
        passed_ss_population_groupEU,
        diff_robust_population_groupEU
    ) |>
    arrange(q_population_groupEU)

phylum_results


write_tsv( 
    phylum_results,
    file = file.path(SUMMARY_TABLE_DIR, paste0(PROJECT_NAME, ".diffAbund.phylum.tsv")),
    quote = "needed"
)



sum(phylum_results$diff_population_groupEU)
sum(phylum_results$diff_robust_population_groupEU)



phylum_results |>
    filter(diff_robust_population_groupEU) |>
    arrange(lfc_population_groupEU)

phylum_structural_zeros <- ancom_phylum$zero_ind |>
    filter(
        `structural_zero (population_group = BF)` |
        `structural_zero (population_group = EU)`
    ) |>
    mutate(
        enriched_in = case_when(
            `structural_zero (population_group = BF)` ~ "EU",
            `structural_zero (population_group = EU)` ~ "BF"
        )
    )


phylum_structural_zeros







# extract abundances table
relabund_phylum <- assay(tse_phylum, "relabundance")

groups <- colData(tse_phylum)$population_group

# Calculate mean relative abundance within each population
phylum_mean_abundance <- data.frame(
    taxon = rownames(relabund_phylum),
    mean_abundance_BF = rowMeans(
        relabund_phylum[, groups == "BF", drop = FALSE]
    ),
    mean_abundance_EU = rowMeans(
        relabund_phylum[, groups == "EU", drop = FALSE]
    )
)

# add abundance information to structural-zero table
phylum_structural_zeros_abundant <- phylum_structural_zeros |>
    left_join(
        phylum_mean_abundance,
        by = "taxon"
    ) |>
    mutate(
        mean_abundance = if_else(
            enriched_in == "BF",
            mean_abundance_BF,
            mean_abundance_EU
        )
    ) |>
    # Keep taxa averaging >= 1% in the group where present
    filter(mean_abundance >= 0.01) |>
    arrange(desc(mean_abundance))

phylum_structural_zeros_abundant







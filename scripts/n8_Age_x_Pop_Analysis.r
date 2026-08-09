suppressPackageStartupMessages({

########### main libraries ###########
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(mia)
library(miaViz) # mia visualizatioon
library(vegan) # ecol stats (beta diversity)
library(ape) # reads trees

########### data manipulation tools ###########
library(tidyverse)
library(tibble)
library(scales)

########### visualization libraries ###########   
library(ggrepel)
library(patchwork)
library(scater)

library(knitr)

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
metadata <- as.data.frame(colData(tse))



# make file output directory if it doesn't already exist
ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
AGE_X_POP_DIR <-  file.path(ANALYSIS_DIR, "age_x_pop")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
FIG_AGE_X_POP_DIR <- file.path(FIG_DIR, "age_x_pop")
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AGE_X_POP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_AGE_X_POP_DIR, recursive = TRUE, showWarnings = FALSE)

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

baseSizeA = 16
titleSizeA = 18
dotSizeA = 4
axisTitleSizeA = 15
axisTextSizeA = 15
figureTitleSizeA = 25
legendTitleSizeA = 11

unique(metadata$age_years)

metadata$age_years

hist(metadata$age_years)


metadata_age <- as.data.frame(colData(tse)) |>
    mutate(
        sample_id = rownames(.),
        age = as.numeric(age_years),
        population_group = factor(population_group)
    )

# metadata_age


# Centered Log-Ratio (CLR)  transform ASVs
tse_clr <- transformAssay(
    tse,
    assay.type = "counts",
    method = "clr",
    pseudocount = TRUE
)

# transpose the matrix 
clr_mat <- t(
    assay(tse_clr, "clr")
)


aitchison_dist <- dist(
    clr_mat,
    method = "euclidean"
)

# make sure metadata are in the same order as the distance matrix
metadata_age <- metadata_age[
    attr(aitchison_dist, "Labels"),
]

age_permanova <- adonis2(
    aitchison_dist ~ population_group * age,
    data = metadata_age,
    permutations = 9999,
    by = "margin"
)

age_permanova

# convert distance ob to matrix
aitchison_mat <- as.matrix(aitchison_dist)

# make sure metadata match matrix order
metadata_age <- metadata_age[
    rownames(aitchison_mat),
    ,
    drop = FALSE
]

# for each sample, calculate mean distance to samples
#                  from the other population


metadata_age$mean_cross_population_distance <- NA

for (i in seq_len(nrow(metadata_age))) {

    other_population <- metadata_age$population_group != metadata_age$population_group[i]

    metadata_age$mean_cross_population_distance[i] <- mean(
        aitchison_mat[i, other_population],
        na.rm = TRUE
    )
}

metadata_age$mean_cross_population_distance

# get label thing for plot 

interaction_results <- age_permanova["population_group:age", ]

interaction_label <- paste0(
    "PERMANOVA\nPopulation × Age\n",
    "R² = ", sprintf("%.3f", interaction_results$R2),
    ", F = ", sprintf("%.2f", interaction_results$F),
    ", p = ", sprintf("%.3f", interaction_results$`Pr(>F)`)
)



metadata_age

options(repr.plot.width = 12, repr.plot.height = 12)
label = "Population × Age\nR² = 0.041, F = 1.34, p = 0.073"


p_age_distance <- ggplot(
    metadata_age,
    aes(
        x = age,
        y = mean_cross_population_distance,
        color = population_group
    )
) +
    geom_point(
        size = 4,
        alpha = 0.8
    ) +
    geom_smooth(
        method = "lm",
        se = TRUE
    ) +
    scale_color_manual(values = population_colors) +
    annotate(
        "text",
        x = Inf,
        y = -Inf,
        label = interaction_label,
        hjust = 1.1,
        vjust = -1,
        size = 5.5,
        color = "black"
    ) +
    labs(
        title = "Population-Associated\nMicrobiome Divergence with Age",
        x = "Age (years)",
        y = "Mean Aitchison Distance to Other Population",
        color = "Population"
    ) +
    theme_minimal(base_size = 16) +
    theme(
        plot.title = element_text(
            face = "bold",
            hjust = 0.5,
            size = 18
        ),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
    )

p_age_distance




ggsave(
    filename = paste0(FIG_AGE_X_POP_DIR, '/', PROJECT_NAME, '.age_x_pop.meanAitchison_distance.png'),
    plot = p_age_distance,
    width = 10,
    height = 10,
    units = "in",
    dpi = 300,
    bg = "white"
)

ggsave(
    filename = paste0(FIG_AGE_X_POP_DIR, '/', PROJECT_NAME, '.age_x_pop.meanAitchison_distance.pdf'),
    plot = p_age_distance,
    width = 12,
    height = 12,
    units = "in",
    bg = "white"
)


# save table for quarto plot
write_tsv( 
    metadata_age,
    file = file.path(QUARTO_PLOT_TABLE_DIR, paste0(PROJECT_NAME, ".age_x_pop.meanAitchison_distance_with_metadata_age.tsv")),
    quote = "needed"
)


















# Test homogeneity of dispersion between populations

pop_dispersion <- betadisper(
    aitchison_dist,
    metadata_age$population_group
)

pop_dispersion_test <- permutest(
    pop_dispersion,
    permutations = 9999
)

dispersion_results <- as.data.frame(
    pop_dispersion_test$tab
)

dispersion_results$term <- rownames(dispersion_results)

dispersion_results <- dispersion_results |>
    select(
        term,
        Df,
        `Sum Sq`,
        `Mean Sq`,
        F,
        `Pr(>F)`
    )

dispersion_results


# save table for quarto plot
write_tsv( 
    dispersion_results,
    file = file.path(QUARTO_PLOT_TABLE_DIR, paste0(PROJECT_NAME, ".age_x_pop.beta_dispersion_results.tsv")),
    quote = "needed"
)












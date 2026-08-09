suppressPackageStartupMessages({

    ########### main libraries ###########
    library(SummarizedExperiment)
    library(TreeSummarizedExperiment)
    library(mia)
    library(miaViz) # mia visualizatioon
    library(tidySingleCellExperiment)
    library(ape) # reads trees
    ########### data manipulation tools ###########
    library(tidyverse)
    library(tibble)
    library(scales)
    ########### visualization libraries ###########   
    library(sechm) # Heatmap library
    library(ggpubr)
    library(pheatmap)
    library(readr)
    library(RColorBrewer)
    library(colorspace)
    library(cowplot)
    library(grid)
    library(patchwork)
    library(gridExtra)
    library(ggrepel)
    library(ggpubr)
    library(ggthemes)
    library(RColorBrewer)
    library(patchwork)
    library(cowplot)
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

TAXA_DIR <- file.path(PROJECT_DIR, "taxonomy_assignment")
tse_rds_file <- file.path(TAXA_DIR, paste0(PROJECT_NAME, ".ASV_with_taxa.TSE.rds"))


ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
ABUND_DIV_DIR <-  file.path(ANALYSIS_DIR, "abund_div")


TREE_DIR <- file.path(PROJECT_DIR, "qiime2_results_exported", "phylo_rooted_tree_export")


tse <- readRDS( file = file.path(ABUND_DIV_DIR, paste0(PROJECT_NAME, ".tse.rds")))
tse_genus <- readRDS( file = file.path(ABUND_DIV_DIR, paste0(PROJECT_NAME, ".tse_genus.rds")))
tse_phylum <- readRDS( file = file.path(ABUND_DIV_DIR, paste0(PROJECT_NAME, ".tse_phylum.rds")))
tree <- read.tree( file = file.path(TREE_DIR, "tree.nwk"))



# make file output directory if it doesn't already exist
ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
ALPHA_DIV_DIR <-  file.path(ANALYSIS_DIR, "alpha_div")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
FIG_ALPHA_DIV_DIR <- file.path(FIG_DIR, "alpha_div")
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ALPHA_DIV_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_ALPHA_DIV_DIR, recursive = TRUE, showWarnings = FALSE)

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

all(rownames(tse) %in% tree$tip.label)
all(tree$tip.label %in% rownames(tse))

rowData(tse)

# map the ASV labels onto the tree
asv_map <- setNames(
    rownames(tse),
    rowData(tse)$original_ASV_id
)
# all(tree$tip.label %in% names(asv_map))
tree$tip.label <- asv_map[tree$tip.label]

# anyNA(tree$tip.label)
# anyDuplicated(tree$tip.label)

rowTree(tse) <- tree
rowTree(tse)

# unagreggated tse
tse <- estimateDiversity(
    tse,
    index = c(
        "observed",
        "chao1",
        "shannon",
        "simpson",
        "gini_simpson",
        "inverse_simpson", 
        "faith"
    )
)

# genus
tse_genus <- estimateDiversity(
    tse_genus,
    index = c(
        "observed",
        "chao1",
        "shannon",
        "simpson",
        "gini_simpson",
        "inverse_simpson", 
        "faith"
    )
)

# phylum
tse_phylum <- estimateDiversity(
    tse_phylum,
    index = c(
        "observed",
        "chao1",
        "shannon",
        "simpson",
        "gini_simpson",
        "inverse_simpson", 
        "faith"
    )
)



# double check for singletons, doubletons
counts <- assay(tse, "counts")

colSums(counts == 1)
colSums(counts == 2)

alpha_div_df1 <- colData(tse)[, c(
    "observed",
    "chao1",
    "shannon",
        "simpson",
        "gini_simpson",
        "inverse_simpson", 
    "faith"
)]
as.data.frame(alpha_div_df1)

alpha_div_all_df <- as.data.frame(colData(tse))
head(as.data.frame(alpha_div_all_df))



make_alpha_plot <- function(
    tse1,
    metric_col1,
    group_col1,
    pop_colors1,
    title1,
    y_label1,
    ySigPush, baseSize1, titleSize1,
    group_levels1 = c("BF", "EU"),
    test_method1 = "wilcox.test"
) {
    
    alpha_div_df1 <- as.data.frame(colData(tse1))
    
    p_alpha <- ggplot(
        alpha_div_df1,
        aes(
            x = .data[[group_col1]],
            y = .data[[metric_col1]],
            fill = .data[[group_col1]]
        )
    ) +
        geom_violin(
            trim = FALSE,
            alpha = 0.7
        ) +
        geom_jitter(
            width = 0.08,
            size = 3,
            alpha = 0.8
        ) +
        stat_compare_means(
            comparisons = list(group_levels1),
            method = test_method1,
            label = "p.signif",
            fontface = "bold",
            # label.y = max(alpha_div_df1[[metric_col1]]) * 1.32,
            # tried to make this dynamic, didn't work consistently
            label.y = ySigPush,
            
            size = 7
        ) +
        scale_fill_manual(values = pop_colors1) +
        labs(
            title = title1,
            x = NULL,
            y = y_label1
        ) +
        theme_minimal(base_size = baseSize1) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = titleSize1,
                hjust = 0.5
            ),
            axis.text.y = element_text(face = "bold"),
            axis.text.x = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.position = "none"
        )
    
    return(p_alpha)
}

titleSizeSet = 18
baseSizeSet = 15


p_div_shannon <- make_alpha_plot(
                    tse1 = tse,
                    metric_col1 = "shannon",
                    group_col1 = "population_group",
                    pop_colors1 = population_colors,
                    title1 = "Shannon",
                    y_label1 = "Shannon diversity",
                    ySigPush = 6, baseSize1 = baseSizeSet, titleSize1 = titleSizeSet
                )




p_div_obs <-  make_alpha_plot(
                        tse1 = tse,
                        metric_col1 = "observed",
                        group_col1 = "population_group",
                        pop_colors1 = population_colors,
                        title1 = "Observed Richness",
                        y_label1 = "Observed ASVs",
                        ySigPush = 250, baseSize1 = baseSizeSet, titleSize1 = titleSizeSet
                    )

p_div_simpson <-  make_alpha_plot(
                    tse1 = tse,
                    metric_col1 = "simpson",
                    group_col1 = "population_group",
                    pop_colors1 = population_colors,
                    title1 = "Simpson",
                    y_label1 = "Simpson diversity",
                    ySigPush = 0.35, baseSize1 = baseSizeSet, titleSize1 = titleSizeSet
                )

p_div_gini_simpson <-  make_alpha_plot(
                    tse1 = tse,
                    metric_col1 = "gini_simpson",
                    group_col1 = "population_group",
                    pop_colors1 = population_colors,
                    title1 = "Gini-Simpson",
                    y_label1 = "Gini-Simpson diversity",
                    ySigPush = 1.2, baseSize1 = baseSizeSet, titleSize1 = titleSizeSet
                )

p_div_inverse_simpson <-  make_alpha_plot(
                    tse1 = tse,
                    metric_col1 = "inverse_simpson",
                    group_col1 = "population_group",
                    pop_colors1 = population_colors,
                    title1 = "Inverse Simpson",
                    y_label1 = "Gini-Simpson diversity",
                    ySigPush = 65, baseSize1 = baseSizeSet, titleSize1 = titleSizeSet
                )


p_div_chao1 <-  make_alpha_plot(
                    tse1 = tse,
                    metric_col1 = "chao1",
                    group_col1 = "population_group",
                    pop_colors1 = population_colors,
                    title1 = "Chao1",
                    y_label1 = "Chao1 diversity",
                    ySigPush = 255, baseSize1 = baseSizeSet, titleSize1 = titleSizeSet
                )

p_faith <- make_alpha_plot(
    tse1 = tse,
    metric_col1 = "faith",
    group_col1 = "population_group",
    pop_colors1 = population_colors,
    title1 = "Faith's Phylogenetic Diversity",
    y_label1 = "Faith's PD",
                    ySigPush = 17, baseSize1 = baseSizeSet, titleSize1 = titleSizeSet
)
# options(repr.plot.width = 16, repr.plot.height = 12)
# (p_div_shannon + p_div_inverse_simpson + p_div_chao1) / (p_div_gini_simpson + p_div_simpson + p_div_obs + p_faith)

options(repr.plot.width = 12, repr.plot.height = 12)

alpha_combined <-
    p_div_shannon +
    p_div_simpson +
    p_div_obs +
    p_faith +
    plot_layout(ncol = 2) +
    plot_annotation(
        title = "Alpha Diversity by Population\n(Non-Agglomerated Taxa)",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 26,
                hjust = 0.5
            )
        )
    )

alpha_combined

options(repr.plot.width = 12, repr.plot.height = 12)
alpha_combined_supp <-

    p_div_gini_simpson +
    p_div_inverse_simpson +
    p_div_chao1 +
    plot_layout(ncol = 2) +
    plot_annotation(
        title = "Alpha Diversity by Population\n(Non-Agglomerated Taxa)",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 24,
                hjust = 0.5
            )
        )
    )


alpha_combined_supp

alpha_div_df1 <- as.data.frame(colData(tse))



# save figure 
# main 
ggsave(
    filename = paste0(FIG_ALPHA_DIV_DIR, '/', PROJECT_NAME, '.alpha_diversity.main.png'),
    plot = alpha_combined,
    width = 12,
    height = 12,
    units = "in",
    dpi = 300,
    bg = "white"
)

ggsave(
    filename = paste0(FIG_ALPHA_DIV_DIR, '/', PROJECT_NAME, '.alpha_diversity.main.pdf'),
    plot = alpha_combined,
    width = 12,
    height = 12,
    units = "in",
    bg = "white"
)
# main 
ggsave(
    filename = paste0(FIG_ALPHA_DIV_DIR, '/', PROJECT_NAME, '.alpha_diversity.supplement.png'),
    plot = alpha_combined_supp,
    width = 12,
    height = 12,
    units = "in",
    dpi = 300,
    bg = "white"
)

ggsave(
    filename = paste0(FIG_ALPHA_DIV_DIR, '/', PROJECT_NAME, '.alpha_diversity.supplement.pdf'),
    plot = alpha_combined_supp,
    width = 12,
    height = 12,
    units = "in",
    bg = "white"
)


# save table for quarto plot
write_tsv( 
    alpha_div_df1,
    file = file.path(QUARTO_PLOT_TABLE_DIR, paste0(PROJECT_NAME, ".alpha_diversity_for_plot.tsv")),
    quote = "needed"
)




alpha_df <- as.data.frame(colData(tse))

alpha_metrics <- c(
    "observed",
    "chao1",
    "shannon",
    "simpson",
    "gini_simpson",
    "inverse_simpson", 
    "faith"
    
)

alpha_stats <- map_dfr(alpha_metrics, function(metric1) {

    # Wilcoxon test
    test1 <- wilcox.test(
        alpha_df[[metric1]] ~ alpha_df$population_group,
        exact = FALSE
    )

    # medians by population
    medians1 <- alpha_df |>
        group_by(population_group) |>
        summarise(
            median = median(.data[[metric1]], na.rm = TRUE),
            .groups = "drop"
        )

    # rank-biserial effect size
    effect1 <- rstatix::wilcox_effsize(
        data = alpha_df,
        formula = as.formula(
            paste(metric1, "~ population_group")
        )
    )

    tibble(
        metric = metric1,

        BF_median = medians1$median[
            medians1$population_group == "BF"
        ],

        EU_median = medians1$median[
            medians1$population_group == "EU"
        ],

        effect_size = effect1$effsize,
        p_value = test1$p.value
    )
}) |>
    mutate(
        p_adj_BH = p.adjust(
            p_value,
            method = "BH"
        )
    )

alpha_stats

main_metrics <- c(
    "observed",
    "shannon",
    "simpson",
    "faith"
)

alpha_stats_main <- alpha_stats |>
    filter(metric %in% main_metrics) |>
    mutate(
        p_adj_BH = p.adjust(
            p_value,
            method = "BH"
        )
    )

alpha_stats_main

alpha_stats_display <- alpha_stats_main |>
    mutate(
        metric = recode(
            metric,
            observed = "Observed Richness",
            shannon  = "Shannon Diversity",
            simpson  = "Simpson Diversity",
            faith    = "Faith's PD"
        ),
        BF_median = round(BF_median, 3),
        EU_median = round(EU_median, 3),
        effect_size = round(effect_size, 3),
        p_value = signif(p_value, 3),
        p_adj_BH = signif(p_adj_BH, 3)
    ) |>
    rename(
        "Metric" = metric,
        "BF Median" = BF_median,
        "EU Median" = EU_median,
        "Effect Size (r)" = effect_size,
        "P-value" = p_value,
        "BH-adjusted P-value" = p_adj_BH
    )



# save table for quarto plot
write_tsv( 
    alpha_stats_display,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".alpha_diversity.stats.tsv")),
    quote = "needed"
)



# save table for making figure in quarto report
saveRDS(tse, file = file.path(ALPHA_DIV_DIR, paste0(PROJECT_NAME, ".tse.withAlphaDiv.rds")))

# save table for making figure in quarto report
saveRDS(tse_genus, file = file.path(ALPHA_DIV_DIR, paste0(PROJECT_NAME, ".tse_genus.withAlphaDiv.rds")))


# save table for making figure in quarto report
saveRDS(tse_phylum, file = file.path(ALPHA_DIV_DIR, paste0(PROJECT_NAME, ".tse_phylum.withAlphaDiv.rds")))



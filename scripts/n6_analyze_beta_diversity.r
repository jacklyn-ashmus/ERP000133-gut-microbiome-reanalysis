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
tse_genus <- readRDS( file = file.path(ALPHA_DIV_DIR, paste0(PROJECT_NAME, ".tse_genus.withAlphaDiv.rds")))
tse_phylum <- readRDS( file = file.path(ALPHA_DIV_DIR, paste0(PROJECT_NAME, ".tse_phylum.withAlphaDiv.rds")))
metadata <- as.data.frame(colData(tse))



# make file output directory if it doesn't already exist
ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
BETA_DIV_DIR <-  file.path(ANALYSIS_DIR, "beta_div")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
FIG_BETA_DIV_DIR <- file.path(FIG_DIR, "beta_div")
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BETA_DIV_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_BETA_DIV_DIR, recursive = TRUE, showWarnings = FALSE)

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

beta_div_metrics = c('BrayCurtis', 'Jaccard', 'Aitchison', 'Unifrac_weighted',  'Unifrac_unweighted')


stats_list = list()
obj_list = list()
stats_list[['BrayCurtis']] = list()
stats_list[['Jaccard']] = list()
stats_list[['Aitchison']] = list()
stats_list[['Unifrac_weighted']] = list()
stats_list[['Unifrac_unweighted']] = list()
plots_list = list()





title_list = list()

title_list[['BrayCurtis']] = "Bray-Curtis PCoA"
title_list[['Jaccard']] = "Jaccard PCoA"
title_list[['Aitchison']] = "Aitchison PCA"
title_list[['Unifrac_weighted']] = "Weighted UniFrac PCoA"
title_list[['Unifrac_unweighted']] = "Unweighted UniFrac PCoA"






### PCoA
run_pcoa <- function(
    distance_object,
    metadata1
) {
    # run pcoa
    pcoa_object <- ape::pcoa(distance_object)
    # extract just the Eigenvalues, Relative_eig
    scores <- as.data.frame(
        pcoa_object$vectors[, 1:2]
    )
    # save rownames as sample, merge scores with metadata
    scores$sample <- rownames(scores)
    scores <- scores %>%
        left_join(
            metadata1 %>%
                rownames_to_column("sample"),
            by = "sample"
        )

    # extract percent variation explained
    percent_variation <- round(
        100 * pcoa_object$values$Relative_eig[1:2],
        1
    )
    return_list = list()
    return_list[['pcoa']] = pcoa_object
    return_list[['scores']] = scores
    return_list[['percent']] = percent_variation
    return(return_list)
}


### PCoA plot function

plot_pcoa <- function(
    pcoa_return1,
    metadata1,
    plot_title,
    group_var ,
    popColors
) {

    pcoa1 = pcoa_return1[['pcoa']]
    scores = pcoa_return1[['scores']]
    percent_variation = pcoa_return1[['percent']]
    
 p1 <-    ggplot(
        scores,
        aes(
            x = Axis.1,
            y = Axis.2,
            color = .data[[group_var]]
        )
    ) +
        geom_point(size = dotSizeA) +
        scale_color_manual(values = popColors) +
        labs(
            title = plot_title,
            x = paste0("PCoA1 (", percent_variation[1], "%)"),
            y = paste0("PCoA2 (", percent_variation[2], "%)"),
            color = "Population"
        ) +
        theme_minimal(base_size = baseSizeA) +
        theme(
            plot.title = element_text(
                face = "bold",
                hjust = 0.5, 
                size = titleSizeA
            ),
            axis.title = element_text(face = "bold", size = axisTitleSizeA),
            legend.title = element_text(face = "bold", size = legendTitleSizeA)
        )

    return(p1)
}


### PCA result extraction
extract_pca <- function(
    tse1,
    metadata1,
    pca_name = "PCA"
) {
    # extract PCA object from reducedDims
    pca_object <- reducedDim(
        tse1,
        pca_name
    )
    # extract PCA scores
    scores <- as.data.frame(
        pca_object
    )
    # save rownames as sample, merge scores with metadata
    scores$sample <- rownames(scores)
    scores <- scores %>%
        left_join(
            metadata1 %>%
                rownames_to_column("sample"),
            by = "sample"
        )
    # extract percent variation explained
    percent_variation <- attr(
        pca_object,
        "percentVar"
    )

    if (!is.null(percent_variation)) {
        percent_variation <- round(
            percent_variation,
            1
        )
    }
    # extract all PCA attributes
    pca_attributes <- attributes(
        pca_object
    )
    return_list = list()

    return_list[['pca']] = pca_object
    return_list[['scores']] = scores
    return_list[['percent']] = percent_variation
    return_list[['attributes']] = pca_attributes

    return(return_list)
}

tse_rel <- transformAssay(
    tse,
    assay.type = "counts",
    method = "relabundance"
)

rel_mat <- t(
    assay(tse_rel, "relabundance")
)

stats_list[['BrayCurtis']][['Dist']]  <- vegdist(
    rel_mat,
    method = "bray"
)


count_mat <- t(
    assay(tse, "counts")
)

# get jaccard distance matrix
stats_list[['Jaccard']][['Dist']] <- vegdist(
    count_mat,
    method = "jaccard",
    binary = TRUE
)


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


stats_list[['Aitchison']][['Dist']] <- dist(
    clr_mat,
    method = "euclidean"
)

# run PCA
tse_clr <- runPCA(
    tse_clr,
    assay.type = "clr",
    name = "PCA",
    ncomponents = 10,
    BSPARAM = BiocSingular::ExactParam()
)



stats_list[['Aitchison']][['PCA']] <- extract_pca(
    tse1 = tse_clr,
    metadata1 = metadata,
    pca_name = "PCA"
) 
stats_list[['Aitchison']][['CLR_TSE']] = tse_clr


# plot PCA (only one, not making a function)
plots_list[['Aitchison']]  <- plotReducedDim(
        tse_clr,
        dimred = "PCA",
        colour_by = "population_group"
    ) +
        labs(
            title = title_list[['Aitchison']],
            color = "Population"
        ) +
        theme_minimal(base_size = baseSizeA) +
        theme(
            plot.title = element_text(
                face = "bold",
                hjust = 0.5, 
                size = titleSizeA
            ),
            axis.title = element_text(face = "bold", size = axisTitleSizeA  ),
            # legend.title = element_text(face = "bold"),
            legend.position = "none"
        ) +
        scale_color_manual(
            values = population_colors,
            name = "Population"
        )

plots_list[["Aitchison"]]$layers[[1]]$aes_params$size <- dotSizeA

# plots_list[['Aitchison']] 

stats_list[['Unifrac_unweighted']][['Dist']] <- calculateUnifrac(
    tse,
    assay.type = "counts",
    weighted = FALSE
)

stats_list[['Unifrac_unweighted']][['Dist']] <- as.dist(stats_list[['Unifrac_unweighted']][['Dist']])

stats_list[['Unifrac_weighted']][['Dist']] <- calculateUnifrac(
    tse,
    assay.type = "counts",
    weighted = TRUE
)

stats_list[['Unifrac_weighted']][['Dist']] <- as.dist(stats_list[['Unifrac_weighted']][['Dist']])

for (beta_div1 in beta_div_metrics) {
    print(beta_div1)

    if (beta_div1 == 'Aitchison'){ # PCA performed elsewhere
        cat('')
    } else {

        # run PCoA, keep the stats
        stats_list[[beta_div1]][['PCoA']] <- run_pcoa(distance_object = stats_list[[beta_div1]][['Dist']] , 
                                              metadata1 = metadata) 

        options(repr.plot.width = 10, repr.plot.height = 8)
        plots_list[[beta_div1]] <-  plot_pcoa(
                        pcoa_return1 = stats_list[[beta_div1]][['PCoA']],
                        metadata1 = metadata,
                        plot_title = title_list[[beta_div1]],
                        group_var = "population_group",
                        popColors = population_colors
                    ) 
        # print(plots_list[[beta_div1]])
    }

    ## PERMANOVA
    set.seed(123)

    stats_list[[beta_div1]][['PERMANOVA']] <- adonis2(
        stats_list[[beta_div1]][['Dist']] ~ population_group,
        data = metadata,
        permutations = 9999
    )

    ### PERMDISP
    stats_list[[beta_div1]][["Dispersion"]] <- betadisper(
        stats_list[[beta_div1]][['Dist']],
        metadata$population_group
    )

    stats_list[[beta_div1]][['PERMDISP']] <- permutest(
        stats_list[[beta_div1]][["Dispersion"]],
        permutations = 9999
    )
}

options(repr.plot.width = 16, repr.plot.height = 16)

p_beta <- (
    plots_list[["BrayCurtis"]] + 
    plots_list[["Jaccard"]] + 
    plots_list[["Aitchison"]] + 
    plots_list[["Unifrac_weighted"]] + 
    plots_list[["Unifrac_unweighted"]]
) +
    plot_layout(
        ncol = 2,
        guides = "collect"
    ) +
    plot_annotation(
        title = "Beta Diversity of Gut Microbiome Communities",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                hjust = 0.5,
                size = figureTitleSizeA
            )
        )
    ) &
    theme(
        axis.title = element_text(
            face = "bold",
            size = axisTitleSizeA
        ),
        axis.text = element_text(
            size = axisTextSizeA
        ),
        legend.title = element_text(
            face = "bold"
        ),
        panel.grid.minor = element_blank()
    )

p_beta


# save figure 
# main 
ggsave(
    filename = paste0(FIG_BETA_DIV_DIR, '/', PROJECT_NAME, '.beta_diversity.PCA_PCoA.png'),
    plot = p_beta,
    width = 16,
    height = 16,
    units = "in",
    dpi = 300,
    bg = "white"
)

ggsave(
    filename = paste0(FIG_BETA_DIV_DIR, '/', PROJECT_NAME, '.beta_diversity.PCA_PCoA.pdf'),
    plot = p_beta,
    width = 16,
    height = 16,
    units = "in",
    bg = "white"
)


saveRDS(plots_list, file = file.path(QUARTO_PLOT_TABLE_DIR, paste0(PROJECT_NAME, ".beta_diversity.PCA_PCoA.plotList.tsv")))




extract_permanovas <- function(
    stats_list,
    beta_div_metrics
) {

    permanova_table <- do.call(
        rbind,
        lapply(beta_div_metrics, function(metric) {

            x <- as.data.frame(
                stats_list[[metric]][["PERMANOVA"]]
            )

            x$Term <- rownames(x)
            rownames(x) <- NULL

            x$Metric <- metric

            # Put identifying columns first
            x <- x[, c(
                "Metric",
                "Term",
                setdiff(
                    colnames(x),
                    c("Metric", "Term")
                )
            )]

            return(x)
        })
    )

    rownames(permanova_table) <- NULL

    return(permanova_table)
}


extract_permanova_results <- function(
    stats_list,
    beta_div_metrics
) {

    results <- do.call(
        rbind,
        lapply(beta_div_metrics, function(metric) {

            x <- as.data.frame(
                stats_list[[metric]][["PERMANOVA"]]
            )

            data.frame(
                Metric = metric,
                Df = x$Df[1],
                SumOfSqs = x$SumOfSqs[1],
                R2 = x$R2[1],
                F = x$F[1],
                p_value = x$`Pr(>F)`[1]
            )
        })
    )

    rownames(results) <- NULL

    results$p_adj_BH <- p.adjust(
        results$p_value,
        method = "BH"
    )

    return(results)
}

permanova_table <- extract_permanovas(
    stats_list = stats_list,
    beta_div_metrics = beta_div_metrics
)

permanova_table


write_tsv( 
    permanova_table,
    file = file.path(SUMMARY_TABLE_DIR, paste0(PROJECT_NAME, ".beta_diversity.permanova.stats.tsv")),
    quote = "needed"
)


permanova_results <- extract_permanova_results(
    stats_list,
    beta_div_metrics
)

permanova_results

# write to quarto folder so i can make figure for main
write_tsv( 
    permanova_results,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".beta_diversity.permanova.statsShort.tsv")),
    quote = "needed"
)


extract_permdisp_results <- function(
    stats_list,
    beta_div_metrics
) {
    results <- do.call(
        rbind,
        lapply(beta_div_metrics, function(metric) {
            x <- stats_list[[metric]][["PERMDISP"]]
            data.frame(
                Metric = metric,
                Df = x$tab[1, "Df"],
                F = x$tab[1, "F"],
                p_value = x$tab[1, "Pr(>F)"]
            )
        })
    )
    rownames(results) <- NULL
    # mult test correction across distance metrics
    results$p_adj_BH <- p.adjust(
        results$p_value,
        method = "BH"
    )

    return(results)
}

permdisp_results <- extract_permdisp_results(
    stats_list = stats_list,
    beta_div_metrics = beta_div_metrics
)

permdisp_results


write_tsv( 
    permdisp_results,
    file = file.path(SUMMARY_TABLE_DIR, paste0(PROJECT_NAME, ".beta_diversity.permdisp.stats.tsv")),
    quote = "needed"
)




extract_ordination_variance <- function(
    stats_list,
    beta_div_metrics
) {

    results <- do.call(
        rbind,
        lapply(beta_div_metrics, function(metric) {

            if (metric == "Aitchison") {

                # PCA variance explained
                pca <- stats_list[[metric]][["PCA"]]

                data.frame(
                    Metric = metric,
                    Ordination = "PCA",
                    Axis1 = pca$percent[1],
                    Axis2 = pca$percent[2],
                    Axis1_2 = sum(pca$percent[1:2])
                )

            } else {

                # PCoA variance explained
                pcoa <- stats_list[[metric]][["PCoA"]]

                data.frame(
                    Metric = metric,
                    Ordination = "PCoA",
                    Axis1 = pcoa$percent[1],
                    Axis2 = pcoa$percent[2],
                    Axis1_2 = sum(pcoa$percent[1:2])
                )
            }
        })
    )

    rownames(results) <- NULL

    return(results)
}

ordination_variance <- extract_ordination_variance(
    stats_list = stats_list,
    beta_div_metrics = beta_div_metrics
)

ordination_variance


write_tsv( 
    ordination_variance,
    file = file.path(SUMMARY_TABLE_DIR, paste0(PROJECT_NAME, ".beta_diversity.PCoA.ordination_variance.tsv")),
    quote = "needed"
)




p_permanova <- ggplot(
    permanova_results,
    aes(
        x = reorder(Metric, R2),
        y = R2
    )
) +
    geom_col(
        width = 0.7
    ) +
    geom_text(
        aes(
            label = paste0(
                "R² = ",
                sprintf("%.3f", R2)
            )
        ),
        hjust = -0.1,
        size = 4
    ) +
    coord_flip() +
    labs(
        title = "Effect of Population on Microbial Community Composition",
        x = NULL,
        y = expression(PERMANOVA~R^2)
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(
            face = "bold",
            hjust = 0.5
        ),
        axis.title = element_text(
            face = "bold"
        ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    expand_limits(
        y = max(permanova_results$R2) * 1.2
    )

p_permanova

# permanova_results$significance <- cut(
#     permanova_results$p_adj_BH,
#     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
#     labels = c("***", "**", "*", "ns")
# )

# geom_text(
#     aes(
#         label = paste0(
#             "R² = ", sprintf("%.3f", R2),
#             "  ", significance
#         )
#     ),
#     hjust = -0.1,
#     size = 4
# )



extract_dispersion_distances <- function(
    stats_list,
    beta_div_metrics,
    metadata,
    group_var = "population_group"
) {

    results <- do.call(
        rbind,
        lapply(beta_div_metrics, function(metric) {

            dispersion <- stats_list[[metric]][["Dispersion"]]

            sample_ids <- names(dispersion$distances)

            data.frame(
                Metric = metric,
                Sample = sample_ids,
                Population = metadata[
                    sample_ids,
                    group_var
                ],
                Distance = as.numeric(
                    dispersion$distances
                ),
                stringsAsFactors = FALSE
            )
        })
    )

    rownames(results) <- NULL

    return(results)
}

dispersion_distances <- extract_dispersion_distances(
    stats_list = stats_list,
    beta_div_metrics = beta_div_metrics,
    metadata = metadata,
    group_var = "population_group"
)

head(dispersion_distances)


write_tsv( 
    dispersion_distances,
    file = file.path(SUMMARY_TABLE_DIR, paste0(PROJECT_NAME, ".beta_diversity.Dispersion.Dists.tsv")),
    quote = "needed"
)


permdisp_annot <- permdisp_results %>%
    dplyr::mutate(
        p_label = dplyr::case_when(
            p_adj_BH < 0.001 ~ "***",
            p_adj_BH < 0.01  ~ "**",
            p_adj_BH < 0.05  ~ "*",
            TRUE             ~ "ns"
        )
    ) %>%
    dplyr::left_join(
        dispersion_distances %>%
            dplyr::group_by(Metric) %>%
            dplyr::summarise(
                y.position = max(Distance, na.rm = TRUE) * 1.10,
                .groups = "drop"
            ),
        by = "Metric"
    ) %>%
    dplyr::mutate(
        group1 = "BF",
        group2 = "EU"
    )

p_dispersion <- ggplot(
    dispersion_distances,
    aes(
        x = Population,
        y = Distance,
        fill = Population
    )
) +
    geom_violin(trim = FALSE) +
    geom_jitter(
        width = 0.08,
        size = 2,
        alpha = 0.7
    ) +
    geom_text(
        data = permdisp_annot,
        aes(
            x = 1.5,
            y = y.position,
            label = p_label
        ),
        inherit.aes = FALSE,
        size = 5
    ) +
    facet_wrap(
        ~ Metric,
        scales = "free_y"
    ) +
    scale_fill_manual(
        values = population_colors
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.15))
    ) +
    labs(
        title = "Dispersion  Composition by Population",
        x = NULL,
        y = "Distance to Group Centroid"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        strip.text = element_text(
            size = 15,
            face = "bold"
        ),
          plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5
    ),
        axis.title = element_text(
        size = 13,
        face = "bold"
    )
    ) 

p_dispersion




# save figure 
# main 
ggsave(
    filename = paste0(FIG_BETA_DIV_DIR, '/', PROJECT_NAME, '.beta_diversity.dispersionComposition.png'),
    plot = p_dispersion,
    width = 16,
    height = 16,
    units = "in",
    dpi = 300,
    bg = "white"
)

ggsave(
    filename = paste0(FIG_BETA_DIV_DIR, '/', PROJECT_NAME, '.beta_diversity.dispersionComposition.pdf'),
    plot = p_dispersion,
    width = 16,
    height = 16,
    units = "in",
    bg = "white"
)





write_tsv( 
    permdisp_annot,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".beta_diversity.PermDist_byPop.tsv")),
    quote = "needed"
)












make_beta_stats_table <- function(
    permanova_results,
    permdisp_results
) {

    beta_stats <- permanova_results %>%
        dplyr::select(
            Metric,
            PERMANOVA_R2 = R2,
            PERMANOVA_F = F,
            PERMANOVA_p = p_value,
            PERMANOVA_p_adj_BH = p_adj_BH
        ) %>%
        dplyr::left_join(
            permdisp_results %>%
                dplyr::select(
                    Metric,
                    PERMDISP_F = F,
                    PERMDISP_p = p_value,
                    PERMDISP_p_adj_BH = p_adj_BH
                ),
            by = "Metric"
        ) %>%
        dplyr::mutate(
            Metric = dplyr::recode(
                Metric,
                "BrayCurtis" = "Bray–Curtis",
                "Jaccard" = "Jaccard",
                "Aitchison" = "Aitchison",
                "Unifrac_weighted" = "Weighted UniFrac",
                "Unifrac_unweighted" = "Unweighted UniFrac"
            ),
            dplyr::across(
                c(PERMANOVA_R2, PERMANOVA_F, PERMDISP_F),
                ~ round(.x, 3)
            )
        )

    return(beta_stats)
}

beta_stats_table <- make_beta_stats_table(
    permanova_results = permanova_results,
    permdisp_results = permdisp_results
)

beta_stats_table



write_tsv( 
    beta_stats_table,
    file = file.path(REPORT_TABLE_DIR, paste0(PROJECT_NAME, ".beta_diversity.CombinedStatistics_PERM.tsv")),
    quote = "needed"
)
















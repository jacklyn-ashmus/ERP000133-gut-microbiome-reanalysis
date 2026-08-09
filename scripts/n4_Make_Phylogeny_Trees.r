suppressPackageStartupMessages({

    ########### main libraries ###########
    library(SummarizedExperiment)
    library(TreeSummarizedExperiment)
    library(mia)
    library(miaViz) # mia visualizatioon
    library(vegan) # ecol stats?

    
    ########### data manipulation tools ###########
    library(tidyverse)
    
    ########### visualization libraries ###########   
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


# read rds file
tse <- readRDS(tse_rds_file)


# tse

# make file output directory if it doesn't already exist
ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
PHYLO_DIR <- file.path(FIG_DIR, "phylo")
FIG_PHYLO_DIR <- file.path(FIG_DIR, "phylo")

dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PHYLO_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_PHYLO_DIR, recursive = TRUE, showWarnings = FALSE)

tse_genus <- agglomerateByRank(
    tse,
    rank = "Genus",
    update.tree = FALSE
)

tse_genus

tse_genus <- addHierarchyTree(tse_genus)

# getTaxonomyRanks()

# test plot row tree 
plotRowTree(tse_genus)

# top 50 features
top <- getTop(tse_genus, top = 50L)

# Plot
plotRowTree(
    tse_genus[top, ],
    tip.colour.by = "Phylum"
)

# Calculate mean abundance
tse_genus <- addPerFeatureQC(tse_genus)
rowData(tse_genus)[["log_mean"]] <- log(rowData(tse_genus)[["mean"]])

# basic dendrogram
# # Plot
# plotRowTree(
#     tse_genus[top, ],
#     layout = "dendrogram",
#     edge.colour.by = "Phylum",
#     tip.colour.by = "log_mean"
# )


options(repr.plot.width = 18, repr.plot.height = 13)

p <- plotRowTree(
    tse_genus[top, ],
    layout = "rectangular",
    edge.colour.by = "Phylum",
    tip.colour.by = "log_mean",
    show.label = TRUE,
    relabel.tree = TRUE,
    levels.rm = TRUE
)

# Increase the interactive text layer sizes
p$layers[[5]]$aes_params$size <- 5 # changes dots 
p$layers[[6]]$aes_params$size <- 5 # changes leaves
p$layers[[7]]$aes_params$size <- 5 # changes branches
p$layers[[7]]$position <- position_nudge(y = 0.5, x = -0.058)

p <- p +
    coord_cartesian(clip = "off") +
    theme(
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        legend.box.margin = margin(
            t = 10,
            r = 0,
            b = 0,
            l = 50 # space on LEFT of legend = pushes legend right
        ) ,
        # Actual space outside panel for labels to spill into
        plot.margin = margin(
            t = 10,
            r = 250,
            b = 10,
            l = 30
        )
    )

p

ggsave(
    filename = paste0(FIG_PHYLO_DIR, '/', 'genus_taxonomy_tree.png'),
    plot = p,
    width = 16.4,
    height = 11,
    units = "in",
    dpi = 300,
    bg = "white"
)

ggsave(
    filename = paste0(FIG_PHYLO_DIR, '/', 'genus_taxonomy_tree.pdf'),
    plot = p,
    width = 15,
    height = 11,
    units = "in",
    bg = "white"
)



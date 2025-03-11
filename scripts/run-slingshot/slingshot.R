library(slingshot)
library(rjson)
library(magrittr)
source("scripts/run-slingshot/utils.R")

# # these are the arguments that are passed to the script
# background_embedding <- "outputs/seurat_integration/groups/four_trees_batch1/rs17/umap2d/projections/rs17.csv"
# umap_projection <- "outputs/seurat_integration/groups/four_trees_batch1/rs17/umap2d/projections/rs17_egr_nextseq_batch1_all.csv"
# cell_identity <- "outputs/seurat_integration/groups/four_trees_batch1/rs17/cell_identities_2d_17_egr_nextseq_batch1_all.csv"
# # lineages <- "3,4,8/6,5,2,1/6,5,2,7"
# lineages <- "3,4,8/6,5,2,7/6,5,2,1"
# palette <- "configs/palette_all.json"
# thresh <- 0.001
# output_pdf <- "outputs/seurat_integration/groups/four_trees_batch1/rs17/umap2d/slingshot/rs17_egr_nextseq_batch1_all.pdf"
# output_dir <- "outputs/seurat_integration/groups/four_trees_batch1/rs17/umap2d/slingshot/rs17_egr_nextseq_batch1_all/"

# parse command line arguments and load data
parse_arguments()
background_embedding <- read.csv(background_embedding, row.names = 1) # Barcode, UMAP.1, UMAP.2
umap_projection <- read.csv(umap_projection, row.names = 1) # Barcode, UMAP.1, UMAP.2
cell_identity <-  read.csv(cell_identity, row.names = 1) # Barcode, Cluster, Color
lineages_list <- strsplit(lineages, "/")[[1]] %>% lapply(function(x) strsplit(x, ",")[[1]])
palette <- fromJSON(paste(readLines(palette), collapse=""))
ratio <- as.numeric(ratio)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
} else {
    message("Warning: Output directory already exists. Deleting and recreating.")
    unlink(output_dir, recursive = TRUE)
    dir.create(output_dir, recursive = TRUE)
}

# lineages are cluster color codes separated by ','
# multiple lineages are separated by '/'

# create cluster label matrix
# cell_identity$Color <- sapply(cell_identity$Cluster, function(x) palette[[x]]) %>% unlist
clusters <- cell_identity[rownames(umap_projection), "Cluster"] %>% factor
clusterLabels <- sapply(
    as.numeric(clusters),
    function(x, max_level) {
        out <- rep(0, max_level)
        out[x] <- 1
        return(out)
    },
    max_level = max(as.numeric(clusters))
) %>% t
rownames(clusterLabels) <- rownames(umap_projection)
colnames(clusterLabels) <- levels(clusters)

# plot slingshot
# open device based on suffix
width <- 7
height <- width * ratio
if (grepl(".pdf", output_pdf)) {
    pdf(output_pdf, width = width, height = height)
} else if (grepl(".svg", output_pdf)) {
    svg(output_pdf, width = width, height = height)
} else if (grepl(".eps", output_pdf)) {
    postscript(output_pdf, width = width, height = height)
}
# plot dots
par(mar = c(0, 0, 0, 0))
par(mai = c(0, 0, 0, 0))
xlim <- c(min(background_embedding[, 1]), max(background_embedding[, 1]))
ylim <- c(min(background_embedding[, 2]), max(background_embedding[, 2]))
plot(
    umap_projection[, 1],
    umap_projection[, 2],
    col = "#D4D4D4",
    pch = 20,
    cex = 0.3,
    ann = FALSE,
    xaxt = "n", 
    yaxt = "n",
    xlim = xlim,
    ylim = ylim,
    axes = FALSE
)

names(lineages_list) <- paste0("lineage_", seq(length(lineages_list)))
for (n in seq_len(length(lineages_list))) {
    lineages <- lineages_list[n]
    # build adjacency matrix
    adjacency = matrix(0,max(as.numeric(clusters)),max(as.numeric(clusters)))
    rownames(adjacency) = levels(clusters)
    colnames(adjacency) = levels(clusters)
    for(L in lineages){
        for(Ci in seq(length(L)-1)){
            adjacency[L[Ci],L[Ci+1]] = 1
            adjacency[L[Ci+1],L[Ci]] = 1
        }
    }

    # get reduced dimension
    reducedDim <- as.matrix(umap_projection[, c("UMAP.1", "UMAP.2")])

    # run slingshot
    lineage <- new("SlingshotDataSet")
    lineage@reducedDim <- reducedDim
    lineage@clusterLabels <- clusterLabels
    lineage@lineages <- lineages
    lineage@adjacency <- adjacency
    curves <- getCurves(
        lineage,
        extend = "n",
        thresh = 0.05,
        smoother = "smooth.spline",
        stretch = 0,
        approx_points = 600
    ) %>% SlingshotDataSet()

    # current best
    # extend = "n",
    # thresh = 0.05,
    # smoother = "smooth.spline",
    # stretch = 0


    # plot lineages
    for (lineage_n in names(slingLineages(curves))) {
        curve_points <- slingCurves(curves)[[lineage_n]]$s
        lineage <- slingLineages(curves)[[lineage_n]] # named character vector
        lineage_color <- palette[lineage]
        lineage_cluster_centers <- get_cluster_centers(
            projections = reducedDims(curves),
            clusters = clusters,
            target_clusters = lineage
        )
        col <- get_curve_color(
            curve_points = curve_points,
            lineage = lineage,
            lineage_color = lineage_color,
            lineage_cluster_centers = lineage_cluster_centers,
            smooth_factor=.95
        )
        points(
            curve_points,
            col = col,
            pch = 20,
            cex = 2
        )
    }
    saveRDS(list(
        curve_points = curve_points,
        color = col
    ), file = file.path(output_dir, paste0("lineage_", n, ".rds")))
}
# get lineage 
# 
dev.off()
library(magrittr)

parse_arguments <- function(args) {
    # parse command line arguments
    args <- R.utils::commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]
    # write to global environment
    keys <- names(args)
    for (key in keys) {
        assign(key, args[[key]], envir = globalenv())
    }
    message(
        "Assigned global variables:",
        str(mget(keys, envir = globalenv()))
    )
}

get_curve_color <- function(
    curve_points,
    lineage,
    lineage_color,
    lineage_cluster_centers,
    smooth_factor
) {
    # curve_points: points composing the principle curve estimated by Slingshot
    # lineage: list of clusters in the lineage (character vector)
    # lineage_color: the color for each cluster in the lineage
    # lineage_cluster_center: the center of each cluster in the lineage
    # smooth_factor: the smooth factor (between 0-1)

    # get the closest cluster center to each point in the curve
    distances <- rbind(
        lineage_cluster_centers,
        curve_points
    ) %>% dist %>% as.matrix
    dist2centers <- distances[
        (nrow(lineage_cluster_centers) + 1):nrow(distances),
        1:nrow(lineage_cluster_centers)
    ]
    closest_center <- apply(
        dist2centers,
        1,
        function(row) {
            which.min(row) %>% names
        }
    )
    colors <- sapply(closest_center, function(name) {
        lineage_color[[name]]
    })

    # smooth out colors
    color_function <- colorRampPalette(
        colors[seq(
            from = 1,
            to = length(colors), 
            length.out = length(colors) * (1 - smooth_factor)
        )]
    )
    return(color_function(length(colors)))   
}

get_cluster_centers <- function(
    projections,
    clusters,
    target_clusters
) {
    # projections: the projections of each cell
    # clusters: the cluster of each cell
    # target_clusters: the clusters to calculate centers for
    x <- sapply(target_clusters, function(cluster_n) {
        mean(projections[clusters == cluster_n, 1])
    })
    y <- sapply(target_clusters, function(cluster_n) {
        mean(projections[clusters == cluster_n, 2])
    })
    return(cbind(x, y))
}

# run_slingshot <- function(
#     reducedDim,
#     clusterLabels,
#     lineages,
#     adjacency
# ) {
#     # create a new SlingshotDataSet object
#     lineage <- new("SlingshotDataSet")
#     lineage@reducedDim <- reducedDim
#     lineage@clusterLabels <- clusterLabels
#     lineage@lineages <- lineages
#     lineage@adjacency <- adjacency

#     # run the slingshot algorithm
#     lineage <- getCurves(lineage, extend="n")

#     return(lineage)
# }

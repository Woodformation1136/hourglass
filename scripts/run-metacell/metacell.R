# import
library(metacell)
library(ggplot2)
library(magrittr)
library(dplyr) # to join color information to
library(rjson)
library(Matrix)
# source("scripts/run-metacell/utils.R", chdir = TRUE)
source("scripts/run-metacell/utils.R", chdir = TRUE)

# configuration
work_dir <- getwd()

# parse arguments
parse_arguments()

# TODO: remove the following variables that should be set by parse_arguments
# matrix_path <- "outputs/seurat_integration/groups/four_trees/rs17/corrected_expression_matrices/ptr_tenx_batch1_rs17.csv"
# cell_ident <- "outputs/seurat_integration/groups/four_trees/rs17/cell_identities.csv"
# palette <- "configs/palette.json"
# scdb_tmpdir <- "outputs/seurat_integration/groups/four_trees/rs17/metacell/ptr_tenx_batch1_tmp"
# output_pdf <- "outputs/seurat_integration/groups/four_trees/rs17/metacell/ptr_tenx_batch1.pdf"
# scm_n_downsamp_gstat <- "300"

# create tmp directory (clear if exists)
if (dir.exists(scdb_tmpdir)) {
    message("Removing existing tmp directory", scdb_tmpdir)
    unlink(scdb_tmpdir, recursive = TRUE)
}
dir.create(scdb_tmpdir, recursive = TRUE)
if (dir.exists(scdb_tmpdir)) {
    message("Created tmp directory ", scdb_tmpdir)
} else {
    stop("Failed to create tmp directory ", scdb_tmpdir)
}
on.exit({
    message("Exit...Removing tmp directory ", scdb_tmpdir)
    unlink(scdb_tmpdir, recursive = TRUE)
})

# manually setting scm_n_downsamp_gstat
message("Setting scm_n_downsamp_gstat to ", scm_n_downsamp_gstat)
tgconfig::set_param(
    "scm_n_downsamp_gstat",
    as.integer(scm_n_downsamp_gstat),
    "metacell"
)
if (tgconfig::get_param("scm_n_downsamp_gstat", "metacell") != scm_n_downsamp_gstat) {
    stop("Failed to set scm_n_downsamp_gstat to ", scm_n_downsamp_gstat)
} else {
    message("Set scm_n_downsamp_gstat to ", scm_n_downsamp_gstat)
}


# Prepare for metacell analyses
# -------------------------------------------------------------
import.umi(
    scdb_tmpdir=scdb_tmpdir,
    matrix_path=matrix_path
)
# -------------------------------------------------------------

# MetaCell analysis
# -------------------------------------------------------------
custom_metacell_pipeline(
    mat_id="mat",
    gstat_id="gstat",
    gset_id="gset",
    graph_id="graph",
    coc_id="coc",
    mc_id="mc",
    mc2d_id="mc2d",
    scm_n_downsamp_gstat=scm_n_downsamp_gstat
)
# -------------------------------------------------------------

# Visualization
# -------------------------------------------------------------
# get the 2d projection
mc2d <- scdb_mc2d("mc2d")

# get single cell projections
projections <- data.frame(
    "X" = mc2d@sc_x,
    "Y" = mc2d@sc_y
)
projections$Barcode <- rownames(projections)

# get metacell projections
mc_projections <- data.frame(
    "X" = mc2d@mc_x,
    "Y" = mc2d@mc_y
)

# give cell identity
cell_ident <- read.csv(cell_ident, header = T, stringsAsFactors = F)
projections %<>% dplyr::inner_join(cell_ident, by = c("Barcode" = "Barcode"))
projections$Cluster <- as.factor(projections$Cluster)

# make palette
if (file.exists(palette)) {
    palette <- fromJSON(paste(readLines(palette), collapse=""))
} else {
    c2c_mapping <- projections[, c("Cluster", "Color")] %>% unique
    palette <- c2c_mapping$Color
    names(palette) <- c2c_mapping$Cluster %>% as.character
}

# # get edges (use NA to break the line)
# fr <- mc2d@graph$mc1
# to <- mc2d@graph$mc2
# # NOTE: plot only metacell number: 1,2,3,4,5,6,7,8,9,10,18,19,20
# edges <- data.frame(
#     X_start = mc2d@mc_x[fr],
#     Y_start = mc2d@mc_y[fr],
#     X_end = mc2d@mc_x[to],
#     Y_end = mc2d@mc_y[to]
# )
# edges[fr %in% c(1,2,3,4,5,6,7,8,9,10,18,19,20) & to %in% c(1,2,3,4,5,6,7,8,9,10,18,19,20),] %<>% 
#     dplyr::mutate(
#         X_end = NA,
#         Y_end = NA
#     )
# write.csv(mc2d@graph, "test/metacell/graph.csv", row.names = F)

# plot
if (!dir.exists(dirname(output_pdf))) {
    dir.create(dirname(output_pdf))
}
pdf(file = output_pdf)
ggplot(
    projections,
    aes(X,Y)
) + geom_point(
        data = projections,
        aes(X, Y, color = Cluster),
        size = 0.4,
        show.legend = FALSE
    ) + scale_color_manual(
        values = palette
    ) + coord_fixed(
        ratio = 1
    )

# the following lines add metacell and edges to the plot
# +
# geom_label(
#     data = mc_projections,
#     aes(X, Y, label = rownames(mc_projections)),
#     show.legend = FALSE,
#     fill = NA,
#     size = 5
# ) +
# geom_segment(
#     data = edges,
#     aes(
#         x = X_start,
#         y = Y_start,
#         xend = X_end,
#         yend = Y_end
#     ),
#     linewidth = 0.1
# )

# close device
dev.off()

# clean up
unlink(scdb_tmpdir, recursive = TRUE)

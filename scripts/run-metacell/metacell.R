# import
library(metacell)
library(ggplot2)
library(magrittr)
library(dplyr) # to join color information to
library(rjson)
library(Matrix)
# source("scripts/run-metacell/utils.R", chdir=TRUE)
source("scripts/run-metacell/utils.R", chdir=TRUE)

# parse arguments
parse_arguments()

# TODO: remove the following variables that should be set by parse_arguments
# matrix_path <- "outputs/seurat_integration/groups/four_trees/rs17/corrected_expression_matrices/ptr_tenx_batch1_rs17.csv"
# cell_ident <- "outputs/seurat_integration/groups/four_trees/rs17/cell_identities.csv"
# palette <- "configs/palette.json"
# scdb_tmpdir <- "outputs/seurat_integration/groups/four_trees/rs17/metacell/ptr_tenx_batch1_tmp"
# output_pdf <- "outputs/seurat_integration/groups/four_trees/rs17/metacell/ptr_tenx_batch1.pdf"
# scm_n_downsamp_gstat <- "300"
# gset_out <- "outputs/seurat_integration/groups/four_trees/rs17/metacell/ptr_tenx_batch1_gset.rds"

# matrix_path <- "outputs/seurat_integration/groups/four_trees/rs17/corrected_expression_matrices/egr_nextseq_batch1_rs17.csv"
# cell_ident <- "outputs/seurat_integration/groups/four_trees/rs17/cell_identities.csv"
# palette <- "configs/palette.json"
# scdb_tmpdir <- "outputs/seurat_integration/groups/four_trees/rs17/metacell/egr_nextseq_batch1_tmp"
# output_pdf <- "outputs/seurat_integration/groups/four_trees/rs17/metacell/egr_nextseq_batch1.pdf"
# scm_n_downsamp_gstat <- "30"
# gset_in <- "outputs/seurat_integration/groups/four_trees/rs17/metacell/ptr_tenx_batch1_gset.rds"
# T_tot=20
# T_top3=1
# T_szcor=-0.01
# T_vm=0.2
# T_niche=0.05
# K=50
# min_mc_size=50
# alpha=10

# create tmp directory (clear if exists)
if (dir.exists(scdb_tmpdir)) {
    message("Removing existing tmp directory", scdb_tmpdir)
    unlink(scdb_tmpdir, recursive=TRUE)
}
dir.create(scdb_tmpdir, recursive=TRUE)
if (dir.exists(scdb_tmpdir)) {
    message("Created tmp directory ", scdb_tmpdir)
} else {
    stop("Failed to create tmp directory ", scdb_tmpdir)
}


# # manually setting scm_n_downsamp_gstat
# message("Setting scm_n_downsamp_gstat to ", scm_n_downsamp_gstat)
# tgconfig::set_param(
#     "scm_n_downsamp_gstat",
#     as.integer(scm_n_downsamp_gstat),
#     "metacell"
# )
# if (tgconfig::get_param("scm_n_downsamp_gstat", "metacell") != scm_n_downsamp_gstat) {
#     stop("Failed to set scm_n_downsamp_gstat to ", scm_n_downsamp_gstat)
# } else {
#     message("Set scm_n_downsamp_gstat to ", scm_n_downsamp_gstat)
# }

# Prepare for metacell analyses
# -------------------------------------------------------------
success <- import.umi(
    scdb_tmpdir=scdb_tmpdir,
    matrix_path=matrix_path
)
if (success) {
    message("Imported UMI matrix")
} else {
    stop("Failed to import UMI matrix")
}
print(scdb_mat("mat"))
# -------------------------------------------------------------

# MetaCell analysis
# -------------------------------------------------------------
# if (exists("gset_in")) {
#     message("Loading gene set from ", gset_in)
#     gset <- readRDS(gset_in)
#     scdb_add_gset("gset", gset)
# }

# run
custom_metacell_pipeline(
    mat_id="mat",
    gstat_id="gstat",
    gset_id="gset",
    graph_id="graph",
    coc_id="coc",
    mc_id="mc",
    mc2d_id="mc2d",
    T_tot=as.numeric(T_tot),
    T_top3=as.numeric(T_top3),
    T_szcor=as.numeric(T_szcor),
    T_vm=as.numeric(T_vm),
    T_niche=as.numeric(T_niche),
    K=as.numeric(K),
    min_mc_size=as.numeric(min_mc_size),
    alpha=as.numeric(alpha)
)
if (exists("gset_out")) {
    message("Saving gene set to ", gset_out)
    gset <- scdb_gset("gset")
    saveRDS(gset, gset_out)
}
# T_tot: total down sampled coverage thresholds 
# (genes with tot UMIs < T_tot are filtered out)

# T_top3: threshold value for the third highest umi count for the gene 
# (genes with top3<T_top3 are filtered out)

# T_szcor: threshold value for the normalized size correlation statistic 
# (only genes with sz_cor < T_szcor are selected). If you use this, consider values around -0.1 - but evaluate carefully your decision using the gstat empirical data

# T_vm: the threshold value for the normalized var/mean 
# (only genes with varmean > T_vm are selected) Recommended values are usually around 0.2, but this may vary with the data. Not recommended for datasets with hihgly heterogeneous cell sizes (e.g. in whole-organisms datasets)

# T_niche: threshold value for the normalized niche score statistic 
# (only genes with niche_norm > T_niche are selected). Recommended to use in combination with szcor to add genes with strongly restricted expression patterns. Consider using values around 0.05

# K: the number of top-coclustering neighbors to consider per node

# alpha: the relexation paramter to apply for filtering coclustering neighbors below the top K ones. 
# A pair (n1,n2) with weight w will be filtered if knn(n1,K) > w*alpha or knn(n2,K) > w*alpha
# -------------------------------------------------------------

# Visualization
# -------------------------------------------------------------
# plot
p <- plot_metacell(
    mc2d_id="mc2d",
    cell_ident=read.csv(cell_ident, header=T, stringsAsFactors=F),
    palette=fromJSON(paste(readLines(palette), collapse="")),
)
# record params to figure title

# save plot
if (!dir.exists(dirname(output_pdf))) {
    dir.create(dirname(output_pdf))
}
if (file.exists(output_pdf)) {
    removed <- file.remove(output_pdf)
    if (removed) message("Removed existing file ", output_pdf)
}
pdf(file=output_pdf)
print(p)
dev.off()
# -------------------------------------------------------------

# Clean up
# -------------------------------------------------------------
unlink(scdb_tmpdir, recursive=TRUE)
# -------------------------------------------------------------

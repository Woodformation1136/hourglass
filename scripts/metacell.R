# import
library(metacell)
library(ggplot2)
library(magrittr)
library(dplyr) # to join color information to 
library(Matrix)
 

# configuration
work_dir <- getwd()

# get snakemake variables
# matrix_dir <- snakemake@input$matrix_dir
# barcode2color_csv <- snakemake@input$barcode2color_csv
# scdb_tmpdir <- snakemake@params$scdb_tmpdir
# output_pdf <- snakemake@output$pdf


sample <- "egr_nextseq_batch1_rs17"
# sample <- "tung_batch1_rs17"
# sample <- "lch_tenx_batch1_rs17"
# sample <- "tar_nextseq_batch1_rs17"

# raw UMI count
# matrix_dir <- sprintf("outputs/cellranger/reanalyze/%s/outs/filtered_feature_bc_matrix", sample)
# output_pdf <- sprintf("test/metacell/%s_raw.pdf", sample)

# corrected orthogroup gene expression count
matrix_dir <- sprintf("outputs/seurat_integration/group1/rs17/corrected_expression_matrices/%s.csv", sample)
output_pdf <- sprintf("test/metacell/%s_seurat_ortho.pdf", sample)

# un-corrected orthogroup gene expression count
# matrix_dir <- sprintf("outputs/seurat_integration/ortho_umi_count/%s.csv", sample)
# output_pdf <- sprintf("test/metacell/%s_uc_ortho.pdf", sample)

barcode2color_csv <- sprintf("test/barcode2color/%s.csv", sample)
scdb_tmpdir <- sprintf("metacell_tmp_%s", sample)

# Prepare for metacell analyses
# -------------------------------------------------------------
# prepare a directory for scdb
if (!dir.exists(scdb_tmpdir)){
    dir.create(scdb_tmpdir)
} else {
    unlink(scdb_tmpdir, recursive = TRUE)
    dir.create(scdb_tmpdir)
} 
scdb_init(scdb_tmpdir, force_reinit = T)
# read the UMI table from 10X format
# mat <- scmat_read_scmat_10x(
#     matrix_fn = file.path(matrix_dir, "matrix.mtx.gz"),  
#     genes_fn = file.path(matrix_dir, "features.tsv.gz"),
#     cells_fn = file.path(matrix_dir, "barcodes.tsv.gz"),
#     paralogs_policy = "remove"
# )
# import the matrix into scdb
# scdb_add_mat("mat", mat)

# read the UMI table from csv format
read.csv(
    matrix_dir,
    header = FALSE,
    stringsAsFactors = FALSE,
    row.names = NULL,
) %>% t %>% 
write.table(
    file=sprintf("%s.T", matrix_dir),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
)
mcell_import_scmat_tsv(
    mat_nm = "mat",
    fn = sprintf("%s.T", matrix_dir),
    dset_nm = "dset"
)
unlink(sprintf("%s.T", matrix_dir))
# -------------------------------------------------------------

# MetaCell analysis
# -------------------------------------------------------------
# select gene markers for metacell grouping
mcell_add_gene_stat(
    gstat_id = "gene_stat",
    mat_id = "mat",
    force = T
) # compute gene stats
mcell_gset_filter_multi(
    gstat_id = "gene_stat",
    gset_id = "clustering_markers",
    T_tot = 20,
    T_top3 = 1,
    T_szcor = -0.01,
    T_vm = 0.2,
    T_niche = 0.05,
    force_new = T
) # select gene markers
# -------------------------------------------------------------
# Manually set Seurat corrected high-variabel orthologs as clustering markers
# gene_set_manual <- rep(c(1), length(scdb_mat("mat")@genes))
# names(gene_set_manual) <- scdb_mat("mat")@genes
# scdb_add_gset(
#     id="clustering_markers",
#     gset=gset_new_gset(
#         sets=gene_set_manual,
#         desc="clustering_markers"
#     )
# )
# -------------------------------------------------------------
# group cells
mcell_add_cgraph_from_mat_bknn(
    mat_id = "mat",
    gset_id = "clustering_markers",
    graph_id = "graph",
    K = 40,
    dsamp = F
) # build graph
mcell_coclust_from_graph_resamp(
    coc_id = "coc",
    graph_id = "graph",
    min_mc_size = 50,
    p_resamp = 0.75,
    n_resamp = 1000
) # cocluster by bootstraping
source("scripts/MY_mcell_mc_from_coclust_balanced.R")
MY_mcell_mc_from_coclust_balanced(
    mat_id = "mat",
    coc_id = "coc",
    mc_id = "mc",
    K = 40,
    min_mc_size = 50,
    alpha = 3
) # final grpah
# -------------------------------------------------------------
# 3. 2D projection of the k-nn graph
# force-directed projection of the k-nn graph
mcell_mc2d_force_knn(
    mc2d_id = "mc2d",
    mc_id = "mc",
    graph_id = "graph"
)
# get the 2d projection
mc2d <- scdb_mc2d("mc2d")
projections <- data.frame(
    "X" = mc2d@sc_x,
    "Y" = mc2d@sc_y
)
projections$Barcode <- rownames(projections)
mc_projections <- data.frame(
    "X" = mc2d@mc_x,
    "Y" = mc2d@mc_y
)
# get barcode to color mapping
barcode2color <- read.csv(barcode2color_csv, header = T, stringsAsFactors = F)
# merge color information to projections
projections %<>% dplyr::inner_join(barcode2color, by = c("Barcode" = "Barcode"))
projections$Cluster <- as.factor(projections$Cluster)
# create cluster to color mapping
c2c_mapping <- projections[, c("Cluster", "Color")] %>% unique
palette <- c2c_mapping$Color
names(palette) <- c2c_mapping$Cluster %>% as.character
# get edges (use NA to break the line)
fr <- mc2d@graph$mc1
to <- mc2d@graph$mc2
edges <- data.frame(
    X_start = mc2d@mc_x[fr],
    Y_start = mc2d@mc_y[fr],
    X_end = mc2d@mc_x[to],
    Y_end = mc2d@mc_y[to]
)
write.csv(mc2d@graph, "test/metacell/graph.csv", row.names = F)

# plot
if (!dir.exists(dirname(output_pdf))) dir.create(dirname(output_pdf))
pdf(file = output_pdf)
ggplot(
    projections,
    aes(X,Y)
) + 
    geom_point(
        data = projections,
        aes(X, Y, color = Cluster),
        size = 0.4,
        show.legend = FALSE
    ) + 
    scale_color_manual(
        values = palette
    ) + 
    coord_fixed(ratio = 1) +
    geom_label(
        data = mc_projections,
        aes(X, Y, label = rownames(mc_projections)),
        show.legend = FALSE,
        fill = NA,
        size = 5
    ) +
    geom_segment(
        data = edges,
        aes(
            x = X_start,
            y = Y_start,
            xend = X_end,
            yend = Y_end
        ),
        linewidth = 0.1
    )
dev.off()

# clean up
unlink(scdb_tmpdir, recursive = TRUE)

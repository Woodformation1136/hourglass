library(metacell)
library(magrittr)
library(rjson)

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


import.umi <- function(
    scdb_tmpdir,
    matrix_path,
    mat_name="mat"
) {
    scdb_init(scdb_tmpdir, force_reinit = T)

    # if matrix_path is a directory, read 10X format
    if (file_test("-d", matrix_path)) {
        mat <- scmat_read_scmat_10x(
            matrix_fn = file.path(matrix_path, "matrix.mtx.gz"),  
            genes_fn = file.path(matrix_path, "features.tsv.gz"),
            cells_fn = file.path(matrix_path, "barcodes.tsv.gz"),
            paralogs_policy = "remove"
        )
        scdb_add_mat(mat_name, mat)
        return(TRUE)
    }

    # if matrix_path is a CSV file, read CSV format
    if (file_test("-f", matrix_path)) {
        # transpose the matrix and write to a temporary file
        # TODO: add 1 to all counts to avoid 0 counts
        # mat <- read.csv(
        #     matrix_path,
        #     header = TRUE,
        #     stringsAsFactors = FALSE,
        #     row.names = 1,
        # )
        # mat <- mat + 1
        # mat <- t(mat)
        # mat <- cbind(
        #     row.names(mat),
        #     mat
        # )
        # mat <- rbind(
        #     c("Barcode", colnames(mat)[-1]),
        #     mat
        # )
        # mat %>% write.table(
        #     file=sprintf("%s.T", matrix_path),
        #     row.names = FALSE,
        #     col.names = FALSE,
        #     quote = FALSE,
        #     sep = "\t"
        # )

        # Do not add 1 to all counts
        read.csv(
            matrix_path,
            header = FALSE,
            stringsAsFactors = FALSE,
            row.names = NULL,
        ) %>% t %>%  write.table(
            file=sprintf("%s.T", matrix_path),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
        mcell_import_scmat_tsv(
            mat_nm = mat_name,
            fn = sprintf("%s.T", matrix_path),
            dset_nm = "dset"
        )
        unlink(
            sprintf("%s.T", matrix_path)
        )
        return(TRUE)
    }
    stop("Matrix path is not a directory or a file")
}

custom_metacell_pipeline <- function(
    mat_id="mat",
    gstat_id="gstat",
    gset_id="gset",
    graph_id="graph",
    coc_id="coc",
    mc_id="mc",
    mc2d_id="mc2d",
    scm_n_downsamp_gstat=NULL,
    T_tot=NULL, 
    T_top3=NULL,
    T_szcor=NULL,
    T_vm=NULL,
    T_niche=NULL,
    K=NULL,
    min_mc_size=NULL,
    alpha=NULL
)   {
    # # select gene markers for metacell grouping
    # if (!is.null(scm_n_downsamp_gstat)) {
    #     message(
    #         "custom_metacell_pipeline: ",
    #         "Setting n_downsamp_gstat to ", scm_n_downsamp_gstat
    #     )
    #     tgconfig::set_param(
    #         "scm_n_downsamp_gstat",
    #         as.integer(scm_n_downsamp_gstat),
    #         "metacell"
    #     )
    #     if (scm_n_downsamp_gstat != tgconfig::get_param("scm_n_downsamp_gstat", "metacell")) {
    #         stop(
    #             "custom_metacell_pipeline: ",
    #             "Failed to set n_downsamp_gstat to ", scm_n_downsamp_gstat
    #         )
    #     } 
    # }
    
    # get marker genes if not already computed
    if (is.null(scdb_gset(gset_id))) {
        message("custom_metacell_pipeline: Computing gene markers")
        # calculate gene stats
        mcell_add_gene_stat(
            gstat_id = gstat_id,
            mat_id = mat_id,
            force = T
        ) # compute gene stats
        print(summary(scdb_gstat(gstat_id)))
        mcell_gset_filter_multi(
            gstat_id = gstat_id,
            gset_id = gset_id,
            T_tot = T_tot,
            T_top3 = T_top3,
            T_szcor = T_szcor,
            T_vm = T_vm,
            T_niche = T_niche,
            force_new = T
        ) # select gene markers
    }
    # group cells
    mcell_add_cgraph_from_mat_bknn(
        mat_id = mat_id,
        gset_id = gset_id,
        graph_id = graph_id,
        K = K,
        dsamp = F
    ) # build graph (K is the target number of edges)
    mcell_coclust_from_graph_resamp(
        coc_id = coc_id,
        graph_id = graph_id,
        min_mc_size = min_mc_size,
        p_resamp = 0.75,
        n_resamp = 1000
    ) # cocluster by bootstraping
    source("scripts/run-metacell/MY_mcell_mc_from_coclust_balanced.R")
    MY_mcell_mc_from_coclust_balanced(
        mat_id = mat_id,
        coc_id = coc_id,
        mc_id = mc_id,
        K = K,
        min_mc_size = min_mc_size,
        alpha = alpha
    ) # final grpah
    # force-directed projection of the k-nn graph
    mcell_mc2d_force_knn(
        mc2d_id = mc2d_id,
        mc_id = mc_id,
        graph_id = graph_id
    )
}

plot_metacell <- function(
    mc2d_id = "mc2d",
    cell_ident,
    palette = NULL,
    mc_projections = NULL
) {
    # get sc_projections
    sc_projections <- data.frame(
        "X" = scdb_mc2d(mc2d_id)@sc_x,
        "Y" = scdb_mc2d(mc2d_id)@sc_y
    )
    sc_projections$Barcode <- rownames(sc_projections)

    # merge cell_ident to sc_projections
    sc_projections %<>% dplyr::inner_join(cell_ident, by = c("Barcode" = "Barcode"))
    sc_projections$Cluster <- as.factor(sc_projections$Cluster)

    # create palette if not provided
    if (is.null(palette)) {
        c2c_mapping <- sc_projections[, c("Cluster", "Color")] %>% unique
        palette <- c2c_mapping$Color
        names(palette) <- c2c_mapping$Cluster %>% as.character
    }
    
    # plot
    p <- ggplot(
        sc_projections,
        aes(X,Y)
    ) + geom_point(
            data = sc_projections,
            aes(X, Y, color = Cluster),
            size = .8,
            show.legend = FALSE
        ) + scale_color_manual(
            values = palette
        ) + coord_fixed(
            ratio = 1
        )
    
    if (!is.null(mc_projections)) {
        # get metacell coordinate and edges
        fr <- mc2d@graph$mc1
        to <- mc2d@graph$mc2
        # NOTE: plot only metacell number: 1,2,3,4,5,6,7,8,9,10,18,19,20
        edges <- data.frame(
            X_start = mc2d@mc_x[fr],
            Y_start = mc2d@mc_y[fr],
            X_end = mc2d@mc_x[to],
            Y_end = mc2d@mc_y[to]
        )
        edges[fr %in% c(1,2,3,4,5,6,7,8,9,10,18,19,20) & to %in% c(1,2,3,4,5,6,7,8,9,10,18,19,20),] %<>% 
            dplyr::mutate(
                X_end = NA,
                Y_end = NA
            )
        # add metacell and edges to plot
        p <- p + geom_label(
            data = mc_projections,
            aes(X, Y, label = rownames(mc_projections)),
            show.legend = FALSE,
            fill = NA,
            size = 5
        ) + geom_segment(
            data = edges,
            aes(
                x = X_start,
                y = Y_start,
                xend = X_end,
                yend = Y_end
            ),
            linewidth = 0.1
        )
    }
    return(p)
}
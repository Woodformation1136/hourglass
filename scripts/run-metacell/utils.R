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
    scm_n_downsamp_gstat=NULL
)   {
    # select gene markers for metacell grouping
    if (!is.null(scm_n_downsamp_gstat)) {
        message(
            "custom_metacell_pipeline: ",
            "Setting n_downsamp_gstat to ", scm_n_downsamp_gstat
        )
        tgconfig::set_param(
            "scm_n_downsamp_gstat",
            as.integer(scm_n_downsamp_gstat),
            "metacell"
        )
        if (scm_n_downsamp_gstat != tgconfig::get_param("scm_n_downsamp_gstat", "metacell")) {
            stop(
                "custom_metacell_pipeline: ",
                "Failed to set n_downsamp_gstat to ", scm_n_downsamp_gstat
            )
        } 
    }
    # calculate gene stats
    gstat <- scm_gene_stat(
        mat_id,
        niche_quantile = 0.2,
        downsample_n = NULL,
        K_std_n = 1
    )
    mcell_add_gene_stat(
        gstat_id = gstat_id,
        mat_id = mat_id,
        force = T
    ) # compute gene stats
    mcell_gset_filter_multi(
        gstat_id = gstat_id,
        gset_id = gset_id,
        T_tot = 20,
        T_top3 = 1,
        T_szcor = -0.01,
        T_vm = 0.2,
        T_niche = 0.05,
        force_new = T
    ) # select gene markers
    # group cells
    mcell_add_cgraph_from_mat_bknn(
        mat_id = mat_id,
        gset_id = gset_id,
        graph_id = graph_id,
        K = 40,
        dsamp = F
    ) # build graph
    mcell_coclust_from_graph_resamp(
        coc_id = coc_id,
        graph_id = graph_id,
        min_mc_size = 50,
        p_resamp = 0.75,
        n_resamp = 1000
    ) # cocluster by bootstraping
    source("scripts/run-metacell/MY_mcell_mc_from_coclust_balanced.R")
    MY_mcell_mc_from_coclust_balanced(
        mat_id = mat_id,
        coc_id = coc_id,
        mc_id = mc_id,
        K = 40,
        min_mc_size = 50,
        alpha = 3
    ) # final grpah
    # force-directed projection of the k-nn graph
    mcell_mc2d_force_knn(
        mc2d_id = mc2d_id,
        mc_id = mc_id,
        graph_id = graph_id
    )
}

plot_metacell <- function(
    sc_projections,
    mc_projections=NULL,
    palette=NULL,
    output_pdf
) {
    # create output directory 
    if (!dir.exists(dirname(output_pdf))) {
        dir.create(dirname(output_pdf), recursive = TRUE)
    }
    
    # load p
    if (is.null(palette)) {
        c2c_mapping <- sc_projections[, c("Cluster", "Color")] %>% unique
        palette <- c2c_mapping$Color
        names(palette) <- c2c_mapping$Cluster %>% as.character
    } else {
        palette <- fromJSON(palette)
    }
    
    if (is.null(mc_projections)) {
    
    }
}
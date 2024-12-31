library(Matrix)
library(magrittr)
library(R.utils)

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

build_seurat_sample_tree <- function(integration_order) {
    sample.tree <- integration_order %>%
        strsplit(",") %>%
        unlist() %>%
        as.numeric() %>%
        matrix(ncol = 2, byrow = TRUE)
    return(sample.tree)
}

csv2seurat <- function(
    csv_file,
    project_name,
    sep = ",",
    header = TRUE,
    row.names = row.names,
    min.cells = 3,
    min.features = 200,
    ...,
    verbose = TRUE
) {
    # csv_file: path to the csv file, first row is gene names, first column is cell names
    if (verbose) {
        message("Reading ", csv_file)
    }
    dat <- read.csv(csv_file, sep = sep, header = header, row.names = row.names, ...)
    if (verbose) {
        message("Converting to Seurat object")
    }
    mat <- Matrix::Matrix(
        data = as.matrix(dat[, -1]),
        sparse = TRUE
    )
    rownames(mat) <- paste(
        project_name, dat[, 1], sep = ":"
    )
    seurat_obj <- Seurat::CreateSeuratObject(
        counts = t(mat),
        project = project_name,
        min.cells = min.cells,
        min.features = min.features
    )
    return(seurat_obj)
}


qc_seurat_object <- function(
    seurat_object,
    n_feature_cutoff_on_sample,
    n_barcode_cutoff_on_sample
) {
    message("QC on Seurat object: ", seurat_object@project.name)
    n_feature <- Features(seurat_object) %>% length()
    n_barcode <- Cells(seurat_object) %>% length()
    is_pass_feature <- (n_feature >= n_feature_cutoff_on_sample)
    is_pass_barcode <- (n_barcode >= n_barcode_cutoff_on_sample)
    is_pass_both <- is_pass_feature & is_pass_barcode
    message(
        "n_feature: ", n_feature, ": ",
        ifelse(is_pass_feature, "pass", "failed")
    )
    message(
        "n_barcode: ", n_barcode, ": ",
        ifelse(is_pass_barcode, "pass", "failed")
    )
    return(is_pass_both)
}

custom_integration_pipeline_1 <- function(
    object.list,
    sample.tree,
    n_pcs = 30,
    nfeatures = 2000,
    verbose = TRUE
) {
    message("Running custom integration pipeline")
    message("Reference:", "https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette")   
    # sctranform
    message("Running SCTransform on each object using `sapply`")
    object.list <- sapply(
        object.list,
        SCTransform,
        vst.flavor = "v2",
        verbose = verbose
    )
    # integrate
    features <- SelectIntegrationFeatures(
        object.list,
        nfeatures = nfeatures,
        verbose = verbose
    )
    object.list <- PrepSCTIntegration(
        object.list = object.list,
        anchor.features = features,
        verbose = verbose
    )
    anchors <- FindIntegrationAnchors(
        object.list = object.list,
        normalization.method = "SCT",
        anchor.features = features,
        verbose = verbose
    )
    combined.sct <- IntegrateData(
        anchorset = anchors,
        sample.tree = sample.tree,
        normalization.method = "SCT",
        verbose = verbose
    )
    combined.sct <- RunPCA(
        combined.sct,
        npcs = n_pcs,
        verbose = verbose
    )
    return(combined.sct)
}

# SCTransform & pca pipeline for single sample
SCTransform_and_PCA <- function(
    seurat_object,
    n_pcs = 30,
    verbose = TRUE
) {
    message("Perform SCTransform and PCA on the Seurat object")
    message("Reference:", "https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette")   
    message("Running SCTransform on the object")
    seurat_object <- SCTransform(
        seurat_object,
        vst.flavor = "v2",
        verbose = verbose
    )
    message("Running PCA on the object")
    seurat_object <- RunPCA(
        seurat_object,
        npcs = n_pcs,
        verbose = verbose
    )
    return(seurat_object)
}
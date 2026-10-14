library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)
library(tools)
library(R.utils)


# parse arguments
args <- R.utils::commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]
keys <- attachLocally(args)
# cat("Command-line arguments attached to global environment:\n");
# print(str(mget(keys, envir=globalenv())))

# # for debugging
# samples2integrate <- "outputs/seurat_integration/ortho_umi_count/tung_batch1_rs17.csv,outputs/seurat_integration/ortho_umi_count/egr_nextseq_batch1_rs17.csv"
# integration_order <- "-1,-2"
# k_nearest_neighbor <- "3"
# n_feature_cutoff_on_sample <- "2000"
# n_barcode_cutoff_on_sample <- "100"
# n_pcs <- "30"                      
# output_pca_projection_csv <- "outputs/seurat_integration/group1/rs17/pca_projection.csv"
# output_clusters_csv <- "outputs/seurat_integration/group1/rs17/clusters.csv"
# output_corrected_expression_matrix_csv <- "outputs/seurat_integration/group1/rs17/corrected_expression_matrix.csv"

# process arguments
samples2integrate <- samples2integrate %>% strsplit(",") %>% unlist()
integration_order <- integration_order %>% 
    strsplit(",") %>% 
    unlist() %>%
    as.numeric() %>% 
    matrix(ncol = 2, byrow = TRUE)
k_nearest_neighbor %<>% as.integer()
n_feature_cutoff_on_sample %<>% as.integer()
n_barcode_cutoff_on_sample %<>% as.integer()
n_pcs %<>% as.integer()
output_pca_projection_csv %<>% as.character()
output_clusters_csv %<>% as.character()
output_corrected_expression_matrix_csv %<>% as.character()
write(
    paste0(
        "Params:\n",
        str(mget(keys, envir=globalenv()))
    )
    , stderr()
)


# Define functions =============================================================
# Input and preprocess data
preprocess_data <- function(input_UMI_csv, project_name) {
    input_UMI_df <- read.csv(input_UMI_csv)
    input_UMI_matrix <- Matrix::Matrix(as.matrix(input_UMI_df[, -1]))
    rownames(input_UMI_matrix) <-
        paste0(project_name, "_", input_UMI_df$Barcode)

    out <-
        # Create Seurat Object
        Seurat::CreateSeuratObject(
            counts = t(input_UMI_matrix),
            project = project_name,
            min.cells = 3,
            min.features = 200
        ) %>%
        # Normalize the data
        Seurat::NormalizeData(
            normalization.method = "LogNormalize",
            scale.factor = 10000
        ) %>%
        # Identify the highly variable features (feature selection)
        Seurat::FindVariableFeatures(
            selection.method = "vst",
            nfeatures = 2000
        )

    return(out)
}

# Integrate different SCseq samples
integrate_and_clustering <- function(
    object_list,
    integration_order,
    k_nearest_neighbor,
    n_pcs
) {
    #Find integration anchors and integrate data
    integration_anchors <- FindIntegrationAnchors(
        object.list = object_list,
        anchor.features = 2000,
        scale = TRUE,
        reduction = "cca",
        l2.norm = TRUE,
        k.anchor = 5
    )

    #Run the standard workflow for visualization and clustering
    combined_object <- integration_anchors %>%
        IntegrateData(
            sample.tree = integration_order,
            preserve.order = TRUE,
            k.weight = 100
        ) %>%
        ScaleData() %>%
        RunPCA(
            npcs = n_pcs
        ) %>%
        FindNeighbors(
            reduction = "pca",
            dims = 1:n_pcs,
            k.param = k_nearest_neighbor
        ) %>%
        FindClusters(
            resolution = 0.5
        )

    return(combined_object)
}


# ========================= main =========================
# Create output directory
if (!dir.exists(dirname(output_pca_projection_csv))) {
    dir.create(dirname(output_pca_projection_csv), recursive = TRUE)
}

# Input and pre-process UMI matrix into a list
UMI_list <- sapply(
    samples2integrate,
    function(x) {
        out <- preprocess_data(
            input_UMI_csv = x,
            project_name = file_path_sans_ext(basename(x))
        )
        return(out)
    }
)


# QC samples. Abort if any sample does not pass QC.
# message("Sample QC:")
qc_sample <- function(seurat_object) {
    message("> ", seurat_object@project.name)
    n_feature <- seurat_object %>%  Features() %>% length()
    n_barcode <- seurat_object %>% Cells() %>% length()
    is_pass_feature <- (n_feature >= n_feature_cutoff_on_sample)
    is_pass_barcode <- (n_barcode >= n_barcode_cutoff_on_sample)
    message(
        ">> n_feature: ", n_feature, ": ",
        ifelse(is_pass_feature, "pass", "failed")
    )
    message(
        ">> n_barcode: ", n_barcode, ": ",
        ifelse(is_pass_barcode, "pass", "failed")
    )
    out <- is_pass_feature & is_pass_barcode
    return(out)
}
passed <- sapply(UMI_list, qc_sample)
if (any(!passed)) stop("Sample(s) failed QC.")


# 
combined_object <- integrate_and_clustering(
    object_list = UMI_list,
    integration_order = integration_order,
    k_nearest_neighbor = k_nearest_neighbor,
    n_pcs = n_pcs
)

# write PCs to csv
write.csv(
    x = Embeddings(
        object = combined_object,
        reduction = "pca"
    ),
    file = output_pca_projection_csv,
    row.names = TRUE,
    quote = FALSE
)

# write clusters to csv
cluster_df <- as.data.frame(combined_object@meta.data$seurat_clusters)
cluster_df <- cbind(
    Barcode = Cells(combined_object),
    cluster_df
)
colnames(cluster_df) <- c("Barcode", "Cluster")
write.csv(
    x = cluster_df,
    file = output_clusters_csv,
    row.names = TRUE,
    quote = FALSE
)

# write batch-corrected UMI matrix to csv
write.csv(
    x = GetAssayData(
        combined_object,
        "integrated",
        "data"
    ) %>% t,
    file = output_corrected_expression_matrix_csv,
    row.names = TRUE,
    quote = FALSE
)
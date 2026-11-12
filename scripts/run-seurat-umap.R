# Run with docker image: satijalab/seurat:5.0.0
library(Seurat)
library(magrittr)

# define functions
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
    # csv_file: path to the csv file, first row is cell names, first column is gene names
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
        counts = mat,
        project = project_name,
        min.cells = min.cells,
        min.features = min.features
    )
    return(seurat_obj)
}

# load Seurat data from CSV (first row is cell names, first column is gene names)
seurat_obj <- csv2seurat(
    csv_file = "umi_table.csv",
    project_name = "project_name"
)

# run NormalizeData > ScaleData > PCA > UMAP
random_seed <- 17
seurat_obj %<>% NormalizeData
seurat_obj %<>% ScaleData
seurat_obj %<>% FindVariableFeatures 
seurat_obj %<>% RunPCA
seurat_obj %<>% RunUMAP(
    reduction = "pca",
    dims = 1:10,
    n.neighbors = 30,
    n.components = 2L,
    metric = "correlation",
    min.dist = 0.5,
    seed.use = random_seed
)

# plot UMAP
svg("umap.svg")
DimPlot(seurat_obj)
dev.off()

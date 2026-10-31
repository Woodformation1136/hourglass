suppressPackageStartupMessages({
  library(sccore)
  library(slingshot)
  library(dtw)
})

# Pairwise similarity
# ------------------------------------------------------------------------------
# log2CPM on a gene expression matrix (normalized each column)
logcpm <- function(x) {
  log2(1 + 1e6 * sweep(as.matrix(x), 2, colSums(x), "/"))
}

# Calculate Jensen-Shannon similarity; return a matrix similar to `cor`
calculate_js_similarity <- function(mat1, mat2) {
  out <- sccore::jsDist(cbind(mat1, mat2))[
    1:ncol(mat1),
    (ncol(mat1)+1):(ncol(mat1)+ncol(mat2))
  ]
  rownames(out) <- colnames(mat1)
  colnames(out) <- colnames(mat2)
  return(1 - out)
}

# Calculate pairwise similarity between pseudobulk profiles
# Available metrics are: `spearman`, `pearson`, and `Jensen-Shannon`
calculate_similarity <- function(mat1, mat2, metrics) {
  stopifnot(metrics %in% c("spearman", "pearson", "Jensen-Shannon"))
  if (metrics %in% c("spearman", "pearson")) {
    cor(
      x = logcpm(mat1),
      y = logcpm(mat2),
      method = metrics
    )
  } else if (metrics == "Jensen-Shannon") {
    calculate_js_similarity(
      as.matrix(mat1),
      as.matrix(mat2)
    )
  }
}
# ------------------------------------------------------------------------------

# Run slingshot
# ------------------------------------------------------------------------------
slingshot_manual <- function(seurat_obj) {
  lineages <- new("SlingshotDataSet")
  
  # add cell embedding
  lineages@reducedDim <- seurat_obj %>% 
    Embeddings(reduction = "pca") %>% 
    extract(1:dim(seurat_obj)[2], 1:10) %>%
    as.matrix
  
  # add cluster matrix
  clusters <- seurat_obj$cell_label # factor
  mtx_clus <- matrix(0, nrow=length(clusters), ncol=length(levels(clusters)))
  row.names(mtx_clus) <- names(clusters)
  colnames(mtx_clus) <- levels(clusters)
  for (i in seq(length(clusters))) {
    mtx_clus[names(clusters)[i], clusters[i]] <- 1
  }
  lineages@clusterLabels <- mtx_clus
  
  # add lineages
  lineages@lineages <- lineage_cell_types
  
  # add adjacency matrix
  mtx_adj <- matrix(0, nrow = length(levels(clusters)), ncol = length(levels(clusters)))
  row.names(mtx_adj) <- colnames(mtx_adj) <- levels(clusters)
  for (l_i in lineages@lineages) {
    for (c_i in seq(length(l_i) - 1)) {
      mtx_adj[l_i[c_i], l_i[c_i + 1]] <- 1
      mtx_adj[l_i[c_i + 1], l_i[c_i]] <- 1
    }
  }
  lineages@adjacency <- mtx_adj
  
  # fit curves
  curves <- getCurves(
    lineages,
    approx_points = 100,
    reassign = FALSE
  ) %>% SlingshotDataSet
  
  # estimate pseudotime
  pseudotime <- slingPseudotime(
    curves,
    na = FALSE
  )
  branch_id <- slingBranchID(
    curves
  )
  return(
    list(
      SlingshotDataSet = curves,
      pseudotime = pseudotime,
      branch_id = branch_id
    )
  )
}
# ------------------------------------------------------------------------------

# Bin expression by pseudotime
# ------------------------------------------------------------------------------
bin_expression <- function(
    seurat_obj,
    pseudotime_mtx, 
    branch_id, 
    lineage,
    exclude = c(),
    prefix, 
    n_bins
) {
  # get pseudotime values
  pt <- pseudotime_mtx[, lineage]
  
  # keep only cells belongs to the lineage
  keep <- branch_id %>% grepl(lineage, .)
  keep <- keep & (!seurat_obj$cell_label %in% exclude)
  pt <- pt[keep]
  
  # set bins 
  bins <- cut(
    pt, 
    breaks = seq(
      min(pt),
      max(pt),
      length.out = n_bins + 1    
    ),
    labels = FALSE, 
    include.lowest = TRUE
  )
  names(bins) <- names(pt)
  
  # subset lineage cells
  obj <- subset(seurat_obj, cells = names(pt))
  obj$bin <- sprintf("%s-%02d", prefix, bins[colnames(obj)])
  
  # aggregate pseudotime
  pb <- AggregateExpression(obj, assays = "RNA", group.by = "bin")$RNA
  pb <- pb[, order(colnames(pb)), drop = FALSE]
  
  # report cell label composition by bins
  comp <- table(
    bin = bins,
    cell_type = seurat_obj$cell_label[keep]
  )
  return(list(Pseudobulk = pb, Composition = comp))
}
# ------------------------------------------------------------------------------

# Find middle module by marker gene expression
# ------------------------------------------------------------------------------
get_markers <- function(orthogroup, ptr_genes, sample_id) {
  if (!"ptr" %in% sample_id) {
    # get orthogroup id
    cid <- orthogroup %>% 
      dplyr::filter(gene_name %in% ptr_genes) %>%
      dplyr::pull(cluster_id)
    
    # get species gene id
    species_marker <- orthogroup %>%
      dplyr::filter(
        cluster_id %in% cid,
        species_name == switch(
          strsplit(sample_id, "_")[[1]][1],
          "ptr" = "PoT",
          "egr" = "EuG",
          "lch" = "LiC",
          "tar" = "TrA"
        )
      ) %>%
      dplyr::pull(gene_name)
  } else {
    species_marker <- ptr_genes
  }
}
get_marker_profile <- function(pb_mat, marker, norm_func = logcpm) {
  colMeans(
    norm_func(pb_mat)[marker, , drop = FALSE]
  )
}
get_middle_module <- function(pb_mat, marker) {
  marker <- intersect(rownames(pb_mat), marker)
  marker_profile <- get_marker_profile(pb_mat, marker, logcpm) 
  middle_peak <- marker_profile %>% which.max
  return(
    list(
      early = 1:(middle_peak-2),
      middle = (middle_peak-1):(middle_peak+1),
      late = (middle_peak+2):ncol(pb_mat)
    )
  )
}
# ------------------------------------------------------------------------------

# Align binned-stages by DTW
# ------------------------------------------------------------------------------
align_similarity_matrix <- function(
    similarity,
    window_size = ceiling(as.numeric(params$window_dtw) * max(dim(similarity)))
) {
  local_cost <- 1 - similarity
  fit <- dtw(
    x = local_cost,
    y = NULL,
    step.pattern = dtw::asymmetric, # Each ray bin is matched exactly once.
    window.type = "slantedband",
    window.size = window_size,
    open.begin = TRUE,
    open.end = TRUE,
    keep.internals = TRUE,
    distance.only = FALSE
  )
  
  path <- tibble::tibble(
    L1_i = fit$index1,
    L2_i = fit$index2,
    L1_bin = rownames(similarity)[fit$index1],
    L2_bin = colnames(similarity)[fit$index2],
    similarity = similarity[
      cbind(fit$index1, fit$index2)
    ]
  )
  list(
    fit = fit,
    path = path,
    local_cost = local_cost
  )
}
# ------------------------------------------------------------------------------

# Summarize cell type alignment correspondence
# ------------------------------------------------------------------------------
composition_proportions <- function(comp, states) {
  comp <- as.matrix(comp)
  out <- matrix(
    0, nrow = nrow(comp), ncol = length(states), 
    dimnames = list(rownames(comp), states)
  )
  shared_states <- intersect(colnames(comp), states)
  
  out[, shared_states] <- comp[, shared_states, drop = FALSE]
  
  totals <- rowSums(out)
  
  sweep(out, 1, totals, "/")
}
extract_bin_id <- function(x) {
  x %>% 
    sub("^.*-(\\d+)$", "\\1", .) %>%
    as.integer %>%
    as.character
}
alignment_label_mass <- function(
    aligned,
    comp_query,
    comp_reference,
    query_states,
    reference_states
) {
  # get cell state proportion in
  query_prop <- composition_proportions(comp_query, query_states)
  reference_prop <- composition_proportions(comp_reference, reference_states)
  
  # get alignment path
  path <- aligned$path
  query_bin <- extract_bin_id(path$L1_bin)
  reference_bin <- extract_bin_id(path$L2_bin)
  
  
  mass <- matrix(
    0, nrow = length(query_states), ncol = length(reference_states),
    dimnames = list(query_state = query_states, reference_state = reference_states)
  )
  
  for (k in seq_len(nrow(path))) {
    mass <- mass + tcrossprod(
      query_prop[query_bin[k], , drop = TRUE],
      reference_prop[reference_bin[k], , drop = TRUE]
    )
  }
  
  # Equal weight for each Ray pseudotime position.
  mass <- mass / nrow(path)
  
  stopifnot(
    isTRUE(all.equal(sum(mass), 1, tolerance = 1e-8))
  )
  
  mass
}
# ------------------------------------------------------------------------------

# Calculate hourglass score
# ------------------------------------------------------------------------------
hourglass_score <- function(S, modules) {
  E <- mean(S[modules$early], na.rm = TRUE)
  M <- mean(S[modules$middle], na.rm = TRUE)
  L <- mean(S[modules$late], na.rm = TRUE)
  c(
    Early = E,
    Middle = M,
    Late = L,
    H = min(M - E, M - L)
  )
}
# ------------------------------------------------------------------------------

# Bootstrapping
# ------------------------------------------------------------------------------
# Re-sampling cells from a SeuratObject with replacement
make_bootstrap_seurat <- function(
    seurat_obj,
    seed,
    assay = "RNA"
) {
  set.seed(seed)
  
  n_cells <- ncol(seurat_obj)
  
  # Bootstrap draw: integer indices preserve repeated selections.
  sampled_idx <- sample.int(
    n = n_cells,
    size = n_cells,
    replace = TRUE
  )
  source_cells <- Cells(seurat_obj)[sampled_idx]
  
  # Extract raw counts, preserving duplicated columns.
  counts_boot <- LayerData(
    seurat_obj,
    assay = assay,
    layer = "counts"
  )[, sampled_idx, drop = FALSE]
  
  # Every bootstrap draw requires its own unique object-level ID.
  bootstrap_cells <- sprintf(
    "Boot-%05d__%s",
    seq_len(n_cells),
    source_cells
  )
  colnames(counts_boot) <- bootstrap_cells
  
  # Copy metadata according to the same sampled indices.
  meta_boot <- seurat_obj[[]][
    sampled_idx,
    ,
    drop = FALSE
  ]
  meta_boot$source_cell <- source_cells
  meta_boot$bootstrap_draw <- seq_len(n_cells)
  rownames(meta_boot) <- bootstrap_cells
  
  # Let CreateSeuratObject recalculate these fields.
  qc_columns <- intersect(
    c(
      paste0("nCount_", assay),
      paste0("nFeature_", assay)
    ),
    colnames(meta_boot)
  )
  meta_boot <- meta_boot[
    ,
    setdiff(colnames(meta_boot), qc_columns),
    drop = FALSE
  ]
  object_boot <- CreateSeuratObject(
    counts = counts_boot,
    assay = assay,
    meta.data = meta_boot,
    project = paste0(Project(seurat_obj), "_bootstrap")
  )
  
  # Preserve biological labels.
  if ("cell_label" %in% colnames(object_boot[[]])) {
    Idents(object_boot) <- "cell_label"
  }
  
  list(
    object = object_boot,
    sampled_idx = sampled_idx,
    source_cells = source_cells
  )
}
# ------------------------------------------------------------------------------


# Pairwise marker discovery and curation: an exploratory Seurat workflow.
# Input: a Seurat object with log-normalized, uncorrected RNA data and cell labels.
# This does not establish independent annotation validity or biological replication.
# References:
# https://bioconductor.org/books/release/OSCA.basic/marker-detection.html
# https://www.bioconductor.org/books/release/SingleRBook/classic-mode.html
# https://satijalab.org/seurat/reference/findmarkers

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# Test each unordered pair once; reuse its two-sided Wilcoxon result in reverse.
# No marker-list or single-copy prefilter is used.
pairwise_cluster_markers <- function(obj, label_col = "cell_label", min_pct = 0.10) {
  if (!label_col %in% colnames(obj[[]])) stop("Missing metadata: ", label_col)
  if (anyDuplicated(colnames(obj))) stop("Cell names are not unique.")
  if (!"RNA" %in% SeuratObject::Assays(obj)) stop("An RNA assay is required.")
  if (length(min_pct) != 1L || !is.finite(min_pct) || min_pct < 0 || min_pct > 1) {
    stop("min_pct must be a number between zero and one.")
  }
  DefaultAssay(obj) <- "RNA"
  lab <- as.character(obj[[label_col, drop = TRUE]])
  keep <- !is.na(lab) & nzchar(trimws(lab)) & !lab %in% c("NA", "Unassigned", "Unknown")
  if (!any(keep)) stop("No annotated cells remain.")
  # Reset locally before subsetting to avoid the missing-active-identity issue.
  Idents(obj) <- rep("plot", ncol(obj))
  obj <- subset(obj, cells = colnames(obj)[keep])
  Idents(obj) <- label_col
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj, assay = "RNA")
  x <- GetAssayData(obj, assay = "RNA", layer = "data")
  if (nrow(x) == 0L || ncol(x) == 0L) stop("Normalize the RNA assay before testing.")
  rm(x)
  clusters <- sort(unique(as.character(Idents(obj))))
  if (length(clusters) < 2L) stop("At least two annotated clusters are needed.")
  sizes <- table(factor(as.character(Idents(obj)), levels = clusters))
  small <- names(sizes[sizes < 3L])
  if (length(small)) stop("Fewer than three cells in: ", paste(small, collapse = ", "))

  pairs <- combn(clusters, 2L, simplify = FALSE)
  result <- bind_rows(lapply(pairs, function(z) {
    message(z[1], " versus ", z[2])
    raw <- FindMarkers(
      object = obj, assay = "RNA", slot = "data",
      ident.1 = z[1], ident.2 = z[2], test.use = "wilcox",
      only.pos = FALSE, logfc.threshold = 0, min.pct = min_pct,
      base = 2, fc.name = "avg_log2FC", verbose = FALSE
    )
    if (!nrow(raw)) stop("No results for ", z[1], " versus ", z[2])
    forward <- raw %>% rownames_to_column("gene") %>% transmute(
      cluster = z[1], comparator = z[2], gene,
      log2FC = avg_log2FC, p_value = p_val, padj = p_val_adj,
      pct_target = pct.1, pct_other = pct.2
    )
    reverse <- forward %>% transmute(
      cluster = z[2], comparator = z[1], gene,
      log2FC = -log2FC, p_value, padj,
      pct_target = forward$pct_other, pct_other = forward$pct_target
    )
    bind_rows(forward, reverse)
  }))
  attr(result, "clusters") <- clusters
  attr(result, "settings") <- list(
    assay = "RNA", layer = "data", test = "wilcox", min_pct = min_pct,
    label_col = label_col, Seurat = as.character(packageVersion("Seurat"))
  )
  result
}

# Summarize one target's contrasts without rerunning any tests.
# Default cutoffs deliberately retain the original FC >= 1 requirement.
# pct_other is reported, not subjected to the earlier near-exclusivity cap.
inspect_marker_overlap <- function(pairwise, target, min_log2fc = 1,
                                   min_target_pct = 0.25, alpha = 0.05,
                                   top_n = 3L) {
  required <- c("cluster", "comparator", "gene", "log2FC", "padj", "pct_target", "pct_other")
  if (!all(required %in% colnames(pairwise))) stop("Unexpected pairwise table columns.")
  if (!target %in% pairwise$cluster) stop("Target not found: ", target)
  if (length(top_n) != 1L || !is.finite(top_n) || top_n < 1 || top_n != as.integer(top_n)) {
    stop("top_n must be a positive integer.")
  }
  if (!is.finite(min_log2fc) || min_log2fc <= 0 ||
      !is.finite(min_target_pct) || min_target_pct < 0 || min_target_pct > 1 ||
      !is.finite(alpha) || alpha <= 0 || alpha > 1) stop("Invalid screening thresholds.")
  clusters <- attr(pairwise, "clusters")
  if (is.null(clusters)) clusters <- sort(unique(c(pairwise$cluster, pairwise$comparator)))
  others <- setdiff(clusters, target)
  k <- length(others)
  tab <- pairwise %>% filter(.data$cluster == .env$target)
  if (anyDuplicated(tab[c("gene", "comparator")])) stop("Duplicate gene-comparator rows.")
  tab <- tab %>% select(all_of(setdiff(required, "cluster"))) %>%
    complete(gene, comparator = others) %>% mutate(
      cluster = target,
      tested = is.finite(log2FC) & is.finite(padj) &
        is.finite(pct_target) & is.finite(pct_other),
      pass = if_else(tested,
        log2FC >= min_log2fc & padj < alpha & pct_target >= min_target_pct,
        NA),
      detection_gap = pct_target - pct_other
    )
  finite_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  }
  summary <- tab %>% group_by(gene) %>% summarise(
    n_tested = sum(tested), n_supported = sum(pass, na.rm = TRUE),
    min_log2FC = if (all(tested)) min(log2FC) else NA_real_,
    target_detection = finite_max(pct_target),
    max_offtarget_detection = finite_max(pct_other),
    weakest_comparator = if (all(tested)) comparator[which.min(log2FC)] else NA_character_,
    supported_vs = paste(comparator[pass %in% TRUE], collapse = "; "),
    unsupported_vs = paste(comparator[pass %in% FALSE], collapse = "; "),
    untested_vs = paste(comparator[!tested], collapse = "; "), .groups = "drop"
  ) %>% mutate(cluster = target, n_comparators = k,
    all_comparators = n_tested == k & n_supported == k) %>%
    arrange(desc(n_supported), desc(min_log2FC), max_offtarget_detection)

  # All passing genes define overlap; top_n is only a review-list size cap.
  sets <- setNames(lapply(others, function(z) {
    tab$gene[tab$comparator == z & tab$pass %in% TRUE]
  }), others)
  overlap <- outer(others, others, Vectorize(function(a, b) {
    length(intersect(sets[[a]], sets[[b]]))
  }))
  dimnames(overlap) <- list(others, others)
  selected <- tab %>% filter(pass %in% TRUE) %>% group_by(comparator) %>%
    arrange(desc(log2FC), desc(detection_gap), gene, .by_group = TRUE) %>%
    slice_head(n = top_n) %>% ungroup() %>% distinct(gene) %>% pull(gene)
  candidates <- summary %>% filter(gene %in% selected) %>% mutate(
    marker_role = if_else(all_comparators, "All-comparator", "Complementary"),
    keep_after_review = NA, curation_note = ""
  )
  membership <- tab %>% filter(gene %in% selected) %>%
    select(gene, comparator, pass) %>% mutate(pass = as.integer(pass)) %>%
    pivot_wider(names_from = comparator, values_from = pass)
  membership <- membership[match(candidates$gene, membership$gene), , drop = FALSE]
  coverage <- tab %>% group_by(comparator) %>% summarise(
    n_passing = sum(pass, na.rm = TRUE),
    n_in_review_panel = sum(pass %in% TRUE & gene %in% selected), .groups = "drop"
  )
  list(
    summary = summary,
    intersection = summary %>% filter(all_comparators),
    sets = sets, overlap = overlap, membership = membership,
    coverage = coverage, candidates = candidates, pairwise = tab
  )
}

# Example (run explicitly after sourcing; these lines are intentionally commented):
# pw <- pairwise_cluster_markers(merged_obj)
# saveRDS(pw, "pairwise-cluster-markers.rds")
# rp <- inspect_marker_overlap(pw, "Ray Precursor", top_n = 3)
# rp$overlap       # diagonal = set size, off-diagonal = intersection size
# rp$membership    # 1 = passes, 0 = tested but fails, NA = missing/invalid
# rp$coverage      # zero passing genes identifies an unresolved comparison
# rp$intersection  # may legitimately contain zero genes
# write.csv(rp$candidates, "RP-marker-review.csv", row.names = FALSE)
# Keep only genes accepted after within-sample expression review.
# Annotate single-copy status afterward using the FULL orthology membership map,
# not the old single_copy_marker_sc list, which omits newly discovered candidates.

library(metacell)
library(tgstat)
library(tgconfig)

MY_mcell_coclust_filt_by_k_deg <- function(coc_id, K, alpha) {
	coc = scdb_coclust(coc_id)
	if(is.null(coc)) {
		stop("MC-ERR: coclust ", coc_id , " is missing when trying to derive K-deg filter")
	}

	edges = coc@coclust
	
	deg_wgt = as.matrix(table(c(edges$node1, edges$node2), c(edges$cnt,edges$cnt)))
	deg_cum = t(apply(deg_wgt, 1, function(x) cumsum(rev(x))))
	thresh_Kr = rowSums(deg_cum > K)
	lev2int = 1:length(levels(edges$node1))
	thresh_K = rep(NA, length(levels(edges$node1)))
	names(thresh_K) = levels(edges$node1)
	if (is.na(sum(as.numeric(names(thresh_Kr))))){
		thresh_K[names(thresh_Kr)] = thresh_Kr
	}
	else {
		thresh_K[as.numeric(names(thresh_Kr))] = thresh_Kr
	}

	filt_edges = thresh_K[edges$node1] < edges$cnt * alpha | 
							thresh_K[edges$node2] < edges$cnt * alpha

	return(filt_edges)
}

MY_mcell_mc_from_coclust_balanced <- function(
    mc_id,
    coc_id,
	mat_id,
    K,
    min_mc_size,
    alpha=2
){
    # copied un-exported function ".set_seed", ".restore_seed" from metacell
    .set_seed = function(rseed = 1) {
        if (exists(".Random.seed", .GlobalEnv)) {
            oldseed = .GlobalEnv$.Random.seed
        }
        else {
            oldseed = NULL
        }
        #message("will set seed")
        set.seed(seed=rseed)

        return(oldseed)
    }
    .restore_seed = function(oldseed) {
        if (!is.null(oldseed)) {
            .GlobalEnv$.Random.seed = oldseed
        }
        else {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }
    old_seed = .set_seed(get_param("mc_rseed", "metacell"))

	tgs_clust_cool = get_param("scm_tgs_clust_cool", "metacell")
	tgs_clust_burn = get_param("scm_tgs_clust_burn_in", "metacell")

	coc = scdb_coclust(coc_id)
	if(is.null(coc)) {
		stop("MC-ERR: coclust ", coc_id , " is missing when running mc from coclust")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR: mat id ", mat_id, " is missing when running add_mc_from_graph")
	}

	edges = coc@coclust
	filt_edges = MY_mcell_coclust_filt_by_k_deg(coc_id, K, alpha)

	message("filtered ", nrow(edges) - sum(filt_edges), " left with ", sum(filt_edges), " based on co-cluster imbalance")
	edges = edges[filt_edges,]

	colnames(edges) = c("col1", "col2", "weight")
	edges$weight = edges$weight/max(edges$weight)

	edges = edges[edges$col1 != edges$col2,]
	redges = edges[,c("col1", "col2", "weight")]
	redges$col1= edges$col2
	redges$col2= edges$col1
	edges = rbind(edges,redges)

	node_clust = tgs_graph_cover(edges, min_mc_size,
					cooling = tgs_clust_cool, burn_in = tgs_clust_burn)

	f_outlier = (node_clust$cluster == 0)

	outliers = colnames(mat@mat)[node_clust$node[f_outlier]]
	mc = as.integer(as.factor(node_clust$cluster[!f_outlier]))
	names(mc) = colnames(mat@mat)[!f_outlier]
	message("building metacell object, #mc ", max(mc))
	cell_names = colnames(mat@mat)
	scdb_add_mc(mc_id, tgMCCov(mc, outliers, mat))
	message("reordering metacells by hclust and most variable two markers")

    
	.restore_seed(old_seed)

	mcell_mc_reorder_hc(mc_id)
}
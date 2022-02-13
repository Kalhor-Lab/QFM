library(BiocParallel)
simulate_bulk <- function(type_graph, mut_p, sample_size) {
        suppressMessages(devtools::load_all("."))
        suppressMessages(require(extraDistr))
        suppressMessages(require(matrixStats))
        produce_data <- function(type_graph, mut_p, sample_size) {
                # set.seed(73)
                mut_param = do.call(make_mut_param_by_rate, mut_p)
                gens0 = make_gens(type_graph = type_graph)
                gens1 = sample_size_gens(gens0, type_graph, sample_size = sample_size)
                gens1 = simulate_muts(gens1, type_graph, mut_param = mut_param)
                allele_mats = extract_allele_matrix(gens1[type_graph$tip_id])
                allele_mats_collapsed = lapply(allele_mats, collapse_reucr_ver)
                return(allele_mats_collapsed)
        }
        return(produce_data(type_graph = type_graph,
                            mut_p = mut_p,
                            sample_size = sample_size))
}
# parallele wrapper
# bp_param = SnowParam(workers = 12, progressbar = T)
simulate_parallel <- function(sim_func, type_graph, mut_p, sample_size, num_iters, bp_param) {
        return(bplapply(1:num_iters, function(i, run_func, type_graph, mut_p, sample_size) {
                run_func(type_graph, mut_p, sample_size)
        }, BPPARAM = bp_param,
        run_func = sim_func,
        type_graph = type_graph,
        mut_p = mut_p,
        sample_size = sample_size))
}
default_dist <- function(mat, params) {
        dist(mat, method = params$dist_method)
}
custom_dist = function(mat, params) {
        1 - cor(t(scale(mat)))
}
# evaluation
produce_dist <- function(data, params, dist_func = default_dist) {
        res = lapply(data, function(x) {
                mat = bulk_pipeline(x, params)
                if (is.null(mat)) {
                        return(NULL)
                }
                return(dist_func(mat, params))
        })
        res
}
produce_phy <- function(dist_list, phy_func = nj) {
        res = lapply(dist_list, function(x) {
                if (length(x) == 0) {
                        return(NULL)
                }
                if (any(is.nan(x))) {
                        return(NULL)
                }
                if (is.null(x)) {
                        return(NULL)
                }
                phy = phy_func(x)
                # phy = phytools::midpoint.root(phy)
                phy
        })
        res
}
eval_phy <- function(phy_list, type_graph, eval_fun = phangorn::RF.dist) {
        res = do.call(rbind, lapply(phy_list, function(phy) {
                if (is.null(phy)) {
                        return(NA)
                }
                eval_fun(phy, as.phylo(type_graph))
        }))
        res
}
#' get edge length statistics to distinguish between g0 and g1
get_edge_len <- function(phy_list) {
        do.call(rbind, lapply(phy_list, function(phy) {
                # plot(phy)
                # edgelabels(format(phy$edge.length, digit = 3), 1:10, cex = 1.5)
                node_x = getMRCA(phy, tip = c("-1", "-2", "-3", "-4"))
                node_y1 = getMRCA(phy, tip = c("-1", "-2"))
                node_y2 = getMRCA(phy, tip = c("-3", "-4"))
                c(phy$edge.length[which(phy$edge[, 1] == node_x & phy$edge[, 2] == node_y1)],
                  phy$edge.length[which(phy$edge[, 1] == node_x & phy$edge[, 2] == node_y2)],
                  phy$edge.length[which(phy$edge[, 2] == node_x)])
        }))
}

sc_to_bulk <- function(muts, cell_type) {
        out = do.call(rbind, lapply(split(muts, cell_type), function(x) table(x)[unique(muts)]))
        out[is.na(out)] = 0
        out
}
bulk_to_sc <- function(mut_mat) {
        cbind(do.call(c, lapply(1:nrow(mut_mat), function(i) {
                out = rep(colnames(mut_mat), times = mut_mat[i, ])
                names(out) = paste0("type_", rownames(mut_mat)[i], "_cell_", 1:length(out))
                out
        })))
}
score_single_element <- function(y) {
        tip_id = rownames(y)
        sim_mat = matrix(0,
                         nrow = length(tip_id),
                         ncol = length(tip_id))
        rownames(sim_mat) = colnames(sim_mat) = tip_id
        if (ncol(y) == 0) {
                return(sim_mat)
        }
        for (j in 1:ncol(y)) {
                x = y[, j]
                p_score = 1/sum(x > 0)
                if(sum(x > 0) > 1 & sum(x > 0) < nrow(y)) {
                        all_comb = combn(names(x)[x>0], m = 2)
                        for (i in 1:ncol(all_comb)) {
                                sim_mat[all_comb[1, i], all_comb[2, i]] = sim_mat[all_comb[1, i], all_comb[2, i]] + p_score
                                sim_mat[all_comb[2, i], all_comb[1, i]] = sim_mat[all_comb[2, i], all_comb[1, i]] + p_score
                        }
                }
        }
        return(sim_mat)
}


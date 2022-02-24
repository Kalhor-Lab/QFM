library(qfm)
output_dir = "./intermediate_data/panel/"
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))
exp_params = readRDS(paste0(output_dir, "exp_data_10rep_wtree.rds"))

# exp_params$tr = map(exp_params$)

library(furrr)
plan(multisession, workers = 6)
total_time = 15.0 - 0.6
exp_params$gr3 = future_map(exp_params$tr3, function(tr) {
        if (is.null(tr)) {
                return(NULL)
        }
        sc_celltypes = get_type_from_id(tr$tip.label)
        names(sc_celltypes) = tr$tip.label
        gr = reconstruct_graph(tr, sc_celltypes, total_time, stat_func = mean, theta = 0.)
        gr
}, .progress = T, .options = furrr_options(seed = T))
exp_params = mutate(exp_params,
                    gr3_eval = map2(gr3, big_graph_id, function(x, i) {
                            if (is.null(x)) {
                                    return(NULL)
                            }
                            evalute_gr(x, qfm::as.phylo.type_graph(all_graphs[[i]]))
                    }))
exp_params$tr3_eval = future_map2(exp_params$tr, exp_params$tr3, function(a, b) {
        if (!is.null(a) & !is.null(b)) {
                return(evalute_gr(b, a))
        }
})
# End Phylotime part
saveRDS(exp_params, paste0(output_dir, "exp_data_10rep_proc.rds"))

#### UPGMA Manhattan ####
# parallelDist needed for Hamming example
# library(parallelDist)
# conflicted::conflict_prefer("dist", "stats")
# exp_params$tr5 = map(1:nrow(exp_params), function(j) {
#         message(j)
#         x = exp_params$data[[j]]
#         if (is.null(x$sc)) {
#                 return(NULL)
#         }
#         sc_onehot = barcode_allele_onehot_new(x$sc)
#         sc_onehot = sc_onehot[, colSums(sc_onehot) > 1]
#         sc_onehot = sc_onehot[sample(nrow(sc_onehot), replace = F), ]
#
#         dmat = proxyC::dist(sc_onehot, method = "manhattan")
#         tr = phangorn::upgma(dmat)
#
#         total_depth = max(node.depth.edgelength(tr))
#         tr$edge.length = tr$edge.length / total_depth * total_time
#         tr = name_nodes(tr)
#         tr
# })
# exp_params$tr5_eval = future_map2(exp_params$tr, exp_params$tr5, function(a, b) {
#         if (!is.null(a) & !is.null(b)) {
#                 return(evalute_gr(b, a))
#         }
# })
# library(furrr)
# plan(multisession, workers = 6)
# total_time = 15.0 - 0.6
# exp_params$gr5 = future_map(exp_params$tr5, function(tr_r) {
#         if (is.null(tr_r)) {
#                 return(NULL)
#         }
#         total_depth = max(node.depth.edgelength(tr_r))
#         tr_r$edge.length = tr_r$edge.length / total_depth * total_time
#         tr_r = name_nodes(tr_r)
#         sc_celltypes = get_type_from_id(tr_r$tip.label)
#         all(names(sc_celltypes) == tr_r$tip.label)
#         gr = reconstruct_graph(tr_r, sc_celltypes, total_time, stat_func = mean, theta = 1.0)
#         gr
# }, .progress = T, .options = furrr_options(seed = T))
# exp_params = mutate(exp_params,
#                     gr5_eval = map2(gr5, big_graph_id, function(x, i) {
#                             if (is.null(x)) {
#                                     return(NULL)
#                             }
#                             evalute_gr(x, as.phylo(all_graphs[[i]]))
#                     }))
#### END UPGMA Hamming (Manhattan) part ####

# Truth Part #
exp_params = mutate(exp_params, gr = future_map(tr, function(tr_r) {
        if (is.null(tr_r)) {
                return(NULL)
        }
        sc_celltypes = get_type_from_id(tr_r$tip.label)
        names(sc_celltypes) = tr_r$tip.label
        gr = reconstruct_graph(tr_r, sc_celltypes, total_time, stat_func = mean, theta = 1.0)
        gr
}, .progress = T, .options = furrr_options(seed = T)))
exp_params = mutate(exp_params,
                    gr_eval = map2(gr, big_graph_id, function(x, i) {
                            if (is.null(x)) {
                                    return(NULL)
                            }
                            evalute_gr(x, as.phylo(all_graphs[[i]]))
                    }))
# End Truth Part #

#### SPS part ####
exp_params$sps = map(1:nrow(exp_params), function(j) {
        message(j)
        if (is.null(exp_params$tr[[j]])) {
                return(NULL)
        }
        shared_progenitor_score(exp_params$tr[[j]],
                                all_graphs[[exp_params$big_graph_id[j]]]$tip_id)
})
exp_params$sps_dist = map(exp_params$sps, function(sps_mat) {
        if (is.null(sps_mat)) {
                return(NULL)
        }
        dmat = 1 - sps_mat/max(sps_mat)
        dmat
})
exp_params$sps_gr = map(exp_params$sps_dist, function(x) {
        if (is.null(x)) {
                return(NULL)
        }
        gr = phangorn::upgma(as.dist(x))
        total_depth = max(node.depth.edgelength(gr))
        gr$edge.length = gr$edge.length / total_depth * total_time #TODO: what to use for sps as max dist?
        gr = name_nodes(gr)
        gr
})
# exp_params$sps_nj_gr = map(exp_params$sps_dist, function(x) {
#         if (is.null(x)) {
#                 return(NULL)
#         }
#         phytools::midpoint.root(phangorn::NJ(as.dist(x)))
# })
exp_params$sps_gr_eval = map2(exp_params$sps_gr, exp_params$big_graph_id,
                              function(x, i) {
                                      if (is.null(x)) {
                                              return(NULL)
                                      }
                                      evalute_gr(x, as.phylo(all_graphs[[i]])
                                      )
                              })
#### end sps part ####


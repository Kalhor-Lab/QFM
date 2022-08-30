init_cell_barcode <- function(mut_param, num_cell = 1, somatic = F) {
        if (somatic) {
                # for more than one cell, creates alleles for each cell
                out = som_barcode_matrix(num_cell)
        } else {
                out = matrix("0", num_cell, length(mut_param))
        }
        out
}
distribute_barcodes <- function(barcodes, n_vec) {
        if (nrow(barcodes) != sum(n_vec)) {
                print(nrow(barcodes))
                print(n_vec)
        }
        assertthat::assert_that(nrow(barcodes) == sum(n_vec))
        labels = rep(seq_along(n_vec), n_vec)
        if (length(labels) > 1) {
                assign = sample(labels, size = length(labels), replace = F)
        } else {
                assign = labels
        }
        indices = lapply(seq_along(n_vec), function(i) {
                which(assign == i)
        })
        barcodes_split = lapply(indices, function(x) {
                barcodes[x, , drop = F]
        })
        barcodes_split
}
double_barcodes <- function(b1, b2) {
        if (nrow(b2) > 0) {
                b2 = b2[rep(1:nrow(b2), each = 2), , drop = F]
        }
        rbind(b1, b2)
}
get_new_muts_barcode <- function(b0, b1) {
        assertthat::assert_that(all(dim(b0) == dim(b1)))
        lapply(1:ncol(b0), function(j) {
                unique(b1[, j][!(b1[, j] %in% b0[, j])])
        })
}
# # update cell barcode table
# simulate_barcodes_distributed <- function(barcodes, c0, c1, c1_dist) {
#         assertthat::assert_that(nrow(barcodes) == sum(c0[[1]]))
#         assertthat::assert_that(nrow(barcodes) == sum(c1[[1]]))
#         new_muts = get_new_muts(c0, c1)
#         if (length(c1_dist[[1]]) > 1) {
#                 n_all = sapply(c1_dist, function(x) sapply(x, sum))
#                 assertthat::assert_that(all(apply(n_all, 1, function(x) all(x == x[1]))))
#                 n_vec = n_all[, 1]
#                 assertthat::assert_that(nrow(barcodes) == sum(n_vec))
#                 labels = rep(seq_along(n_vec), n_vec)
#                 if (length(labels) > 1) {
#                         assign = sample(labels, size = length(labels), replace = F)
#                 } else {
#                         assign = labels
#                 }
#                 indices = lapply(seq_along(n_vec), function(i) {
#                         which(assign == i)
#                 })
#         } else {
#                 indices = list(1:nrow(barcodes))
#         }
#         barcode_dist = lapply(seq_along(indices), function(i) {
#                 barcode_subset = barcodes[indices[[i]], , drop = F]
#                 if (length(indices[[i]]) > 0) {
#                         for (j in 1:length(mut_param)) {
#                                 for (m in new_muts[[j]]) {
#                                         if (m %in% names(c1_dist[[j]][[i]])){
#                                                 m_count = c1_dist[[j]][[i]][m]
#                                                 if (m_count > 0) {
#                                                         u_indices = which(barcode_subset[, j] == "0")
#                                                         # assertthat::assert_that(length(u_indices) >= m_count)
#                                                         if (length(u_indices) < m_count) {
#                                                                 print(barcode_subset[, j])
#                                                                 print(c0[[j]])
#                                                                 print(c1[[j]])
#                                                                 print(u_indices)
#                                                                 print(m_count)
#                                                         }
#                                                         barcode_subset[sample(u_indices, m_count, replace = F), j] = m
#                                                 }
#                                         }
#                                 }
#                         }
#                 }
#                 barcode_subset
#         })
#         return(list(indices = indices,
#                     barcodes = barcode_dist))
# }
# check_same <- function(x, y) {
#         # x = c0[[5]]
#         # y = c1b[[5]]
#         x = x[x>0]
#         y = y[y>0]
#         if (!(length(x) == length(y))) {
#                 print(x[order(names(x))])
#                 print(y[order(names(y))])
#                 return(F)
#         }
#         if(!all(sort(names(x)) == sort(names(y)))) {
#                 return(F)
#         }
#         if(!all(y[match(names(x), names(y))] == x)) {
#                 return(F)
#         }
#         return(T)
# }
# barcodes0 = x$barcode_history[[1]]
# barcodes_update = simulate_barcodes_distributed(barcodes0, c0, c1, c1_dist)
#
# barcodes0 = double_barcodes(barcodes_update$barcodes[[1]],
#                             barcodes_update$barcodes[[2]])
# c1b = lapply(1:ncol(barcodes0), function(j) table(barcodes0[, j]))
# if (!all(mapply(check_same, c0, c1b))) {
#         sapply(1:length(mut_param), function(i) {
#                 message(i)
#                 if (!check_same(c0[[i]], c1b[[i]])) {
#                         print(c0[[i]][order(names(c0[[i]]))])
#                         print(c1b[[i]][order(names(c1b[[i]]))])
#                         stop()
#                 }
#         })
# }
# assertthat::assert_that(all(mapply(check_same, c0, c1b)))
# barcodes_mode = simulate_barcodes_distributed(barcodes0, c_end, c_end_mut, c_mode)
#
# mapply(check_same,
#        lapply(c_mode, "[[", 1),
#        lapply(1:length(mut_param), function(j) table(barcodes_mode$barcodes[[1]][, j])))
# x$end_barcodes = simulate_barcodes_distributed(barcodes0, c_end, c_end_mut, lapply(c_end_mut, list))$barcodes
# b1 = double_barcodes(y$end_barcodes_mode$barcodes[[3]],
#                      y$end_barcodes_mode$barcodes[[1]])
# b2 = double_barcodes(y$end_barcodes_mode$barcodes[[3]],
#                      y$end_barcodes_mode$barcodes[[2]])
#
# x1$barcode_history = list(b1)
# x2$barcode_history = list(b2)
#
# collect_gens[[type_graph$root_id]]$barcode_history = list(init_cell_barcode(mut_param))
# x$end_barcodes_mode = barcodes_mode
sample_sc_mut_history_gens <- function(x, mut_param, somatic = F, target_time=NA) {
        if (somatic) {
                barcode_mut_func = sample_somatic_barcodes
        } else {
                barcode_mut_func = sample_all_barcodes
        }
        # message(x$cell_type)
        # generate history of mutations from "start_barcodes"
        assertthat::assert_that(nrow(x$start_barcodes) == x$sample_size[1])
        b0 = x$start_barcodes
        b_history = list()
        e_history = list()
        n_history = list()
        # m_history = list()
        if (x$num_gen >= 1) {
                for (i in 1:x$num_gen) {
                        b1 = barcode_mut_func(b0,
                                              mut_param = mut_param,
                                              start_time = x$start_time + x$double_time * (i - 1),
                                              end_time = x$start_time + x$double_time * i)
                        # singleton
                        n0 = x$sample_size[i]*2 - x$sample_size[i+1]
                        # non-singleton
                        n1 = x$sample_size[i] - n0
                        # distributed and mutated counts
                        b1_dist = distribute_barcodes(b1, c(n0, n1))
                        b_history = append(b_history, list(b1_dist))
                        # m_history = append(m_history, list(get_new_muts_barcode(b0, b1)))
                        b0 = double_barcodes(b1_dist[[1]],
                                             b1_dist[[2]])
                        # assign cell id
                        rownames(b0) = paste0("type_", x$cell_type, "_gen_", i, "_", 1:nrow(b0))
                        edges = data.frame(in_node = c(rownames(b1_dist[[1]]), rep(rownames(b1_dist[[2]]), each = 2)),
                                           out_node = rownames(b0), stringsAsFactors = F)
                        e_history = append(e_history, list(edges))
                        nodes = data.frame(id = c(rownames(b1_dist[[1]]),
                                                  rownames(b1_dist[[2]])),
                                           degree = c(rep(2, nrow(b1_dist[[1]])),
                                                     rep(3, nrow(b1_dist[[2]]))),
                                           type = x$cell_type,
                                           gen = i-1,
                                           time = x$start_time + i * x$double_time,
                                           double_time = x$double_time,
                                           leaf = FALSE, stringsAsFactors = F)
                        n_history = append(n_history, list(nodes))
                }
        }
        b_end = b0
        m_end_time = ifelse(!x$active & !is.na(target_time),
                            target_time,
                            x$end_time + x$double_time)
        b_end_mut = barcode_mut_func(b_end,
                                     mut_param = mut_param,
                                     start_time = x$end_time,
                                     end_time = m_end_time)

        # m_history = append(m_history, list(get_new_muts(b_end, b_end_mut)))
        x$barcode_history = b_history
        x$edges_history = e_history
        x$nodes_history = n_history
        x$end_barcodes = b_end_mut
        if (x$active) {
                b_mode_dist = distribute_barcodes(b_end_mut, x$end_mode_sample_size)
                x$end_barcodes_mode = b_mode_dist
        }
        x
}
split_barcodes <- function(y, x1, x2, mut_param, target_time, somatic) {
        message(y$cell_type)
        # b1 = double_barcodes(y$end_barcodes_mode[[3]], y$end_barcodes_mode[[1]])
        # b2 = double_barcodes(y$end_barcodes_mode[[3]], y$end_barcodes_mode[[2]])

        # y = gens1$`5`
        # x1 = gens1$`4`
        # x2 = gens1$`-6`
        # rename distributed barcodes to the daughter types
        old_node_names1 = rownames(y$end_barcodes_mode[[1]])
        old_node_names2 = rownames(y$end_barcodes_mode[[2]])
        new_node_names1 = paste0("type_", x1$cell_type, "_gen_-1_", 1:nrow(y$end_barcodes_mode[[1]]))
        new_node_names2 = paste0("type_", x2$cell_type, "_gen_-1_", 1:nrow(y$end_barcodes_mode[[2]]))
        names(new_node_names1) = old_node_names1
        names(new_node_names2) = old_node_names2
        # rename the old nodes to new ones in edge node history
        if (y$num_gen == 0) {
                y$transition_edges$out_node = c(new_node_names1, new_node_names2)[y$transition_edges$out_node]
        } else {
                y$edges_history[[y$num_gen]]$out_node = c(new_node_names1, new_node_names2)[y$edges_history[[y$num_gen]]$out_node]
        }
        rownames(y$end_barcodes_mode[[1]]) = new_node_names1
        rownames(y$end_barcodes_mode[[2]]) = new_node_names2

        # child 1
        b1_n0_m1 = 2 * y$end_mode_sample_size[1] - x1$start_mode_sample_size[1]
        b1_n1_m1 = y$end_mode_sample_size[1] - b1_n0_m1
        b1_m1_dist = distribute_barcodes(y$end_barcodes_mode[[1]],
                                         c(b1_n0_m1, b1_n1_m1))
        b1_m1 = double_barcodes(b1_m1_dist[[1]], b1_m1_dist[[2]])
        # TODO: need to add singletons for assymetric mode
        b1 = b1_m1

        # child 2
        b2_n0_m1 = 2 * y$end_mode_sample_size[2] - x2$start_mode_sample_size[1]
        b2_n1_m1 = y$end_mode_sample_size[2] - b2_n0_m1
        b2_m1_dist = distribute_barcodes(y$end_barcodes_mode[[2]],
                                         c(b2_n0_m1, b2_n1_m1))
        b2_m1 = double_barcodes(b2_m1_dist[[1]], b2_m1_dist[[2]])
        # TODO: need to add singletons for assymetric mode
        b2 = b2_m1
        nodes1 = data.frame(id = c(rownames(b1_m1_dist[[1]]),
                                   rownames(b1_m1_dist[[2]])),
                            degree = c(rep(2, nrow(b1_m1_dist[[1]])),
                                       rep(3, nrow(b1_m1_dist[[2]]))),
                            type = x1$cell_type,
                            gen = -1,
                            time = y$end_time + y$double_time,
                            double_time = y$double_time,
                            leaf = FALSE, stringsAsFactors = F)
        nodes2 = data.frame(id = c(rownames(b2_m1_dist[[1]]),
                                   rownames(b2_m1_dist[[2]])),
                            degree = c(rep(2, nrow(b2_m1_dist[[1]])),
                                       rep(3, nrow(b2_m1_dist[[2]]))),
                            type = x2$cell_type,
                            gen = -1,
                            time = y$end_time + y$double_time,
                            double_time = y$double_time,
                            leaf = FALSE, stringsAsFactors = F)
        rownames(b1) = paste0("type_", x1$cell_type, "_gen_", 0, "_", 1:nrow(b1))
        rownames(b2) = paste0("type_", x2$cell_type, "_gen_", 0, "_", 1:nrow(b2))
        edges1 = data.frame(in_node = c(rownames(b1_m1_dist[[1]]),
                                        rep(rownames(b1_m1_dist[[2]]), each = 2)),
                            out_node = rownames(b1), stringsAsFactors = F)
        edges2 = data.frame(in_node = c(rownames(b2_m1_dist[[1]]),
                                        rep(rownames(b2_m1_dist[[2]]), each = 2)),
                            out_node = rownames(b2), stringsAsFactors = F)

        assertthat::assert_that(nrow(b1) == x1$sample_size[1])
        assertthat::assert_that(nrow(b2) == x2$sample_size[1])
        x1$start_barcodes = b1
        x2$start_barcodes = b2
        # here transition is during differentiation (not unsampled)
        x1$transition_edges = edges1
        x2$transition_edges = edges2
        x1$transition_nodes = nodes1
        x2$transition_nodes = nodes2
        x1 = sample_sc_mut_history_gens(x1, mut_param = mut_param, target_time = target_time, somatic = somatic)
        x2 = sample_sc_mut_history_gens(x2, mut_param = mut_param, target_time = target_time, somatic = somatic)
        return(list(x1, x2, y))
}
simulate_sc_muts <- function(collect_gens, type_graph, mut_param, somatic = F) {
        collect_gens[[type_graph$root_id]]$start_barcodes = init_cell_barcode(mut_param,
                                                                              num_cell = collect_gens[[type_graph$root_id]]$sample_size[1],
                                                                              somatic = somatic)
        rownames(collect_gens[[type_graph$root_id]]$start_barcodes) =
                paste0("type_", type_graph$root_id, "_gen_", 0, "_", 1:nrow(collect_gens[[type_graph$root_id]]$start_barcodes))
        collect_gens[[type_graph$root_id]] = sample_sc_mut_history_gens(collect_gens[[type_graph$root_id]],
                                                                        mut_param = mut_param,
                                                                        target_time = type_graph$target_time,
                                                                        somatic = somatic)
        for (node_id in forward_merge_sequence(type_graph)) {
                node_dau = as.character(type_graph$merge[node_id, ])
                temp = split_barcodes(collect_gens[[node_id]],
                                      collect_gens[[node_dau[1]]],
                                      collect_gens[[node_dau[2]]],
                                      mut_param = mut_param,
                                      target_time = type_graph$target_time,
                                      somatic = somatic)
                collect_gens[[node_dau[1]]] = temp[[1]]
                collect_gens[[node_dau[2]]] = temp[[2]]
                # rewriting the current gen just to update the cell names
                collect_gens[[as.character(node_id)]] = temp[[3]]
        }
        collect_gens
}

# gens1$`5`$start_barcodes = init_cell_barcode(mut_param)
# rownames(gens1$`5`$start_barcodes) = paste0("type_", type_graph$root_id, "_gen_", 0, "_", 1)
# gens1$`5` = sample_sc_mut_history_gens(gens1$`5`,
#                                        mut_param = mut_param)

# type_graph = type_graph0
# mut_param = do.call(make_mut_param_by_rate, mut_p)
# gens0 = make_gens(type_graph = type_graph)
# gens1 = sample_size_gens(gens0, type_graph, sample_size = 5000)
# gens1 = simulate_sc_muts(gens1, type_graph, mut_param = mut_param)
#

# deprecated
# collect_nodes <- function(gens, tip_id) {
#         # gens = gens1
#         nodes_leaf = lapply(gens[tip_id], function(x) {
#                 nodes = data.frame(id = rownames(x$end_barcodes),
#                                    degree = 1,
#                                    gen = x$num_gen+1,
#                                    type = x$cell_type,
#                                    time = x$end_time,
#                                    double_time = x$double_time,
#                                    leaf = TRUE, stringsAsFactors = F)
#                 nodes
#         })
#         nodes_internal = lapply(gens, function(x) {
#                 rbind(x$transition_nodes,
#                       do.call(rbind, x$nodes_history))
#         })
#         rbind(do.call(rbind, nodes_internal),
#               do.call(rbind, nodes_leaf))
# }
# deprecated
# collect_edges <- function(gens) {
#         out = do.call(rbind, lapply(gens, function(x) {
#                 rbind(x$transition_edges,
#                       do.call(rbind, x$edges_history))
#         }))
# }



# bb = calculate_type_counts(type_graph2)[type_graph2$node_id]
# aa = calculate_type_counts(type_graph2)[type_graph2$tip_id]
# sapply(aa, function(z) {
#         m = ceiling(z/bb[3]*16)
#         n = z - m
#         k = 1000
#         dhyper(0, m = m, n = n, k = k)
# })
#
# library(dendextend)
# cutree(tr_upgma, k = 10)






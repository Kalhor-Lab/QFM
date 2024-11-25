make_edge_tb <- function(tr_r) {
        tr_node_depth = node.depth.edgelength(tr_r)
        names(tr_node_depth) = c(tr_r$tip.label, tr_r$node.label)

        edge_mapper = character(max(tr_r$edge))
        edge_mapper[(1:Nnode(tr_r))+Ntip(tr_r)] = tr_r$node.label
        edge_mapper[1:Ntip(tr_r)] = tr_r$tip.label
        edge_df = tibble(from = edge_mapper[tr_r$edge[, 1]],
                         to = edge_mapper[tr_r$edge[, 2]])
        edge_df$length = tr_r$edge.length
        edge_df$from_time = tr_node_depth[edge_df$from]
        edge_df$to_time = tr_node_depth[edge_df$to]
        edge_df
}
make_gr_tr_data <- function(gr, tr, sc_celltypes) {
        # processed input
        # gr
        gr_node_time = node.depth.edgelength(gr)
        names(gr_node_time) = c(gr$tip.label, gr$node.label)
        out_list = list_dd_and_tips(gr)
        gr_tip_list = out_list$tips
        gr_dd = out_list$dd
        gr_tips = map(gr$tip.label, function(x) x)
        names(gr_tips) = gr$tip.label
        gr_tip_list = c(gr_tip_list, gr_tips)
        gr_dd = c(gr_dd, gr_tips)
        gr_tip_size = map_dbl(gr_tip_list, length)
        gr_edges_tb = make_edge_tb(gr)
        # tr
        tr_node_depth = node.depth.edgelength(tr)
        names(tr_node_depth) = c(tr$tip.label, tr$node.label)
        out_list = list_dd_and_tips(tr)
        tr_dd = out_list$dd
        tr_tip_list = out_list$tips
        tr_tip_type_list = map(tr_tip_list, function(x) {
                cell_types = unique(sc_celltypes[x])
                cell_types[!is.na(cell_types)]
        })
        edge_tb = make_edge_tb(tr)
        gr_root_len = max(tr_node_depth) - max(gr_node_time)

        list(gr = gr,
             gr_node_time = gr_node_time,
             gr_tips = gr_tips,
             gr_tip_list = gr_tip_list,
             gr_tip_size = gr_tip_size,
             gr_edges_tb = gr_edges_tb,
             gr_dd = gr_dd,
             gr_root_len = gr_root_len,
             tr = tr,
             sc_celltypes = sc_celltypes,
             tr_node_time = tr_node_depth,
             tr_tip_list = tr_tip_list,
             tr_tip_type_list = tr_tip_type_list,
             tr_dd = tr_dd,
             tr_edges_tb = edge_tb
        )
}
assign_node_states <- function(data_obj) {
        tr_node_assign = map_chr(data_obj$tr$node.label, function(node) {
                x = data_obj$tr_tip_type_list[[node]]
                # naive state assignment
                candid_states = names(data_obj$gr_tip_list)[map_lgl(data_obj$gr_tip_list, function(y) all(x %in% y))]
                candid_sizes = data_obj$gr_tip_size[candid_states]
                assertthat::assert_that(sum(candid_sizes == min(candid_sizes)) == 1)
                temp_state = candid_states[which.min(candid_sizes)]
                temp_state
        })
        names(tr_node_assign) = data_obj$tr$node.label
        tr_node_assign = c(tr_node_assign, data_obj$sc_celltypes)
        tr_node_assign
}
correct_trans_time <- function(gr_trans_time, gr_tr_data) {
        gr_node_dd = gr_tr_data$gr_dd[gr_tr_data$gr$node.label]
        for (x in names(gr_node_dd)) {
                if(gr_trans_time[gr_node_dd[[x]][1]] < gr_trans_time[x]) {
                        gr_trans_time[gr_node_dd[[x]][1]] = gr_trans_time[x] + pmax(1e-5, rnorm(1, mean = 0.001, sd = 0.0001))
                        # message(paste0('corrected: ', gr_node_dd[[x]][1]))
                }
                if(gr_trans_time[gr_node_dd[[x]][2]] < gr_trans_time[x]) {
                        gr_trans_time[gr_node_dd[[x]][2]] = gr_trans_time[x] + pmax(1e-5, rnorm(1, mean = 0.001, sd = 0.0001))
                        # message(paste0('corrected: ', gr_node_dd[[x]][2]))
                }
        }
        return(gr_trans_time)
}
reconstruct_graph <- function(tr, sc_celltypes, total_time, theta, stat_func = mean) {
        transition_tb = get_transitions(tr, sc_celltypes)
        dmat_nondisjoint = transition2dist(transition_tb, stat_func, total_time)
        dmat_disjoint = transition2dist(filter(transition_tb, disjoint), stat_func, total_time)
        type_names = rownames(dmat_nondisjoint)

        dmat = theta * dmat_disjoint[type_names, type_names] +
                (1 - theta) * dmat_nondisjoint[type_names, type_names]
        ord_indices = sample(nrow(dmat), replace = F)
        dmat = dmat[ord_indices, ord_indices]
        gr = phangorn::upgma(dmat)
        # gr = upgma(dmat)
        gr$root.edge = total_time - max(node.depth.edgelength(tr))
        gr
}
get_transitions <- function(tr, sc_celltypes) {
        out_list = list_dd_and_tips(tr)
        tr_dd = out_list$dd
        tr_tip_list = out_list$tips
        tr_tip_type_list = map(tr_tip_list, function(x) {
                cell_types = unique(sc_celltypes[x])
                cell_types[!is.na(cell_types)]
        })
        tr_tip_type_list = c(tr_tip_type_list, as.list(sc_celltypes))
        tr_node_depth = node.depth.edgelength(tr)
        names(tr_node_depth) = c(tr$tip.label, tr$node.label)
        # transitions where both daughter restricted
        transition_tb = map(names(tr_dd), function(x) {
                d1 = tr_dd[[x]][1]
                d2 = tr_dd[[x]][2]
                # t_x = paste0(sort(as.numeric(tr_tip_type_list[[x]]), decreasing = T), collapse = ", ")
                # t_d1 = paste0(sort(as.numeric(tr_tip_type_list[[d1]]), decreasing = T), collapse = ", ")
                # t_d2 = paste0(sort(as.numeric(tr_tip_type_list[[d2]]), decreasing = T), collapse = ", ")
                t_x = paste0(sort(tr_tip_type_list[[x]], decreasing = T), collapse = ", ")
                t_d1 = paste0(sort(tr_tip_type_list[[d1]], decreasing = T), collapse = ", ")
                t_d2 = paste0(sort(tr_tip_type_list[[d2]], decreasing = T), collapse = ", ")
                transition_list = tibble()
                if (t_x != t_d1 & t_x != t_d2) {
                        transition_list = bind_rows(transition_list,
                                                    tibble(in_node = x,
                                                           in_type = t_x,
                                                           out_node1 = d1,
                                                           out_type1 = t_d1,
                                                           out_node2 = d2,
                                                           out_type2 = t_d2,
                                                           time = tr_node_depth[x]))
                }
                return(transition_list)
        }) %>% bind_rows()

        transition_tb$in_type_vec = strsplit(transition_tb$in_type, ", ")
        transition_tb$out_type1_vec = strsplit(transition_tb$out_type1, ", ")
        transition_tb$out_type2_vec = strsplit(transition_tb$out_type2, ", ")

        # check union
        all(map_lgl(1:nrow(transition_tb), function(i) all(sort(transition_tb$in_type_vec[[i]]) ==
                                                                   sort(union(transition_tb$out_type1_vec[[i]],
                                                                              transition_tb$out_type2_vec[[i]])))))
        # check disjoint
        transition_tb = mutate(transition_tb, disjoint = map2_lgl(out_type1_vec, out_type2_vec, function(x, y) {
                all(!x %in% y)
        }))
        transition_tb
}
transition2dist <- function(transition_tb, stat_func, total_time, replace_na = T) {
        # generate pairs
        transition_tb = mutate(transition_tb, pairs = map2(out_type1_vec, out_type2_vec, function(x, y) {
                tidyr::crossing(pair_x = setdiff(x, y),
                                pair_y = setdiff(y, x))
        }))
        transition_tb_pairs = transition_tb %>% unnest(cols = pairs)
        # make pairs unique
        transition_tb_pairs = mutate(transition_tb_pairs,
                                     pair_x_sorted = map2_chr(pair_x, pair_y, function(x, y) {
                                             # c(x, y)[order(as.numeric(c(x, y)))[2]]
                                             sort(c(x, y))[2]
                                     }),
                                     pair_y_sorted = map2_chr(pair_x, pair_y, function(x, y) {
                                             # c(x, y)[order(as.numeric(c(x, y)))[1]]
                                             sort(c(x, y))[1]
                                     }))
        transition_tb_pairs = transition_tb_pairs %>%
                select(-c(disjoint, pair_x, pair_y)) %>%
                nest(data = - c(pair_x_sorted, pair_y_sorted))
        # check no duplicates
        any(map_lgl(transition_tb_pairs$data, function(x) as.logical(anyDuplicated(x$in_node))))

        # mean restriction time
        transition_tb_pairs = mutate(transition_tb_pairs, time = map_dbl(data, function(x) stat_func(x$time)))

        dist_df = transition_tb_pairs %>%
                mutate(dist = map_dbl(time, function(x) {
                        (total_time - x)*2
                })) %>%
                mutate(N1 = pair_x_sorted, N2 = pair_y_sorted, Var1 = pair_x_sorted, Var2 = pair_y_sorted)
        dmat = dist_df2mat(dist_df)
        diag(dmat) = 0
        if (replace_na) {
                dmat[is.na(dmat)] = total_time * 2
        }
        dmat
}
est_transition_time <- function(data_obj, tr_node_assign, stat_func = mean) {
        trans_assign_tb = get_trans_assign_tb(data_obj, tr_node_assign)
        trans_type_time = trans_assign_tb %>% group_by(in_type) %>% summarize(time = stat_func(time))
        trans_time = trans_type_time$time
        names(trans_time) = trans_type_time$in_type
        trans_time
}
get_parent_node_from_edges <- function(node, edges_tb) {
        edges_tb$from[edges_tb$to == node]
}
update_edge_tb_state <- function(data_obj, tr_node_assign) {
        # assign state to edge_tb
        edge_tb = data_obj$tr_edges_tb
        gr_root_len = max(data_obj$tr_node_time) - max(data_obj$gr_node_time)
        edge_tb$from_type = tr_node_assign[edge_tb$from]
        edge_tb$to_type = tr_node_assign[edge_tb$to]
        # assign type path for each edge
        edge_tb = mutate(edge_tb, type_path = map2(from_type, to_type, function(from, to) {
                if (from == to) {
                        return(to)
                } else {
                        state_path = c(to)
                        cur_state = to
                        while(cur_state != from) {
                                cur_state = get_parent_node_from_edges(cur_state, data_obj$gr_edges_tb)
                                state_path = c(state_path, cur_state)
                        }
                        return(state_path)
                }
        }))
        # extending type path upstream
        gr_root_node = data_obj$gr$node.label[1]
        edge_tb$type_path_ext = map(edge_tb$type_path, function(x) {
                cur = x[length(x)]
                out = x
                while (cur != gr_root_node) {
                        cur_new = get_parent_node_from_edges(cur, data_obj$gr_edges_tb)
                        out = c(out, cur_new)
                        cur = cur_new
                }
                out
        })
        data_obj$tr_edges_state_tb = edge_tb
        data_obj
}
#' reconstructs phylogeny with phylotime, and infers quantitative fate map with ice_fase (deprecated)
#' @param tr a phylogenetic tree of the "phylo" type
#' @param sc_celltypes either a named character vector specifying types for each row in the cell_mat or a function to be applied to the rownames of the character matrix
#' @param total_time time at sample collection since barcode activation
#' @param root_time length of the root edge in tr
#' @param theta depracated parameters always set to zero
#' @return Fitted ICE_FASE results
ice_fase <- function(tr,
                     sc_celltypes,
                     total_time,
                     root_time = 0,
                     theta = 0.0,
                     gr = NULL) {
        tr = name_nodes(tr)
        if (is.character(sc_celltypes)) {
                assertthat::assert_that(!is.null(names(sc_celltypes)),
                                        msg = "cell types must be named.")
                assertthat::assert_that(all(tr$tip.label %in% names(sc_celltypes)))
                cell_type_vec = sc_celltypes
        } else {
                assertthat::assert_that(class(sc_celltypes) == "function")
                cell_type_vec = sc_celltypes(tr$tip.label)
                names(cell_type_vec) = rownames(tr$tip.label)
        }
        if (is.null(gr)) {
                gr = name_nodes(reconstruct_graph(tr,
                                                  sc_celltypes = cell_type_vec,
                                                  theta = theta,
                                                  total_time = total_time,
                                                  stat_func = mean))
        }
        data_obj = make_gr_tr_data(gr, tr, cell_type_vec)
        tr_node_assign = assign_node_states(data_obj)

        gr_trans_time = est_transition_time(data_obj, tr_node_assign, stat_func = mean)
        # cases where no node is assigned, assign parent time
        if (any(!gr$node.label %in% names(gr_trans_time))) {
                na_state = gr$node.label[!gr$node.label %in% names(gr_trans_time)]
                assertthat::assert_that(all(!na_state %in% tr_node_assign))
                gr_trans_time = gr_trans_time[gr$node.label]
                names(gr_trans_time) = gr$node.label
                gr_trans_time[na_state] = 0.
        }
        gr_tip_time = rep(total_time, length(unique(cell_type_vec)))
        names(gr_tip_time) = unique(cell_type_vec)
        gr_trans_time = c(gr_trans_time, gr_tip_time)
        gr_trans_time = correct_trans_time(gr_trans_time, data_obj)
        gr_trans_time_root = gr_trans_time + root_time

        data_obj = update_edge_tb_state(data_obj, tr_node_assign)
        gr_node_sizes = get_node_size(data_obj, tr_node_assign, gr_trans_time)

        type_counts = table(sc_celltypes[tr$tip.label])
        gr_all_nodes = gr$node.label
        gr_node_sizes_in = map_dbl(gr_node_sizes[gr_all_nodes], 1)
        gr_node_collect_size = map_dbl(data_obj$gr_tip_list[gr_all_nodes], function(x) {
                sum(type_counts[x])
        })
        ps_cov = gr_node_collect_size / gr_node_sizes_in

        out = list(tr = tr,
                   sc_celltypes = sc_celltypes,
                   total_time = total_time,
                   root_time = root_time,
                   theta = theta,
                   gr = gr,
                   gr_trans_time = gr_trans_time_root,
                   gr_node_sizes = gr_node_sizes,
                   tr_node_assign = tr_node_assign,
                   gr_tr_data = data_obj,
                   ps_cov = ps_cov)
        class(out) = "ice_fase_results"
        out = set_color_palette(out)
        out
}
ice_fase_mod = ice_fase
output_estimates <- function(res) {
        out_tb = tibble(progenitor_state = res$gr$node.label)
        out_tb$commitment_time = res$gr_trans_time[out_tb$progenitor_state]
        out_tb$population_size = map_dbl(res$gr_node_sizes[out_tb$progenitor_state], 1)
        out_tb$downstream_state1 = map_chr(res$gr_node_sizes[out_tb$progenitor_state], function(x) names(x)[2])
        out_tb$downstream_size1 = map_dbl(res$gr_node_sizes[out_tb$progenitor_state], 2)
        out_tb$downstream_state2 = map_chr(res$gr_node_sizes[out_tb$progenitor_state], function(x) names(x)[3])
        out_tb$downstream_size2 = map_dbl(res$gr_node_sizes[out_tb$progenitor_state], 3)
        out_tb$commitment_bias1 = out_tb$downstream_size1 / (out_tb$downstream_size1 + out_tb$downstream_size2)
        out_tb$commitment_bias2 = out_tb$downstream_size2 / (out_tb$downstream_size1 + out_tb$downstream_size2)
        out_tb
}
get_trans_assign_tb <- function(data_obj, tr_node_assign) {
        map(names(data_obj$tr_dd), function(x) {
                tr_dd = data_obj$tr_dd
                d1 = tr_dd[[x]][1]
                d2 = tr_dd[[x]][2]
                t_x = tr_node_assign[x]
                t_d1 = tr_node_assign[d1]
                t_d2 = tr_node_assign[d2]
                transition_list = tibble()
                if (t_x != t_d1 & t_x != t_d2) {
                        d1_depth = 0
                        d2_depth = 0
                        cur_states = t_x
                        while(!t_d1 %in% cur_states) {
                                d1_depth = d1_depth + 1
                                cur_states = reduce(data_obj$gr_dd[cur_states], c)
                                if (all(cur_states %in% data_obj$gr_tips)) {
                                        assertthat::assert_that(t_d1 %in% cur_states)
                                }
                        }
                        cur_states = t_x
                        while(!t_d2 %in% cur_states) {
                                d2_depth = d2_depth + 1
                                cur_states = reduce(data_obj$gr_dd[cur_states], c)
                                if (all(cur_states %in% data_obj$gr_tips)) {
                                        assertthat::assert_that(t_d2 %in% cur_states)
                                }
                        }
                        transition_list = bind_rows(transition_list,
                                                    tibble(in_node = x,
                                                           in_type = t_x,
                                                           out_node1 = d1,
                                                           out_type1 = t_d1,
                                                           out_depth1 = d1_depth,
                                                           out_node2 = d2,
                                                           out_type2 = t_d2,
                                                           out_depth2 = d2_depth,
                                                           time = data_obj$tr_node_time[x],
                                                    ))
                }
                return(transition_list)
        }) %>% bind_rows()
}
plot_ice_times <- function(res) {
        lout = create_layout(as.igraph(res$gr), layout = "dendrogram")
        lout_node = lout %>% filter(!leaf) %>% arrange(desc(x))
        total_time = max(res$gr_trans_time)
        trans_assign_tb = get_trans_assign_tb(res$gr_tr_data, res$tr_node_assign)
        b_time = trans_assign_tb %>% mutate(in_type = factor(in_type, levels = lout_node$name)) %>%
                ggboxplot(y = "time", x = "in_type", fill = "in_type", ylab = "Time", xlab = "", width = 0.8, size = 0.25) +
                scale_fill_manual(values = res$col_pal, guide = F) + ylim(c(total_time, 0)) +
                theme(axis.title.x = element_text(),
                      axis.text.x = element_text(angle = 30))
        b_time
}
plot_node_sizes <- function(res) {
        lout = create_layout(as.igraph(res$gr), layout = "dendrogram")
        lout_node = lout %>% filter(!leaf) %>% arrange(desc(x))

        gr_node_sizes = res$gr_node_sizes
        node_size_tb = bind_rows(map(names(gr_node_sizes), function(x) {
                tibble(node = x,
                       type = c("Progenitor population sizes", rep("Committed split sizes", 2)),
                       split = names(gr_node_sizes[[x]]),
                       count = gr_node_sizes[[x]])
        }))
        node_size_tb$type = factor(node_size_tb$type, levels = c("Progenitor population sizes", "Committed split sizes"))
        b_size = node_size_tb %>% mutate(node = factor(node, levels = lout_node$name)) %>%
                ggbarplot(x = "node", y = "count", fill = "split", ylab = "Count", xlab = "", col = NA) %>%
                facet(facet.by = "type", ncol = 1) + scale_fill_manual(values = res$col_pal, guide = F) +
                theme(text = element_text(size = 12),
                      axis.text.x = element_text(angle = 30))
        b_size
}

# this is not used
update_gr_trans_time <- function(data_obj, gr_trans_time) {
        # assign estimated trans time to edges_tb
        gr_edges_trans = rlang::duplicate(data_obj$gr_edges_tb)
        gr_edges_trans$from_time = gr_trans_time[gr_edges_trans$from]
        gr_edges_trans$to_time = gr_trans_time[gr_edges_trans$to]
        gr_edges_trans$to_time[gr_edges_trans$to %in% data_obj$gr$tip.label] = max(node.depth.edgelength(data_obj$tr))
        gr_edges_trans$length = gr_edges_trans$to_time - gr_edges_trans$from_time
        data_obj$gr_edges_trans = gr_edges_trans
        return(data_obj)
}









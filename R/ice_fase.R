# exp_params
# big_graph_list
# tr_col = "tr3"
# gr_col = "gr3"
# gr_eval_col = "gr3_eval"
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
fate_to_str <- function(x) {
        paste0(sort(as.numeric(x), decreasing = T), collapse = "_")
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


est_transition_time <- function(data_obj, tr_node_assign, stat_func = median) {
        trans_assign_tb = get_trans_assign_tb(data_obj, tr_node_assign)
        trans_type_time = trans_assign_tb %>% group_by(in_type) %>% summarize(time = stat_func(time))
        trans_time = trans_type_time$time
        names(trans_time) = trans_type_time$in_type
        trans_time
}

# data_obj = update_edge_tb_state(data_obj, tr3_assigned_states)
# x = "Node-3"
# trans_time = gr3_trans_time
# tr_node_assign = tr3_assigned_states

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

get_edge_diff <- function(x, data_obj, trans_time) {
        diff_time = trans_time[x]
        gr_tip_list = data_obj$gr_tip_list
        gr_dd = data_obj$gr_dd
        assertthat::assert_that(!is.null(data_obj$tr_edges_state_tb))
        edge_tb = data_obj$tr_edges_state_tb

        state_include = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[x]]))]
        state_include_d1 = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[gr_dd[[x]][1]]]))]
        state_include_d2 = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[gr_dd[[x]][2]]]))]

        edge_tb_diff = mutate(edge_tb, ind = map_lgl(type_path, function(y) any(state_include %in% y)))
        edge_tb_diff = edge_tb_diff %>% filter(from_time <= diff_time &
                                                       to_time > diff_time &
                                                       ind)
        edge_tb_diff = mutate(edge_tb_diff, d1_ind = map_lgl(type_path, function(y) any(state_include_d1 %in% y))) # old version gr_dd[[x]][1] (old version has been overwritten)
        edge_tb_diff = mutate(edge_tb_diff, d2_ind = map_lgl(type_path, function(y) any(state_include_d2 %in% y)))
        assertthat::assert_that(all(!(edge_tb_diff$d1_ind & edge_tb_diff$d2_ind)))
        edge_tb_diff
}

get_node_size <- function(data_obj, tr_node_assign, trans_time) {
        # data_obj = update_edge_tb_state(data_obj, tr_node_assign)
        node_size = map(data_obj$gr$node.label, function(x) {
                edge_diff = get_edge_diff(x, data_obj, trans_time)
                n0 = length(unique(edge_diff$from))
                n1 = sum(edge_diff$d1_ind)
                n2 = sum(edge_diff$d2_ind)
                out = c(n0, n1, n2)
                names(out) = c(x, data_obj$gr_dd[[x]])
                out
        })
        names(node_size) = data_obj$gr$node.label
        node_size
}

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
                                             c(x, y)[order(as.numeric(c(x, y)))[2]]
                                     }),
                                     pair_y_sorted = map2_chr(pair_x, pair_y, function(x, y) {
                                             c(x, y)[order(as.numeric(c(x, y)))[1]]
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
get_parent_node_from_edges <- function(node, edges_tb) {
        edges_tb$from[edges_tb$to == node]
}
filter_noise_spacers <- function(tb, noise_factor = 4) {
        if (nrow(tb) > 1) {
                tb = arrange(tb, desc(count))
                tb_count = tb$count
                if (tb_count[1] > sum(tb_count[2:length(tb_count)]) * noise_factor) {
                        return(tb[1, ])
                } else {
                        return(tb)
                }
        } else {
                return(tb)
        }
}
correct_trans_time <- function(gr_trans_time, gr_tr_data) {
        gr_node_dd = gr_tr_data$gr_dd[gr_tr_data$gr$node.label]
        for (x in names(gr_node_dd)) {
                if(gr_trans_time[gr_node_dd[[x]][1]] < gr_trans_time[x]) {
                        gr_trans_time[gr_node_dd[[x]][1]] = gr_trans_time[x] + rnorm(1, mean = 0.001, sd = 0.0001)
                        # message(paste0('corrected: ', gr_node_dd[[x]][1]))
                }
                if(gr_trans_time[gr_node_dd[[x]][2]] < gr_trans_time[x]) {
                        gr_trans_time[gr_node_dd[[x]][2]] = gr_trans_time[x] + rnorm(1, mean = 0.001, sd = 0.0001)
                        # message(paste0('corrected: ', gr_node_dd[[x]][2]))
                }
        }
        return(gr_trans_time)
}
#' reconstructs phylogeny with phylotime, and infers quantitative fate map with ice_fase
#' @param cell_mat character matrix of single cell lineage barcodes cell by barcoding sites
#' @param sc_celltypes either a named character vector specifying types for each row in the cell_mat or a function to be applied to the rownames of the character matrix
#' @param total_time time at sample collection since barcode activation
#' @return Fitted ICE_FASE results
ice_fase <- function(cell_mat,
                     sc_celltypes,
                     total_time,
                     root_time = 0,
                     theta = 0.0,
                     nrounds = 20,
                     max_depth = 4,
                     remove_uniform = F,
                     time_stat_func = mean) {
        if (remove_uniform) {
                cell_mat = remove_uniform_id(cell_mat)
        }
        mut_p = estimate_mut_p(cell_mat, total_time)
        mat_im = impute_characters(cell_mat, nrounds = nrounds, max_depth = max_depth)

        tr_upgma = phylotime(mat_im[sample(nrow(mat_im)), ],
                             mut_p = mut_p,
                             total_time)
        if (is.character(sc_celltypes)) {
                assertthat::assert_that(!is.null(names(sc_celltypes)),
                                        msg = "cell types must be named.")
                assertthat::assert_that(all(rownames(cell_mat) %in% names(sc_celltypes)))
                cell_type_vec = sc_celltypes
        } else {
                assertthat::assert_that(class(sc_celltypes) == "function")
                cell_type_vec = sc_celltypes(rownames(mat_im))
                names(cell_type_vec) = rownames(mat_im)
        }
        tr = name_nodes(tr_upgma)
        gr = name_nodes(reconstruct_graph(tr,
                                          sc_celltypes = cell_type_vec,
                                          theta = theta,
                                          total_time = total_time,
                                          stat_func = mean))

        data_obj = make_gr_tr_data(gr, tr, cell_type_vec)
        tr_node_assign = assign_node_states(data_obj)

        gr_trans_time = est_transition_time(data_obj, tr_node_assign, stat_func = time_stat_func)
        gr_tip_time = rep(total_time, length(unique(cell_type_vec)))
        names(gr_tip_time) = unique(cell_type_vec)
        gr_trans_time = c(gr_trans_time, gr_tip_time)
        gr_trans_time = correct_trans_time(gr_trans_time, data_obj)
        gr_trans_time = gr_trans_time + root_time

        data_obj = update_edge_tb_state(data_obj, tr_node_assign)
        gr_node_sizes = get_node_size(data_obj, tr_node_assign, gr_trans_time)

        out = list(mat = cell_mat,
                   impuated = mat_im,
                   gr = gr,
                   gr_trans_time = gr_trans_time,
                   gr_node_sizes = gr_node_sizes,
                   tr_node_assign = tr_node_assign,
                   gr_tr_data = data_obj)
        class(out) = "ice_fase_results"
        out
}
#' plot ice_fase results
#' @param out_data "ice_fase_results" object
plot_gr_data <- function(out_data, target_time, gr_col = NULL, gr_lab = NULL) {
        if (is.null(gr_col)) {
                gr_col = gr_color(out_data$gr)
        }
        plot_gr(out_data$gr,
                out_data$gr_trans_time,
                gr_col,
                target_time = target_time,
                node_label = gr_lab)
}
#' Impute missing characters from cell allele matrix
#' @param mat character matrix with rows corresponding to cells and columns corresponding to barcoding sites
options(na.action='na.pass')
impute_characters <- function(mat, nrounds = 100, max_depth = 4) {
        # separate out part of the matrix that has no variability
        col_div = apply(mat, 2, function(x) length(unique(x[!is.na(x)])))
        im_indices = which(col_div > 1)
        keep_indices = which(col_div == 1)

        mat_keep = rlang::duplicate(mat[, keep_indices, drop =F])
        if (ncol(mat_keep) > 0) {
                for (i in 1:ncol(mat_keep)) {
                        char_vec = mat_keep[, i]
                        char_vec_unique = unique(char_vec[!is.na(char_vec)])
                        assertthat::assert_that(length(char_vec_unique) == 1)
                        mat_keep[is.na(char_vec), i] = char_vec_unique
                }
        }
        assertthat::assert_that(all(!is.na(mat_keep)))

        # impute identical elements with the same allele
        mat_im = as.data.frame(rlang::duplicate(mat[, im_indices, drop =F]))
        element_id = 1:ncol(mat_im)
        na_frac = map_dbl(element_id, function(j) {
                mean(is.na(mat[, j]))
        })
        # imputation order
        element_id = element_id[order(na_frac)]
        for (i in element_id) {
                # message(i)
                design_mat = model.matrix(~., mat_im[, -i])[, -1]
                y_fac = factor(as.matrix(mat_im[i])[, 1])
                y_vec = as.numeric(y_fac)
                suppressWarnings(bst <- xgboost(data = design_mat[!is.na(y_vec), ],
                                                label = y_vec[!is.na(y_vec)] - 1,
                                                max_depth = max_depth,
                                                eta = 1,
                                                nthread = 12,
                                                nrounds = nrounds,
                                                objective = "multi:softmax",
                                                num_class= max(y_vec[!is.na(y_vec)]),
                                                verbose = 0))
                y_pred = levels(y_fac)[predict(bst, design_mat[is.na(y_vec), ]) + 1]
                mat_im[is.na(y_vec), i] = y_pred
        }
        mat_im = as.matrix(mat_im)
        rownames(mat_im) = rownames(mat)

        mat_out = cbind(mat_im, mat_keep)
        mat_out = mat_out[, colnames(mat)]
        mat_out
}
remove_uniform_id <- function(mat, abund_thres = 0) {
        # count diversity ignoring single allele
        id_diversity = apply(mat, 2, function(x) {
                x_tab = table(x[!is.na(x)])
                x_tab = x_tab[x_tab > abund_thres]
                return(length(x_tab))
        })
        mat[, id_diversity > 1]
}
filter_cells <- function(mat, cutoff = 2, abund_thres = 1) {
        # make a list of informative spacers
        spacer_ind_mat = do.call(cbind, map(1:ncol(mat), function(j) {
                x = mat[, j]
                x_tab = table(x[!is.na(x)])
                x_tab = x_tab[x_tab > abund_thres]
                if (length(x_tab) <= 1) {
                        return(rep(F, nrow(mat)))
                }
                allele_freq = sort(table(mat[, j]))
                out_spacers = names(allele_freq)[allele_freq > 1]
                out_spacers = out_spacers[out_spacers != "0"]
                mat[, j] %in% out_spacers
        }))
        mat[rowSums(spacer_ind_mat) > cutoff, ]
}

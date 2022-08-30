#' ICE-FASE inference when topology has multifurcations
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
        assertthat::assert_that(all(edge_tb_diff$d1_ind + edge_tb_diff$d2_ind <= 1))
        if (length(gr_dd[[x]]) == 3) {
                state_include_d3 = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[gr_dd[[x]][3]]]))]
                edge_tb_diff = mutate(edge_tb_diff, d3_ind = map_lgl(type_path, function(y) any(state_include_d3 %in% y)))
                assertthat::assert_that(all(edge_tb_diff$d1_ind + edge_tb_diff$d2_ind + edge_tb_diff$d3_ind <= 1))
        }
        edge_tb_diff
}
get_node_size <- function(data_obj, tr_node_assign, trans_time) {
        # data_obj = update_edge_tb_state(data_obj, tr_node_assign)
        node_size = map(data_obj$gr$node.label, function(x) {
                edge_diff = get_edge_diff(x, data_obj, trans_time)
                n0 = length(unique(edge_diff$from))
                n1 = sum(edge_diff$d1_ind)
                n2 = sum(edge_diff$d2_ind)
                if (length(data_obj$gr_dd[[x]]) == 3) {
                        n3 = sum(edge_diff$d3_ind)
                        out = c(n0, n1, n2, n3)
                } else {
                        out = c(n0, n1, n2)
                }
                names(out) = c(x, data_obj$gr_dd[[x]])
                out
        })
        names(node_size) = data_obj$gr$node.label
        node_size
}
make_gr_tr_data_mod2 <- function(gr, tr, sc_celltypes) {
        # processed input
        # gr
        gr_node_time = node.depth.edgelength(gr)
        names(gr_node_time) = c(gr$tip.label, gr$node.label)
        out_list = list_dd_and_tips_mod2(gr)
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
ice_fase_mutlifurc <- function(tr,
                               sc_celltypes,
                               total_time,
                               root_time = 0,
                               theta = 0.0,
                               gr = NULL) {
        # nrounds = 20,
        # max_depth = 4,
        # remove_uniform = F,
        # time_stat_func = mean) {
        # if (remove_uniform) {
        #         cell_mat = remove_uniform_id(cell_mat)
        # }
        # mut_p = estimate_mut_p(cell_mat, total_time)
        # mat_im = impute_characters(cell_mat, nrounds = nrounds, max_depth = max_depth)
        #
        # tr_upgma = phylotime(mat_im[sample(nrow(mat_im)), ],
        #                      mut_p = mut_p,
        #                      total_time)
        # tr = name_nodes(tr)
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
        # if (is.null(gr)) {
        #         gr = name_nodes(reconstruct_graph(tr,
        #                                           sc_celltypes = cell_type_vec,
        #                                           theta = theta,
        #                                           total_time = total_time,
        #                                           stat_func = mean))
        # }
        data_obj = make_gr_tr_data_mod2(gr, tr, cell_type_vec)
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
# replace gr edge length with those from trans_time
correct_gr_edge_len <- function(res) {
        res$gr_tr_data$gr_edges_tb$from_time = res$gr_trans_time[res$gr_tr_data$gr_edges_tb$from]
        res$gr_tr_data$gr_edges_tb$to_time = res$gr_trans_time[res$gr_tr_data$gr_edges_tb$to]
        res$gr_tr_data$gr_edges_tb$length = res$gr_tr_data$gr_edges_tb$to_time - res$gr_tr_data$gr_edges_tb$from_time
        res$gr$edge.length = res$gr_tr_data$gr_edges_tb$length
        res
}
#' collaspe nodes below certain threshold
collaspe_di <- function(res, tol) {
        res = correct_gr_edge_len(res)
        res$gr = di2multi(res$gr, tol = tol)
        res
}

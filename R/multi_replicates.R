# average dmat
get_dmat <- function(tr, sc_celltypes, total_time, theta, stat_func = mean) {
        all_celltypes = sort(unique(sc_celltypes))

        transition_tb = get_transitions(tr, sc_celltypes)
        dmat_nondisjoint = transition2dist(transition_tb, stat_func, total_time)
        dmat_disjoint = transition2dist(filter(transition_tb, disjoint), stat_func, total_time)
        type_names = rownames(dmat_nondisjoint)

        dmat = theta * dmat_disjoint[type_names, type_names] +
                (1 - theta) * dmat_nondisjoint[type_names, type_names]
        dmat[all_celltypes, all_celltypes]
}
reconstruct_graph_multi <- function(tr_list, sc_celltypes_list, total_time, theta, stat_func = mean) {
        dmat_list = map2(tr_list, sc_celltypes_list, function(tr_r, sc_celltypes) {
                get_dmat(tr_r,
                         sc_celltypes,
                         total_time = total_time,
                         theta = theta,
                         stat_func  = stat_func)
        })
        tr_total_len = mean(map_dbl(tr_list, function(tr_r) {
                max(node.depth.edgelength(tr_r))
        }))
        dmat = purrr::reduce(dmat_list, `+`)
        # shuffle order before applying clustering
        ord_indices = sample(nrow(dmat), replace = F)
        dmat = dmat[ord_indices, ord_indices]
        gr = phangorn::upgma(dmat)
        # gr = upgma(dmat)
        gr$root.edge = total_time - tr_total_len
        name_nodes(gr)
}
# TODO: alternatively, average within each tree first
est_transition_time_multi <- function(data_obj_list, stat_func = mean) {
        trans_assign_all = bind_rows(map(data_obj_list, function(out) {
                get_trans_assign_tb(out$gr_tr_data, out$tr_node_assign)
        }))
        trans_type_time = trans_assign_all %>% group_by(in_type) %>%
                summarize(time = stat_func(time))
        trans_time = trans_type_time$time
        names(trans_time) = trans_type_time$in_type
        trans_time
}
merge_node_size <- function(x, y) {
        map2(x, y, function(a, b) {
                a + b
        })
}
get_node_size_multi = function (data_obj_list, gr_trans_time) {
        node_sizes_list = map(data_obj_list, function(out) {
                get_node_size(out$gr_tr_data, out$tr_node_assign, trans_time = gr_trans_time)
        })
        node_sizes_merged = purrr::reduce(node_sizes_list, merge_node_size)
        node_sizes_merged = map(node_sizes_merged, function(x) {
                x / length(node_sizes_list)
        })
        node_sizes_merged
}
ice_fase_multi <- function(tr_list, sc_celltypes_list,
                           total_time,
                           root_time = 0,
                           theta = 0.0,
                           gr = NULL) {
        all_celltypes = sort(unique(sc_celltypes_list[[1]]))
        if (is.null(gr)) {
                gr = reconstruct_graph_multi(tr_list, sc_celltypes_list, total_time = total_time, theta = theta)
        }
        data_obj_list = map2(tr_list, sc_celltypes_list, function(tr, sc_celltypes) {
                cell_type_vec = sc_celltypes
                gr_tr_data = make_gr_tr_data(gr, tr, cell_type_vec)
                tr_node_assign = assign_node_states(gr_tr_data)
                gr_tr_data = update_edge_tb_state(gr_tr_data, tr_node_assign)
                list(gr_tr_data = gr_tr_data,
                     tr_node_assign = tr_node_assign)
        })
        gr_trans_time = est_transition_time_multi(data_obj_list)

        # cases where no node is assigned, assign parent time
        if (any(!gr$node.label %in% names(gr_trans_time))) {
                na_state = gr$node.label[!gr$node.label %in% names(gr_trans_time)]
                gr_trans_time = gr_trans_time[gr$node.label]
                names(gr_trans_time) = gr$node.label
                gr_trans_time[na_state] = 0.
        }
        gr_tip_time = rep(total_time, length(all_celltypes))
        names(gr_tip_time) = all_celltypes
        gr_trans_time = c(gr_trans_time, gr_tip_time)
        gr_trans_time = correct_trans_time(gr_trans_time, data_obj_list[[1]])
        gr_trans_time_root = gr_trans_time + root_time

        gr_node_sizes = get_node_size_multi(data_obj_list, gr_trans_time)

        out = list(tr = tr_list,
                   sc_celltypes = sc_celltypes_list,
                   total_time = total_time,
                   root_time = root_time,
                   theta = theta,
                   gr = gr,
                   gr_trans_time = gr_trans_time_root,
                   gr_node_sizes = gr_node_sizes,
                   gr_tr_data = data_obj_list[[1]]$gr_tr_data # TODO: placeholder for evaluations
        )
        class(out) = "ice_fase_results"
        out = set_color_palette(out)
        out
}

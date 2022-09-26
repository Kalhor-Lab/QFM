# updated evaluation of ice_fase results
# temporary function to get tip size
get_tip_sample_size <- function(type_graph, sample_size) {
        map_dbl(sample_size_gens_mod2(make_gens_mod2(type_graph),
                                      type_graph = type_graph, sample_size),
                function(x) x$sample_size[1])[type_graph$tip_id]
}
process_split_ratio <- function(eval_tb) {
        # node split is the reference for ordering left and right split
        out_tb = mutate(eval_tb,
                        node_split_order = map_dbl(node_split, function(x) {
                                if (is.null(x)) {
                                        return(NA)
                                }
                                return(x[1])
                        }))
        # mapping sampled ratio to the same order
        out_tb = mutate(out_tb, node_split_sampled_order = map2_dbl(node_split, node_split_sampled, function(x, y) {
                if (is.null(x)) {
                        return(NA)
                }
                return(y[names(x)[1]]/sum(y))
        }))
        out_tb = mutate(out_tb, gr_node_split_order = map2_dbl(node_split, gr_node_size_mapped, function(x, y) {
                if (is.null(y)) {
                        return(NA)
                }
                if (is.null(x)) {
                        return(NA)
                }
                y = y[2:3]
                return((1+y[names(x)[1]])/(2+sum(y)))
        }))
        out_tb
}
process_node_size <- function(eval_tb) {
        out_tb = mutate(eval_tb,
                        gr_node_size_in = map_dbl(gr_node_size, function(x) {
                                if (is.null(x)) {
                                        return(NA)
                                }
                                return(pmax(1, x[1]))
                        }))
        out_tb
}
make_true_gr_data_mod2 <- function(type_graph) {
        gr_tips = type_graph$tip_id
        names(gr_tips) = gr_tips
        gr_true_phy = as.phylo_mod2.type_graph(type_graph)
        gr_true_phy$node.label = str_replace_all(gr_true_phy$node.label, ":", "")
        gr_tip_list_true = list_dd_and_tips_mod2(gr_true_phy)$tips
        gr_tip_list_true = c(gr_tip_list_true, gr_tips)
        gr_dd_true = list_dd_and_tips_mod2(gr_true_phy)$dd

        # TODO: this may not be the correct time
        gr_node_time_truth = node.depth.edgelength(gr_true_phy)
        names(gr_node_time_truth) = c(gr_true_phy$tip.label, gr_true_phy$node.label)

        list(gr = gr_true_phy,
             gr_node_time = gr_node_time_truth,
             gr_tip_list = gr_tip_list_true,
             gr_dd = gr_dd_true)
}
get_node_mapping_mod2 <- function(data_obj, type_graph) {
        gr = data_obj$gr
        true_data_obj = make_true_gr_data_mod2(type_graph)

        in_fates = map_chr(data_obj$gr_tip_list[1:Nnode(gr)], fate_to_str)
        in_fates_true = map_chr(true_data_obj$gr_tip_list[1:Nnode(gr)], fate_to_str)
        in_fates_node = names(in_fates)[in_fates %in% in_fates_true]
        matched_indices = match(in_fates[in_fates %in% in_fates_true], in_fates_true)
        in_fates_true_node = names(in_fates_true)[matched_indices]
        # in_fates_node and in_fates_true_node are the matched nodes by in_fates
        # now checking out_fates
        out_fates = map(data_obj$gr_dd[in_fates_node], function(x) {
                out = c(fate_to_str(data_obj$gr_tip_list[[x[1]]]),
                        fate_to_str(data_obj$gr_tip_list[[x[2]]]))
                names(out) = x
                out
        })
        out_fates_true = map(true_data_obj$gr_dd[in_fates_true_node], function(x) {
                out = c(fate_to_str(true_data_obj$gr_tip_list[[x[1]]]),
                        fate_to_str(true_data_obj$gr_tip_list[[x[2]]]))
                names(out) = x
                out
        })
        out_fates_same = map2_lgl(out_fates, out_fates_true, function(x, y) {
                all(x %in% y)
        })
        # make a final vector with matched nodes from gr and truth
        gr_node = in_fates_node[out_fates_same]
        gr_node_mapped = in_fates_true_node[out_fates_same]
        gr_node_all = names(in_fates)
        gr_node_mapped_all = gr_node_mapped[match(gr_node_all, gr_node)]
        # reverse matched gr and truth
        # gr_node_true = names(in_fates_true)
        # gr_node_mapped_rev = gr_node[match(gr_node_true, gr_node_mapped)]
        names(gr_node_mapped_all) = gr_node_all
        gr_node_mapped_all
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
generate_control <- function(exp_tb, tree_panel, num_sim = 100) {
        exp_tb$gr_control_eval = future_map(exp_tb$big_graph_id, function(graph_id) {
                gr = tree_panel$gr[[graph_id]]
                gr_control_tb = bind_rows(map(1:num_sim, function(i) {
                        phy = rtree(length(gr$tip.label))
                        phy$tip.label = gr$tip.label
                        as_tibble(evalute_gr(gr, phy))
                })) %>% summarise(rf = mean(rf),
                                  kc0 = mean(kc0),
                                  kc1 = mean(kc1),
                                  treevec_cor = mean(treevec_cor))
                as.list(gr_control_tb)
        }, .progress = T, .options = furrr_options(seed=TRUE))
        exp_tb
}
evaluate_qfm_v1 <- function(exp_obj, type_graph, true_sampled_sizes, delay = F) {
        # exp_obj is the ice_fase result now
        tr = exp_obj$tr
        gr = name_nodes(exp_obj$gr)
        sc_celltypes = exp_obj$sc_celltypes
        gr_tr_data = exp_obj$gr_tr_data
        # tr_node_assign = exp_obj$tr_node_assign # does not get used
        gr_trans_time = exp_obj$gr_trans_time
        gr_node_sizes = exp_obj$gr_node_sizes
        gr_eval = exp_obj$gr_eval

        if (is.null(tr)) {
                return(list(recon = NULL,
                            true = NULL))
        }
        #### truth ####
        gr_data_true = make_true_gr_data_mod2(type_graph)
        tip_collect_size = get_tip_sample_size(type_graph, map_dbl(true_sampled_sizes[type_graph$tip_id], 1))
        gr_dd_collect_size = map(gr_data_true$gr_dd, function(x) {
                s1 = ifelse(x[1] %in% gr_data_true$gr$tip.label, tip_collect_size[x[1]],
                            sum(tip_collect_size[gr_data_true$gr_tip_list[[x[1]]]]))
                s2 = ifelse(x[2] %in% gr_data_true$gr$tip.label, tip_collect_size[x[2]],
                            sum(tip_collect_size[gr_data_true$gr_tip_list[[x[2]]]]))
                out = c(s1, s2)
                names(out) = x
                out
        })

        node_time = get_true_diff_time_mod3(type_graph, delay = delay)
        node_size = get_true_size_mod2(type_graph, delay = delay)
        node_size_in = map_dbl(node_size, 1)
        node_split_ratio = map(node_size[type_graph$node_id], function(x) {
                x[2:3] / sum(x[2:3])
        })
        node_split_sampled = map(true_sampled_sizes[type_graph$node_id], function(x) {
                x[2:3] / sum(x[2:3])
        })
        #### end truth ####

        # mapped truth for each node in reconstructed
        gr_node_mapped_all = get_node_mapping_mod2(gr_tr_data, type_graph)
        # reconstructed node for each truth (reverse)
        gr_node_true = names(gr_data_true$gr_tip_list[1:Nnode(gr)])
        gr_node_true_mapped = names(gr_node_mapped_all)[match(gr_node_true, gr_node_mapped_all)]
        names(gr_node_true_mapped) = gr_node_true

        gr_node_size_mapped = map(gr_node_sizes, function(x) {
                names(x) = c(gr_node_mapped_all, gr_tr_data$gr_tips)[names(x)]
                x
        })
        node_size_sampled_true = map_dbl(true_sampled_sizes, 1)

        gr_node_all = names(gr_node_mapped_all)
        node_collect_size = map(gr_data_true$gr_tip_list, function(x) {
                purrr::reduce(tip_collect_size[x], `+`)
        })
        names(node_collect_size) = names(gr_data_true$gr_tip_list)

        out_tb_recon = tibble(# gr_eval
                kc0 = gr_eval$kc0,
                kc1 = gr_eval$kc1,
                rf = gr_eval$rf,
                # reconstructed
                node_gr = gr_node_all,
                # gr_time_est = gr_tr_data$gr_node_time[gr_node_all] + gr_tr_data$gr_root_len, # deprecated
                gr_time_trans = gr_trans_time[gr_node_all],
                gr_node_size = gr_node_sizes[gr_node_all],
                gr_node_size_mapped = gr_node_size_mapped[gr_node_all],
                # gr_node_mean_final = node_mean_final[gr_node_all],
                # gr_node_short_edge = frac_short_edge[gr_node_all],
                # mapped truth
                node = gr_node_mapped_all,
                node_time = node_time[gr_node_mapped_all],
                node_size = node_size_in[gr_node_mapped_all],
                node_size_sampled = node_size_sampled_true[gr_node_mapped_all],
                node_split = node_split_ratio[gr_node_mapped_all],
                node_split_sampled = node_split_sampled[gr_node_mapped_all],
                node_dd_collect = gr_dd_collect_size[gr_node_mapped_all],
                node_size_collect = node_collect_size[gr_node_mapped_all]
        )
        # comparison based on truth
        out_tb_true = tibble(# gr_eval
                kc0 = gr_eval$kc0,
                kc1 = gr_eval$kc1,
                rf = gr_eval$rf,
                # truth
                node = gr_node_true,
                node_time = node_time[gr_node_true],
                node_size = node_size_in[gr_node_true],
                node_size_sampled = node_size_sampled_true[gr_node_true],
                node_split = node_split_ratio[gr_node_true],
                node_split_sampled = node_split_sampled[gr_node_true],
                node_dd_collect = gr_dd_collect_size[gr_node_true],
                node_size_collect = node_collect_size[gr_node_true],
                # mapped reconstructed
                node_gr = gr_node_true_mapped,
                # gr_time_est = gr_tr_data$gr_node_time[gr_node_true_mapped] + gr_tr_data$gr_root_len, # deprecated
                gr_time_trans = gr_trans_time[gr_node_true_mapped],
                gr_node_size = gr_node_sizes[gr_node_true_mapped],
                gr_node_size_mapped = gr_node_size_mapped[gr_node_true_mapped],
                # gr_node_mean_final = node_mean_final[gr_node_true_mapped],
                # gr_node_short_edge = frac_short_edge[gr_node_true_mapped]
        )
        out_tb_true = mutate(out_tb_true, is_resolved = as.numeric(!is.na(node_gr)))
        out_tb_true = process_node_size(out_tb_true)
        out_tb_true = process_split_ratio(out_tb_true)
        out_tb_true = mutate(out_tb_true, log2_node_size = log2(node_size))
        out_tb_true = mutate(out_tb_true, gr_time_trans_error = gr_time_trans - node_time)
        # out_tb_true = mutate(out_tb_true, gr_time_est_error = gr_time_est - node_time)
        out_tb_true = mutate(out_tb_true, gr_node_size_logfc = log2(gr_node_size_in / node_size))
        out_tb_true = mutate(out_tb_true, log2_node_sampled = log2(node_size_sampled / node_size))

        out_tb_recon = mutate(out_tb_recon, is_resolved = !is.na(node))
        out_tb_recon = process_node_size(out_tb_recon)
        out_tb_recon = process_split_ratio(out_tb_recon)
        out_tb_recon = mutate(out_tb_recon, log2_node_size = log2(node_size))
        out_tb_recon = mutate(out_tb_recon, gr_time_trans_error = gr_time_trans - node_time)
        # out_tb_recon = mutate(out_tb_recon, gr_time_est_error = gr_time_est - node_time)
        out_tb_recon = mutate(out_tb_recon, gr_node_size_logfc = log2(gr_node_size_in / node_size))
        out_tb_recon = mutate(out_tb_recon, log2_node_sampled = log2(node_size_sampled / node_size))

        list(recon = out_tb_recon,
             true = out_tb_true)
}

perturb_celltypes <- function(sc_celltypes, frac = 0.05, min_keep = 10) {
        if (frac == 0) {
                return(sc_celltypes)
        }
        sc_celltype_counts = table(sc_celltypes)
        perturb_size = pmin(ceiling(sc_celltype_counts * frac), pmax(sc_celltype_counts - min_keep, 0))
        all(perturb_size <= sc_celltype_counts)
        perturbed = unlist(map(names(perturb_size), function(x) {
                ct_exlucde = sc_celltype_counts[names(sc_celltype_counts) != x]
                out = sc_celltypes[sc_celltypes == x]
                out[sample(1:length(out), size = perturb_size[x], replace = F)] =
                        sample(names(ct_exlucde),
                               size = perturb_size[x],
                               prob = ct_exlucde / sum(ct_exlucde),
                               replace = T)
                out
        }))
        perturbed
}
evaluate_gr_res <- function(res, gr_true) {
        res$gr_eval = evalute_gr(res$gr, gr_true)
        res
}
ice_fase_mod <- function(tr,
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
set_color_palette <- function(res, palette = NULL) {
        if (is.null(palette)) {
                res$col_pal = gr_color_v1(res$gr)
        } else {
                assertthat::assert_that(
                        length(palette) == length(c(res$gr$node.label, res$gr$tip.label))
                )
                names(palette) = c(res$gr$node.label, res$gr$tip.label)
                res$col_pal = palette
        }
        res
}
plot_gr_data_mod1 <- function(out_data, target_time, gr_col = out_data$col_pal, gr_lab = NULL) {
        if (is.null(gr_col)) {
                gr_col = gr_color_v1(out_data$gr)
        }
        plot_gr(out_data$gr,
                out_data$gr_trans_time,
                total_time = target_time,
                type_col = gr_col,
                node_label = gr_lab)
}

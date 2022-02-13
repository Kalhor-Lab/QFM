# exp_params
# big_graph_list = all_graphs
# tr_col = "tr"
# gr_col = "gr"
# gr_eval_col = "gr3_eval"
generate_evaluate_tb <- function(exp_params, big_graph_list, tr_col, gr_col, gr_eval_col, perturb = F) {
        map(1:nrow(exp_params), function(j) {
                message(j)
                if (is.null(exp_params[[tr_col]][[j]])) {
                        return(list(recon = NULL,
                                    true = NULL))
                }
                big_graph_id = exp_params$big_graph_id[j]
                tr = exp_params[[tr_col]][[j]]
                sc_celltypes = get_type_from_id(tr$tip.label)
                names(sc_celltypes) = tr$tip.label

                if (perturb) {
                        sc_celltypes = perturb_celltypes(sc_celltypes)
                }

                gr_eval = exp_params[[gr_eval_col]][[j]]

                gr = name_nodes(exp_params[[gr_col]][[j]])
                gr_tr_data = make_gr_tr_data(gr, tr, sc_celltypes)

                # truth
                gr_data_true = make_true_gr_data(big_graph_list[[big_graph_id]])

                # collected sample sizes
                if (exp_params$sampling[j] == "fixed") {
                        tip_collect_size = rep(exp_params$data[[j]]$sample_size, length(gr_data_true$gr$tip.label))
                        names(tip_collect_size) = gr_data_true$gr$tip.label
                } else {
                        tip_collect_size = exp_params$data[[j]]$sample_size
                }
                gr_dd_collect_size = map(gr_data_true$gr_dd, function(x) {
                        s1 = ifelse(x[1] %in% gr_data_true$gr$tip.label, tip_collect_size[x[1]],
                                    sum(tip_collect_size[gr_data_true$gr_tip_list[[x[1]]]]))
                        s2 = ifelse(x[2] %in% gr_data_true$gr$tip.label, tip_collect_size[x[2]],
                                    sum(tip_collect_size[gr_data_true$gr_tip_list[[x[2]]]]))
                        out = c(s1, s2)
                        names(out) = x
                        out
                })
                # true split ratio
                node_split_ratio = map(big_graph_list[[big_graph_id]]$node_id, function(node) {
                        out = big_graph_list[[big_graph_id]]$diff_mode_probs[[node]][1:2]
                        names(out) = big_graph_list[[big_graph_id]]$merge[node, ]
                        out
                })
                names(node_split_ratio) = big_graph_list[[big_graph_id]]$node_id
                node_split_sampled = exp_params$data[[j]]$sampled_field_split
                names(node_split_sampled) = big_graph_list[[big_graph_id]]$node_id
                # end truth

                # mapped truth for each node in reconstructed
                gr_node_mapped_all = get_node_mapping(gr_tr_data, big_graph_list[[big_graph_id]])
                # reconstructed node for each truth (reverse)
                gr_node_true = names(gr_data_true$gr_tip_list[1:Nnode(gr)])
                gr_node_true_mapped = names(gr_node_mapped_all)[match(gr_node_true, gr_node_mapped_all)]
                names(gr_node_true_mapped) = gr_node_true

                tr_node_assign = assign_node_states(gr_tr_data)
                gr_trans_time = est_transition_time(gr_tr_data, tr_node_assign, stat_func = mean)

                gr_tr_data$tr_edges_tb = mutate(gr_tr_data$tr_edges_tb,
                                                final = (max(from_time) - from_time - length) <= 1e-10,
                                                final_num = as.numeric(final))

                gr_tr_data = update_edge_tb_state(gr_tr_data, tr_node_assign)
                gr_node_size = get_node_size(gr_tr_data, tr_node_assign, gr_trans_time)

                gr_node_size_mapped = map(gr_node_size, function(x) {
                        names(x) = c(gr_node_mapped_all, gr_tr_data$gr_tips)[names(x)]
                        x
                })

                node_size_sampled_true = exp_params$data[[j]]$sampled_field_size

                gr_node_all = names(gr_node_mapped_all)

                node_collect_size = map(gr_data_true$gr_tip_list, function(x) {
                        purrr::reduce(tip_collect_size[x], `+`)
                })
                names(node_collect_size) = names(gr_data_true$gr_tip_list)

                # new edge_len based robust stat
                node_mean_final = map_dbl(gr_tr_data$gr$node.label, function(x) {
                        edge_diff = get_edge_diff(x, gr_tr_data, gr_trans_time)
                        mean(edge_diff$final_num)
                })
                names(node_mean_final) = gr_tr_data$gr$node.label

                # classifying probabilities of undersampled edge
                # gr_tr_data$tr_edges_state_tb = classify_edges(gr_tr_data$tr_edges_state_tb)
                # frac_short_edge = map_dbl(gr_node_all, function(x) {
                #         mean(filter(gr_tr_data$tr_edges_state_tb,
                #                     from_type == x)$class == "low")
                # })
                # names(frac_short_edge) = gr_node_all

                out_tb_recon = tibble(j = j,
                                      # gr_eval
                                      kc0 = gr_eval$kc0,
                                      kc1 = gr_eval$kc1,
                                      rf = gr_eval$rf,
                                      # reconstructed
                                      node_gr = gr_node_all,
                                      gr_time_est = gr_tr_data$gr_node_time[gr_node_all] + gr_tr_data$gr_root_len,
                                      gr_time_trans = gr_trans_time[gr_node_all],
                                      gr_node_size = gr_node_size[gr_node_all],
                                      gr_node_size_mapped = gr_node_size_mapped[gr_node_all],
                                      gr_node_mean_final = node_mean_final[gr_node_all],
                                      # gr_node_short_edge = frac_short_edge[gr_node_all],
                                      # mapped truth
                                      node = gr_node_mapped_all,
                                      node_time = gr_data_true$gr_node_time[gr_node_mapped_all] + big_graph_list[[big_graph_id]]$root_time,
                                      node_size = get_true_size(big_graph_list[[big_graph_id]])[gr_node_mapped_all],
                                      node_size_sampled = node_size_sampled_true[gr_node_mapped_all],
                                      node_split = node_split_ratio[gr_node_mapped_all],
                                      node_split_sampled = node_split_sampled[gr_node_mapped_all],
                                      node_dd_collect = gr_dd_collect_size[gr_node_mapped_all],
                                      node_size_collect = node_collect_size[gr_node_mapped_all]
                )
                # comparison based on truth
                out_tb_true = tibble(j = j,
                                     # gr_eval
                                     kc0 = gr_eval$kc0,
                                     kc1 = gr_eval$kc1,
                                     rf = gr_eval$rf,
                                     # truth
                                     node = gr_node_true,
                                     node_time = gr_data_true$gr_node_time[gr_node_true] + big_graph_list[[big_graph_id]]$root_time,
                                     node_size = get_true_size(big_graph_list[[big_graph_id]])[gr_node_true],
                                     node_size_sampled = node_size_sampled_true[gr_node_true],
                                     node_split = node_split_ratio[gr_node_true],
                                     node_split_sampled = node_split_sampled[gr_node_true],
                                     node_dd_collect = gr_dd_collect_size[gr_node_true],
                                     node_size_collect = node_collect_size[gr_node_true],
                                     # mapped reconstructed
                                     node_gr = gr_node_true_mapped,
                                     gr_time_est = gr_tr_data$gr_node_time[gr_node_true_mapped] + gr_tr_data$gr_root_len,
                                     gr_time_trans = gr_trans_time[gr_node_true_mapped],
                                     gr_node_size = gr_node_size[gr_node_true_mapped],
                                     gr_node_size_mapped = gr_node_size_mapped[gr_node_true_mapped],
                                     gr_node_mean_final = node_mean_final[gr_node_true_mapped],
                                     # gr_node_short_edge = frac_short_edge[gr_node_true_mapped]
                )
                list(recon = out_tb_recon,
                     true = out_tb_true)
        })
}

perturb_celltypes <- function(sc_celltypes, pct = 0.05, min_keep = 10) {
        sc_celltype_counts = table(sc_celltypes)
        perturb_size = pmin(ceiling(sc_celltype_counts * pct), pmax(sc_celltype_counts - min_keep, 0))
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
}

evaluate_node_assign <- function(exp_params, big_graph_list, tr_col, gr_col, perturb = F) {
        map(1:nrow(exp_params), function(j) {
                message(j)
                if (is.null(exp_params[[tr_col]][[j]])) {
                        return(NULL)
                }
                big_graph_id = exp_params$big_graph_id[j]
                tr = exp_params[[tr_col]][[j]]
                sc_celltypes = get_type_from_id(tr$tip.label)
                names(sc_celltypes) = tr$tip.label

                if (perturb) {
                        sc_celltypes = perturb_celltypes(sc_celltypes)
                }
                gr = name_nodes(exp_params[[gr_col]][[j]])
                gr_tr_data = make_gr_tr_data(gr, tr, sc_celltypes)
                tr_node_assign = assign_node_states(gr_tr_data)
                gr_node_mapped_all = get_node_mapping(gr_tr_data, big_graph_list[[big_graph_id]])

                node_size_sampled_true = exp_params$data[[j]]$sampled_field_size
                node_size_true = get_true_size(big_graph_list[[big_graph_id]])

                out_tb_raw = tibble(type_true = get_type_from_id(names(tr_node_assign)),
                                type_assign = gr_node_mapped_all[tr_node_assign]) %>%
                        mutate(correct = (type_true == type_assign)) %>%
                        filter(!is.na(type_assign))
                # out_tb = out_tb_raw %>%
                #         group_by(type_assign) %>%
                #         summarize(frac_correct = mean(correct)) %>%
                #         arrange(type_assign) %>%
                #         mutate(j = j)
                # out_tb$sampling_frac = node_size_sampled_true[out_tb$type_assign] / node_size_true[out_tb$type_assign]
                mutate(out_tb_raw, j = j)
        })
}

get_true_size <- function(type_graph) {
        gens0 = make_gens(type_graph)
        out = map_dbl(gens0[c(type_graph$node_id, type_graph$tip_id)], function(x) {
                ifelse(x$active,
                       x$end_count/2,
                       x$end_count)
        })
        names(out) = c(type_graph$node_id, type_graph$tip_id)
        out
}

make_true_gr_data <- function(type_graph) {
        gr_tips = type_graph$tip_id
        names(gr_tips) = gr_tips
        gr_true_phy = as.phylo(type_graph)
        gr_true_phy$node.label = str_replace_all(gr_true_phy$node.label, ":", "")
        gr_tip_list_true = list_dd_and_tips(gr_true_phy)$tips
        gr_tip_list_true = c(gr_tip_list_true, gr_tips)
        gr_dd_true = list_dd_and_tips(gr_true_phy)$dd

        gr_node_time_truth = node.depth.edgelength(gr_true_phy)
        names(gr_node_time_truth) = c(gr_true_phy$tip.label, gr_true_phy$node.label)

        list(gr = gr_true_phy,
             gr_node_time = gr_node_time_truth,
             gr_tip_list = gr_tip_list_true,
             gr_dd = gr_dd_true)
}
get_node_mapping <- function(data_obj, type_graph) {
        gr = data_obj$gr
        true_data_obj = make_true_gr_data(type_graph)

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

# library(mclust)
# classify_edges <- function(edge_tb) {
#         mcl_out = Mclust(edge_tb$length, G = 2)
#         mcl_mean = mcl_out$parameters$mean
#         message(paste(c("mixture means:", mcl_mean), collapse = "_"))
#         lo_class = names(mcl_mean)[which.min(mcl_mean)]
#         edge_tb$class = factor(mcl_out$classification == lo_class, labels = c("high", "low"))
#         edge_tb
# }
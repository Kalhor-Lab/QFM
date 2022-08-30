# contains functions for running, processing and evaluating simulated experiments
distribute_sample_size <- function(type_graph, sample_size, sampling = c("fixed", "proportional"), stochastic = F) {
        if (sampling == "fixed") {
                ss = sample_size
        }
        if (sampling == "proportional") {
                gens0 = make_gens_mod2(type_graph)
                tip_size = map_dbl(gens0[type_graph$tip_id], "end_count")
                min_size = ceiling(sample_size / 25)
                if (stochastic) {
                        tip_sample = extraDistr::rmvhyper(nn = 1, k = length(type_graph$tip_id) * sample_size, n = tip_size)[1, ]
                } else {
                        tip_sample = round_frac(length(type_graph$tip_id) * sample_size,  tip_size / sum(tip_size))
                }
                tip_sample = pmax(tip_sample, min_size)
                names(tip_sample) = type_graph$tip_id
                ss = tip_sample
        }
        return(ss)
}
# mut_p can be provided optionally
run_experiment <- function(exp_params, tree_panel, i, rep_dir, mut_p = NULL, delay = F) {
        big_graph_id = exp_params$big_graph_id[i]
        sample_size = exp_params$sample_size[i]
        sampling = exp_params$sampling[i]
        if (is.null(mut_p)) {
                mut_p = exp_params$mut_p[[i]]
        }
        ss = distribute_sample_size(tree_panel$type_graph[[big_graph_id]], sample_size, sampling)
        out_file = paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        if (!file.exists(out_file)) {
                try({
                        out = simulate_sc_data_mod2(tree_panel$type_graph[[big_graph_id]], mut_p, sample_size = ss, delay = delay)
                        saveRDS(out, file = out_file)
                })
        }
}
load_exp_data <- function(exp_params, rep_dir) {
        data_list = map(1:nrow(exp_params), function(i) {
                out_file = paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
                try({
                        out = readRDS(out_file)
                        out
                })
        })
        exp_params$sc = map(data_list, "sc")
        exp_params$tr = map(data_list, "tr")
        exp_params$true_sampled_sizes = map(data_list, "true_sampled_sizes")
        exp_params
}
load_true_sampled_sizes <- function(exp_params, data_dir) {
        exp_params$true_sampled_sizes = map(1:nrow(exp_params), function(i) {
                out_file = paste0(data_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
                try({
                        out = readRDS(out_file)
                        return(out$true_sampled_sizes)
                })
        })
        exp_params
}
dump_tr <- function(exp_params, tr_col, out_dir) {
        if (!dir.exists(out_dir)) {
                dir.create(out_dir)
        }
        walk(1:nrow(exp_params), function(i) {
                tr = exp_params[[tr_col]][[i]]
                out_file = paste0(out_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
                saveRDS(tr, out_file)
        })
}
run_ice_fase <- function(exp_params, tree_panel, tr_col = "tr") {
        exp_params$res = future_map2(exp_params[[tr_col]], exp_params$big_graph_id, function(tr, graph_id) {
                tt = tree_panel$type_graph[[graph_id]]$target_time
                res = ice_fase_mod(tr,
                                   get_type_from_id(tr$tip.label),
                                   total_time = tt - 0.6,
                                   root_time = 0.6,
                                   theta = 0.)
                res
        }, .progress = T, .options = furrr_options(seed = T))
        exp_params$res = map2(exp_params$big_graph_id, exp_params$res, function(i, res) {
                res = evaluate_gr_res(res, tree_panel$gr[[i]])
                res
        })
        exp_params$kc0 = map_dbl(exp_params$res, function(x) {
                x$gr_eval$kc0
        })
        exp_params$bsum = tree_panel$bsum[exp_params$big_graph_id]
        exp_params$num_tip = tree_panel$num_tip[exp_params$big_graph_id]
        exp_params
}
walk_ice_fase <- function(exp_params, tree_panel, tr_col, out_dir) {
        tr_list = exp_params[[tr_col]]
        graph_id_vec = exp_params$big_graph_id
        future_walk(1:nrow(exp_params), function(j) {
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                if (!file.exists(out_file)) {
                        tr = tr_list[[j]]
                        graph_id = graph_id_vec[j]
                        tt = tree_panel$type_graph[[graph_id]]$target_time
                        res = ice_fase_mod(tr,
                                           get_type_from_id(tr$tip.label),
                                           total_time = tt - 0.6,
                                           root_time = 0.6,
                                           theta = 0.)
                        saveRDS(res, out_file)
                }
        }, .progress = T, .options = furrr_options(seed = T))
}
walk_ice_fase_data <- function(exp_params, input_dir, tree_panel, out_dir) {
        future_walk(1:nrow(exp_params), function(j) {
                data_out = readRDS(paste0(input_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds"))
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                if (!file.exists(out_file)) {
                        tr = data_out$tr
                        graph_id = exp_params$big_graph_id[j]
                        tt = tree_panel$type_graph[[graph_id]]$target_time
                        res = ice_fase_mod(tr,
                                           get_type_from_id(tr$tip.label),
                                           total_time = tt - 0.6,
                                           root_time = 0.6,
                                           theta = 0.)
                        saveRDS(res, out_file)
                }
        }, .progress = T, .options = furrr_options(seed = T))
}
walk_ice_fase_data_gr <- function(exp_params, input_dir, tree_panel, out_dir) {
        future_walk(1:nrow(exp_params), function(j) {
                data_out = readRDS(paste0(input_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds"))
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                if (!file.exists(out_file)) {
                        tr = data_out$tr
                        graph_id = exp_params$big_graph_id[j]
                        tt = tree_panel$type_graph[[graph_id]]$target_time
                        res = ice_fase_mod(tr,
                                           get_type_from_id(tr$tip.label),
                                           total_time = tt - 0.6,
                                           root_time = 0.6,
                                           theta = 0.,
                                           gr = tree_panel$gr[[graph_id]])
                        saveRDS(res, out_file)
                }
        }, .progress = T, .options = furrr_options(seed = T))
}
walk_ice_fase_data_tr3 <- function(exp_params, input_dir, tree_panel, out_dir) {
        future_walk(1:nrow(exp_params), function(j) {
                tr3 = readRDS(paste0(input_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds"))
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                if (!file.exists(out_file)) {
                        tr = tr3
                        graph_id = exp_params$big_graph_id[j]
                        tt = tree_panel$type_graph[[graph_id]]$target_time
                        res = ice_fase_mod(tr,
                                           get_type_from_id(tr$tip.label),
                                           total_time = tt - 0.6,
                                           root_time = 0.6,
                                           theta = 0.)
                        saveRDS(res, out_file)
                }
        }, .progress = T, .options = furrr_options(seed = T))
}
walk_ice_fase_data_tr3_gr <- function(exp_params, input_dir, tree_panel, out_dir) {
        future_walk(1:nrow(exp_params), function(j) {
                tr3 = readRDS(paste0(input_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds"))
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                if (!file.exists(out_file)) {
                        tr = tr3
                        graph_id = exp_params$big_graph_id[j]
                        tt = tree_panel$type_graph[[graph_id]]$target_time
                        res = ice_fase_mod(tr,
                                           get_type_from_id(tr$tip.label),
                                           total_time = tt - 0.6,
                                           root_time = 0.6,
                                           theta = 0.,
                                           gr = tree_panel$gr[[graph_id]])
                        saveRDS(res, out_file)
                }
        }, .progress = T, .options = furrr_options(seed = T))
}

run_phylotime_data <- function(exp_params, j, input_dir, tree_panel, mut_p, out_dir) {
        data_out = readRDS(paste0(input_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds"))
        out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
        if (!file.exists(out_file)) {
                sc = data_out$sc
                graph_id = exp_params$big_graph_id[j]
                tt = tree_panel$type_graph[[graph_id]]$target_time
                tr = phylotime(sc_mat = sc, t_total = tt - 0.6, mut_p = mut_p, parallel = T)
                saveRDS(tr, out_file)
        }
}
walk_phylotime_data <- function(exp_params, input_dir, tree_panel, mut_p, out_dir) {
        walk(1:nrow(exp_params), function(j) {
                message(j)
                run_phylotime_data(exp_params, j, input_dir, tree_panel, mut_p, out_dir)
        })
}
walk_ice_fase_gr <- function(exp_params, tree_panel, tr_col, out_dir) {
        tr_list = exp_params[[tr_col]]
        graph_id_vec = exp_params$big_graph_id
        future_walk(1:nrow(exp_params), function(j) {
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                if (!file.exists(out_file)) {
                        tr = tr_list[[j]]
                        graph_id = graph_id_vec[j]
                        tt = tree_panel$type_graph[[graph_id]]$target_time
                        res = ice_fase_mod_gr(tr,
                                              get_type_from_id(tr$tip.label),
                                              gr = tree_panel$gr[[graph_id]],
                                              total_time = tt - 0.6,
                                              root_time = 0.6,
                                              theta = 0.)
                        saveRDS(res, out_file)
                }
        }, .progress = T, .options = furrr_options(seed = T))
}
gather_ice_fase <- function(exp_params, tree_panel, out_dir) {
        # exp_params$bsum = tree_panel$bsum[exp_params$big_graph_id]
        # exp_params$num_tip = tree_panel$num_tip[exp_params$big_graph_id]
        exp_params$res = map(1:nrow(exp_params), function(j) {
                message(j)
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                readRDS(out_file)
        })
        exp_params$res = map2(exp_params$big_graph_id, exp_params$res, function(i, res) {
                res = evaluate_gr_res(res, tree_panel$gr[[i]])
                res
        })
        exp_params$kc0 = map_dbl(exp_params$res, function(x) {
                x$gr_eval$kc0
        })
        exp_params
}
run_qfm_eval <- function(exp_params, tree_panel, res_col = "res", out_col = "qfm_eval", delay = F) {
        exp_params[[out_col]] = map(1:nrow(exp_params), function(i) {
                message(i)
                eval_out = evaluate_qfm_v1(exp_params[[res_col]][[i]],
                                           tree_panel$type_graph[[exp_params$big_graph_id[i]]],
                                           exp_params$true_sampled_sizes[[i]],
                                           delay = delay)
                eval_out
        })
        exp_params
}
extract_evals <- function(exp_params, eval_col = "qfm_eval", log2_sampling_cutoff = -2.0) {
        eval_tb = bind_rows(map(1:nrow(exp_params), function(i) {
                mutate(exp_params[[eval_col]][[i]]$recon,
                       exp_id = i,
                       sampling = exp_params$sampling[i],
                       num_tip = exp_params$num_tip[i],
                       bsum = exp_params$bsum[i],
                       target_time = exp_params$target_time[i]
                       )
        }))
        eval_tb$suff_sampled = eval_tb$log2_node_sampled > log2_sampling_cutoff
        eval_tb = eval_tb %>%
                mutate(log2_node_size = log2(node_size),
                       log2_gr_node_size_in = log2(gr_node_size_in))
        eval_tb$gr_time_trans_abs_error = abs(eval_tb$gr_time_trans_error)
        eval_tb$gr_node_size_abs_logfc = abs(eval_tb$gr_node_size_logfc)
        eval_tb$gr_node_split_abs_error = abs(eval_tb$gr_node_split_order - eval_tb$node_split_order)
        eval_tb
}
extract_evals_true <- function(exp_params, eval_col = "qfm_eval", log2_sampling_cutoff = -2.0) {
        eval_tb = bind_rows(map(1:nrow(exp_params), function(i) {
                mutate(exp_params[[eval_col]][[i]]$true,
                       exp_id = i,
                       sampling = exp_params$sampling[i],
                       num_tip = exp_params$num_tip[i],
                       bsum = exp_params$bsum[i],
                       target_time = exp_params$target_time[i]
                )
        }))
        eval_tb$suff_sampled = eval_tb$log2_node_sampled > log2_sampling_cutoff
        eval_tb = eval_tb %>%
                mutate(log2_node_size = log2(node_size),
                       log2_gr_node_size_in = log2(gr_node_size_in))
        eval_tb$gr_time_trans_abs_error = abs(eval_tb$gr_time_trans_error)
        eval_tb$gr_node_size_abs_logfc = abs(eval_tb$gr_node_size_logfc)
        eval_tb$gr_node_split_abs_error = abs(eval_tb$gr_node_split_order - eval_tb$node_split_order)
        eval_tb
}
plot_evals <- function(exp_params, facet.by) {
        eval_tb = extract_evals(exp_params)
        sampling_col = RColorBrewer::brewer.pal(9, "Spectral")
        g1 = (eval_tb %>%
                      # filter(suff_sampled) %>%
                      ggscatter(x = "node_time", y=  "gr_time_trans",
                                color = "log2_node_sampled",
                                size = 0.15,
                                alpha = 0.05,
                                # facet.by = c("sampling", "num_tip"),
                                # xlab = "True commitment time",
                                # ylab = "Inferred commitment time",
                                xlab = "", ylab = ""
                      ) %>% facet(facet.by = facet.by, nrow = 2) +
                      geom_abline(size = 0.25)  +
                      geom_smooth(aes(x =  node_time, y = gr_time_trans),
                                  color = sampling_col[9],
                                  method = "loess",
                                  se = F,
                                  span = 1,
                                  data = filter(eval_tb, suff_sampled == "TRUE")) +
                      geom_smooth(aes(x =  node_time, y = gr_time_trans),
                                  color = sampling_col[1],
                                  method = "loess",
                                  span = 1,
                                  se = F,
                                  data = filter(eval_tb, suff_sampled == "FALSE")) +
                      stat_cor(aes(group = suff_sampled,
                                   label =  ..r.label..),
                               method = "spearman",
                               cor.coef.name = "rho",
                               digits = 3) +
                      scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                                             breaks = (-4):(-1), na.value = sampling_col[1]) +
                      theme(text = element_text(size = 10),
                            # legend.position = "top",
                            legend.position = "none",
                            strip.background = element_blank()))

        # saveRDS(exp_params, file = "./intermediate_data/panel_mod2/exp_data_1rep_proc.rds")
        # eval_tb$log2_node_size = eval_tb$log2_node_size

        g2 = (eval_tb %>%
                      ggscatter(x = "log2_node_size",
                                y = "log2_gr_node_size_in",
                                size = 0.15,
                                alpha = 0.05,
                                # facet.by = c("sampling", "num_tip"),
                                color = "log2_node_sampled") %>%
                      facet(facet.by = facet.by, nrow = 2) +
                      geom_abline() +
                      # geom_smooth() +
                      geom_smooth(aes(x =  log2_node_size, y = log2_gr_node_size_in),
                                  color = sampling_col[9],
                                  method = "loess",
                                  se = F,
                                  span = 1,
                                  data = filter(eval_tb, suff_sampled == "TRUE")) +
                      geom_smooth(aes(x =  log2_node_size, y = log2_gr_node_size_in),
                                  color = sampling_col[1],
                                  method = "loess",
                                  span = 1,
                                  se = F,
                                  data = filter(eval_tb, suff_sampled == "FALSE")) +
                      xlab("") +
                      ylab("") +
                      # xlab("log2(Progenitor population size)") +
                      # ylab("log2(Inferred progenitor population size)") +
                      stat_cor(aes(group = suff_sampled,
                                   label =  ..r.label..),
                               method = "spearman",
                               cor.coef.name = "rho",
                               digits = 3) +
                      scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                                             breaks = (-4):(-1), na.value = sampling_col[1]) +
                      theme(text = element_text(size = 10), legend.position = "none", strip.background = element_blank()))


        g3 = (eval_tb %>%
                      # filter(suff_sampled) %>%
                      ggscatter(x = "node_split_order", y = "gr_node_split_order",
                                # facet.by = c("sampling", "num_tip"),
                                color = "log2_node_sampled",
                                alpha = 0.05,
                                size = 0.15,
                                # xlab = "True commitment ratio",
                                # ylab = "Inferred commitment ratio"
                                xlab = "",
                                ylab = ""
                      ) %>%
                      facet(facet.by = facet.by, nrow = 2) +
                      geom_abline() +
                      geom_smooth(aes(x =  node_split_order, y = gr_node_split_order),
                                  color = sampling_col[9],
                                  method = "loess",
                                  se = F,
                                  span = 1,
                                  data = filter(eval_tb, suff_sampled == "TRUE")) +
                      geom_smooth(aes(x =  node_split_order, y = gr_node_split_order),
                                  color = sampling_col[1],
                                  method = "loess",
                                  span = 1,
                                  se = F,
                                  data = filter(eval_tb, suff_sampled == "FALSE")) +
                      stat_cor(aes(group = suff_sampled,
                                   label =  ..r.label..),
                               method = "spearman",
                               cor.coef.name = "rho",
                               digits = 3) +
                      scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                                             breaks = (-4):(-1), na.value = sampling_col[1]) +
                      theme(text = element_text(size = 10), legend.position = "none", strip.background = element_blank()))
        (g1 / g2 / g3)
}

plot_eval_v1 <- function(exp_eval_all) {
        g1 = exp_eval_all %>%
                filter(!is.na(suff_sampled)) %>%
                ggboxplot(x = "sample_size",
                          y = "gr_time_trans_abs_error",
                          ylab = "abs(Time error)",
                          fill = "sample_size",
                          facet.by = c("sampling", "suff_sampled"))
        g2 = exp_eval_all %>%
                filter(!is.na(suff_sampled)) %>%
                ggboxplot(x = "sample_size",
                          y = "gr_node_size_abs_logfc",
                          ylab = "abs(Size logfc)",
                          fill = "sample_size",
                          facet.by = c("sampling", "suff_sampled"))
        g3 = exp_eval_all %>%
                filter(!is.na(suff_sampled)) %>%
                ggboxplot(x = "sample_size",
                          y = "gr_node_split_abs_error",
                          ylab = "abs(Bias error)",
                          fill = "sample_size",
                          facet.by = c("sampling", "suff_sampled"))
        g1 + g2 + g3 + plot_layout(guide = "collect") & theme(legend.position = 'top')

        # g1 = exp_eval_all %>% ggscatter(x = "log2_node_sampled", y = "gr_time_trans_abs_error",
        #                                 alpha = 0.1,
        #                                 ylab = "abs(Time error)",
        #                                 color = "sample_size",
        #                                 facet.by = c("sampling")) +
        #         geom_smooth(aes(color = sample_size, group = sample_size), se = F)
        # g2 = exp_eval_all %>% ggscatter(x = "log2_node_sampled", y = "gr_node_size_abs_logfc",
        #                                 alpha = 0.1,
        #                                 ylab = "abs(Size logfc)",
        #                                 color = "sample_size",
        #                                 facet.by = c("sampling")) +
        #         geom_smooth(aes(color = sample_size, group = sample_size), se = F)
        # g3 = exp_eval_all %>% ggscatter(x = "log2_node_sampled", y = "gr_node_split_abs_error",
        #                                 alpha = 0.1,
        #                                 ylab = "abs(Bias error)",
        #                                 color = "sample_size",
        #                                 facet.by = c("sampling")) +
        #         geom_smooth(aes(color = sample_size, group = sample_size), se = F)
        # g1 + g2 + g3 + plot_layout(guide = "collect") & theme(legend.position = 'top')

}
run_fase_eval <- function(exp_params, tree_panel) {
        exp_params$fase_eval = future_map2(exp_params$tr, exp_params$big_graph_id, function(tr, graph_id) {
                summarize_fase_results(tr, tree_panel$type_graph[[graph_id]])
        }, .progress = T)
        exp_params$fase_eval = map2(exp_params$fase_eval, exp_params$data, function(eval_tb, x) {
                sample_size_vec = map_dbl(x$true_sampled_sizes, 1)
                out = left_join(eval_tb,
                                tibble(state = names(sample_size_vec),
                                       log2_sample_size = log2(sample_size_vec)),
                                by = "state")
                out = mutate(out, log2_sample_frac = log2_sample_size - log2_state_size)
                out
        })
        exp_params$mean_n_fase = map_dbl(exp_params$fase_eval, function(x) {
                mean(x$n_fase)
        })
        exp_params$rmse_fase_time = map_dbl(exp_params$fase_eval, function(x) {
                sqrt(mean(x$error^2))
        })
        exp_params
}
produce_strat_time_sum <- function(tb, group_by) {
        tb  %>%
                group_by(.dots = map(group_by, as.symbol)) %>%
                filter(!is.na(gr_time_trans) & !is.na(node_time)) %>%
                summarise(cc = cor(node_time, gr_time_trans, method = "pearson"),
                          cs = cor(node_time, gr_time_trans, method = "spearman"),
                          rmse = sqrt(mean((node_time - gr_time_trans)^2)),
                          total = n())
}
produce_strat_size_sum <- function(tb, group_by) {
        tb %>%
                group_by(.dots = map(group_by, as.symbol)) %>%
                filter(!is.na(log2_gr_node_size_in) & !is.na(log2_node_size)) %>%
                summarise(cc = cor(log2_node_size, log2_gr_node_size_in, method = "pearson"),
                          cs = cor(log2_node_size, log2_gr_node_size_in, method = "spearman"),
                          rmse = sqrt(mean((log2_node_size - log2_gr_node_size_in)^2)),
                          total = n())
}
produce_strat_split_sum <- function(tb, group_by) {
        tb %>%
                group_by(.dots = map(group_by, as.symbol)) %>%
                filter(!is.na(gr_node_split_order) & !is.na(node_split_order)) %>%
                summarise(cc = cor(node_split_order, gr_node_split_order, method = "pearson"),
                          cs = cor(node_split_order, gr_node_split_order, method = "spearman"),
                          rmse = sqrt(mean((node_split_order - gr_node_split_order)^2)),
                          total = n())
}
generate_panel_gr <- function(tree_panel) {
        tree_panel$gr = map(tree_panel$type_graph, as.phylo_mod2.type_graph)
        tree_panel
}
summarize_all_rmse <- function(eval_tb, group_by) {
        list(produce_strat_time_sum(eval_tb, group_by) %>%
                     summarise(mean_rmse = mean(rmse),
                               mode = "time"),
             produce_strat_size_sum(eval_tb, group_by) %>%
                     summarise(mean_rmse = mean(rmse),
                               mode = "size"),
             produce_strat_split_sum(eval_tb, group_by) %>%
                     summarise(mean_rmse = mean(rmse),
                               mode = "split")) %>%
                bind_rows()
}



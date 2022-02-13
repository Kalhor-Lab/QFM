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

all_graphs = readRDS("../LTModelData/all_graphs_v0.rds")
pec_graphs = readRDS("../LTModelData/pec_graphs_v1.rds")
pec_graphs500 = readRDS("../LTModelData/pec_graphs500_v1.rds")
all_graphs = c(all_graphs, pec_graphs, pec_graphs500)

g_a = plot_type_graph(all_graphs[[2]],
                node_col_mapper = function(x) gr_color(all_graphs[[2]])[x],
                show_node_text = F) + ylab("") +
        theme(text = element_text(size = 6s))
g_b = plot_type_graph(all_graphs[[51]],
                      node_col_mapper = function(x) gr_color(all_graphs[[52]])[x],
                      show_node_text = F) + ylab("") +
        theme(text = element_text(size = 6))
g_c = plot_type_graph(all_graphs[[61]],
                      node_col_mapper = function(x) gr_color(all_graphs[[62]])[x],
                      show_node_text = F) + ylab("") +
        theme(text = element_text(size = 6))

(g_a + g_b + g_c) %>%
        push_pdf(file_name = "graph_type_example",
                 w = 8.5, h = 3.,
                 dir = "../LTModelPlots/panel/")


e1 = readRDS("../LTModelData/all_graphs_v1_exp_data_proc.rds")
e2 = readRDS("../LTModelData/all_graphs_v1prop_exp_data_proc.rds")
e1$sampling = "fixed"
e2$sampling = "proportional"
e_all = bind_rows(e1, e2)

big_graph_list = all_graphs
total_time = 15.0 - 0.6
library(furrr)
plan(multisession, workers = 12)
# runing alternative methods
exp_params = mutate(exp_params, gr = future_map(tr, function(tr) {
        sc_celltypes = get_type_from_id(tr$tip.label)
        names(sc_celltypes) = tr$tip.label
        gr = reconstruct_graph(tr, sc_celltypes, total_time, stat_func = mean, theta = 1.0)
        gr
}, .progress = T))
exp_params = mutate(exp_params,
                    gr_eval = map2(gr, big_graph_id, function(x, i) evalute_gr(x, as.phylo(all_graphs[[i]]))))

exp_params = mutate(exp_params, gr3alt0 = future_map(tr3, function(tr) {
        sc_celltypes = get_type_from_id(tr$tip.label)
        names(sc_celltypes) = tr$tip.label
        gr = reconstruct_graph(tr, sc_celltypes, total_time, stat_func = mean, theta = 1.0)
        gr
}, .progress = T))
exp_params = mutate(exp_params,
                    gr3alt0_eval = map2(gr3alt0, big_graph_id, function(x, i) evalute_gr(x, as.phylo(all_graphs[[i]]))))

exp_params = mutate(exp_params, gr3alt1 = future_map(tr3, function(tr) {
        sc_celltypes = get_type_from_id(tr$tip.label)
        names(sc_celltypes) = tr$tip.label
        gr = reconstruct_graph(tr, sc_celltypes, total_time, stat_func = mean, theta = 0.5)
        gr
}, .progress = T))
exp_params = mutate(exp_params,
                    gr3alt1_eval = map2(gr3alt1, big_graph_id, function(x, i) evalute_gr(x, as.phylo(all_graphs[[i]]))))




# metric = "kc0"
# e_all %>% transmute(big_graph_id, num_tip, graph_type, sampling,
#                     gr3_val = map_dbl(gr3_eval, metric),
#                     # gr3a_val = map_dbl(gr3a_eval, metric),
#                     sps3_gr_val = map_dbl(sps3_gr_eval, metric)) %>%
#         gather(key = "method", value = "val", -c(big_graph_id, num_tip, graph_type, sampling)) %>%
#         ggboxplot(x = "num_tip", y = "val", color = "method", add = "jitter", numeric.x.axis = T) %>%
#         facet(facet.by = c("graph_type", "sampling"), nrow = 1)

# Here, ignoring the alternatives and get
eval_tb_gr3 = generate_evaluate_tb(exp_params, all_graphs, tr_col = "tr3", gr_col = "gr3", gr_eval_col = "gr3_eval")
eval_tb_gr3alt0 = generate_evaluate_tb(exp_params, all_graphs, tr_col = "tr3", gr_col = "gr3alt0", gr_eval_col = "gr3alt0_eval")
eval_tb_gr3alt1 = generate_evaluate_tb(exp_params, all_graphs, tr_col = "tr3", gr_col = "gr3alt1", gr_eval_col = "gr3alt1_eval")
save(eval_tb_gr3, eval_tb_gr3alt0, eval_tb_gr3alt1, file = "../LTModelData/temp_all_graphs_eval.rda")

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
                y = y[2:3]
                return(y[names(x)[1]]/sum(y))
        }))
        out_tb
}
process_node_size <- function(eval_tb) {
        out_tb = mutate(eval_tb,
                        gr_node_size_in = map_dbl(gr_node_size, function(x) {
                                if (is.null(x)) {
                                        return(NA)
                                }
                                return(x[1])
                        }))
        out_tb
}

process_res <- function(eval_tb_list, return_true = F) {
        if (return_true) {
                tb_true_all = bind_rows(map(eval_tb_list, "true"))
                tb_true_all = mutate(tb_true_all, is_resolved = as.numeric(!is.na(node_gr)))
                tb_true_all = process_node_size(tb_true_all)
                tb_true_all = process_split_ratio(tb_true_all)
                tb_true_all = mutate(tb_true_all, log2_node_size = log2(node_size))
                tb_true_all = mutate(tb_true_all, gr_time_trans_error = gr_time_trans - node_time)
                tb_true_all = mutate(tb_true_all, gr_time_est_error = gr_time_est - node_time)
                tb_true_all = mutate(tb_true_all, gr_node_size_logfc = log2(gr_node_size_in / node_size))
                return(tb_true_all)
        } else {
                tb_recon_all = bind_rows(map(eval_tb_list, "recon"))
                tb_recon_all = mutate(tb_recon_all, is_resolved = !is.na(node))
                tb_recon_all = process_node_size(tb_recon_all)
                tb_recon_all = process_split_ratio(tb_recon_all)
                tb_recon_all = mutate(tb_recon_all, log2_node_size = log2(node_size))
                tb_recon_all = mutate(tb_recon_all, gr_time_trans_error = gr_time_trans - node_time)
                tb_recon_all = mutate(tb_recon_all, gr_time_est_error = gr_time_est - node_time)
                tb_recon_all = mutate(tb_recon_all, gr_node_size_logfc = log2(gr_node_size_in / node_size))
                return(tb_recon_all)
        }
}

load("../LTModelData/temp_all_graphs_eval.rda")
eval_tb = eval_tb_gr3alt0
eval_tb_all = process_res(eval_tb, return_true = T)

# Calculating COSAR Error for each node
# library(furrr)
# plan(multisession, workers = 12)
# node_stat_tb = eval_tb_all %>%
#         select(j, node,
#                node_size, node_split, node_size_collect) %>% mutate(
#                         mean_cosar_error = future_pmap_dbl(., function(node_size_collect, node_split, node_size, ...) {
#                                 if (is.na(node_size)) {
#                                         return(NA)
#                                 }
#                                 if (node_size > 2048 & (abs(node_split[1] - 0.5) <= 0.7)) {
#                                         return(-Inf)
#                                 }
#                                 assertthat::assert_that(all(names(node_size_collect) == names(node_split)))
#                                 mean(
#                                         compute_mean_cosar_error(size = node_size,
#                                                                  split_ratio = node_split,
#                                                                  collect_size = node_size_collect)
#                                 )
#         }, .progress = T)
# )
# saveRDS(node_stat_tb, file = "../LTModelData/temp_all_graphs_node_stat.rds")

# add experiement level information TODO: this needs upadting
exp_params$j = 1:nrow(exp_params)
exp_params = mutate(exp_params, num_tip = map_dbl(big_graph_id, function(bid) {
        length(all_graphs[[bid]]$tip_id)
}))

eval_tb_all = left_join(eval_tb_all, select(exp_params, j, big_graph_id, num_tip, sampling, graph_type), by = "j")
eval_tb_all = left_join(eval_tb_all, select(node_stat_tb, j, node, mean_cosar_error), by = c("j", "node"))
robust_cutoff =  -0.5
eval_tb_all = mutate(eval_tb_all, robust = mean_cosar_error > robust_cutoff)

eval_tb_all_nest = nest(eval_tb_all, data = -j)
qvals = map_dbl(1:nrow(eval_tb_all_nest), function(j) {
        if(all(!is.na(eval_tb_all_nest$data[[j]]$node_gr))) {
                big_graph_id = eval_tb_all_nest$data[[j]]$big_graph_id[1]
                cor_null_vec = map_dbl(1:n_null, function(i) {
                        node_seq = sample_node_ranking(as.phylo(all_graphs[[big_graph_id]]))
                        node_rank = 1:length(node_seq)
                        names(node_rank) = node_seq
                        cor(eval_tb_all_nest$data[[j]]$node_time,
                            node_rank[eval_tb_all_nest$data[[j]]$node], method = "spearman")
                })
                cor_obs = cor(eval_tb_all_nest$data[[j]]$node_time,
                              eval_tb_all_nest$data[[j]]$gr_time_trans, method = "spearman")
                out = qvalue::empPvals(cor_obs, cor_null_vec)
                return(out)
        } else {
                return(NA)
        }
})

compute_total_ranking <- function(gr) {
        list_out = list_dd_and_tips(name_nodes(gr))
        gr_dd = list_out$dd
        gr_tip_lists = list_out$tips
        prod(map_dbl(gr_dd, function(x) {
                s1 = pmax(1, length(gr_tip_lists[[x[1]]])) - 1
                s2 = pmax(1, length(gr_tip_lists[[x[2]]])) - 1
                factorial(s1 + s2) / factorial(s1) / factorial(s2)
        }))
}
sample_node_ranking <- function(gr) {
        list_out = list_dd_and_tips(name_nodes(gr))
        gr_dd = list_out$dd
        gr_tip_lists = list_out$tips
        # genrate the draw for each node first
        node_shuffle = map(gr_dd, function(x) {
                s1 = pmax(1, length(gr_tip_lists[[x[1]]])) - 1
                s2 = pmax(1, length(gr_tip_lists[[x[2]]])) - 1
                sample(c(rep(x[1], s1), rep(x[2], s2)))
        })
        root_node = names(gr_dd)[1]
        ranking = rep(root_node, gr$Nnode)
        assign_dd <- function(dd) {
                if (dd %in% gr$tip.label) {
                        return()
                } else {
                        dd_indices = which(ranking == dd)
                        if (length(dd_indices) <= 1) {
                                return()
                        } else {
                                ranking[dd_indices[-1]] <<- node_shuffle[[dd]]
                                assign_dd(gr_dd[[dd]][1])
                                assign_dd(gr_dd[[dd]][2])
                        }
                }
        }
        assign_dd(root_node)
        ranking
}


all_graph_rankings = map_dbl(all_graphs, function(x) compute_total_ranking(as.phylo(x)))
log2(all_graph_rankings[52])
log2(all_graph_rankings[6:10])



# set this tibble to be plotted
eval_tb_plot = eval_tb_all

eval_tb_plot %>%
        filter(!is.na(mean_cosar_error)) %>%
        filter(robust) %>%
        ggscatter(x = "node_split_order", y = "gr_node_split_order",
                  facet.by = c("sampling", "graph_type"),
                  color = "mean_cosar_error",
                  xlab = "True split ratio", ylab = "Estimated split ratio") + geom_abline() + geom_smooth(method = "lm")+
        scale_colour_gradientn(limits = c(-4, 0),
                               colors = c("red", "white", "blue"),
                               breaks = c(-4, -0.5, 0),
                               na.value = "red")

eval_tb_plot %>%
        filter(!is.na(mean_cosar_error)) %>%
        ggscatter(x = "node_split_sampled_order", y = "gr_node_split_order",
                  facet.by = c("robust", "graph_type"),
                  color = "mean_cosar_error",
                  xlab = "Sampled split ratio", ylab = "Estimated split ratio"
        ) + geom_abline() + geom_smooth()+
        scale_colour_gradientn(limits = c(-4, 0),
                               colors = c("red", "white", "blue"),
                               breaks = c(-4, -0.5, 0),
                               na.value = "red")

eval_tb_plot %>%
        filter(!is.na(mean_cosar_error)) %>%
        ggscatter(x = "node_time", y=  "gr_time_trans",
                  color = "mean_cosar_error",
                  facet.by = c("robust", "graph_type"),
                  xlab = "True Time", ylab = "Estimated Time",
        ) +
        geom_abline()  + scale_colour_gradientn(limits = c(-4, 0), colors = c("red", "white", "blue"),
                                                breaks = c(-4, -0.5, 0),na.value = "red")

eval_tb_plot %>%
        filter(!is.na(mean_cosar_error)) %>%
        ggscatter(x = "node_size", y=  "gr_node_size_in",
                  facet.by = c("robust", "graph_type"),
                  color = "mean_cosar_error", xlim = c(0, 4096)) + geom_abline()  +
        xlab("Node Size Total") + ylab("Node Size Estimated") +
        scale_colour_gradientn(limits = c(-4, 0), colors = c("red", "white", "blue"), breaks = c(-4, -0.5, 0),na.value = "red")

eval_tb_plot %>%
        filter(!is.na(mean_cosar_error)) %>%
        ggscatter(x = "node_size_sampled", y=  "gr_node_size_in",
                  facet.by = c("robust", "graph_type"),
                  color = "mean_cosar_error") + geom_abline() +
        xlab("Node Size Sampled") + ylab("Node Size Estimated") + scale_colour_gradientn(limits = c(-4, 0), colors = c("red", "white", "blue"), breaks = c(-4, -0.5, 0),na.value = "red")

# summary plots
(eval_tb_plot %>%
        filter(!is.na(mean_cosar_error)) %>%
        ggboxplot(x = "num_tip",  y = "gr_node_size_logfc", facet.by = c("robust", "graph_type")) +
        geom_hline(yintercept = 0)) %>%
        push_pdf(file_name = "node_size", w = 7.5, h = 2.5, dir = "../LTModelPlots/temp_panel_eval/")

(eval_tb_plot %>%
        filter(!is.na(mean_cosar_error)) %>%
        ggboxplot(x = "num_tip",  y = "gr_time_trans_error", facet.by = c("robust", "graph_type")) +
        geom_hline(yintercept = 0)) %>%
        push_pdf(file_name = "time_error", w = 7.5, h = 2.5, dir = "../LTModelPlots/temp_panel_eval/")

(eval_tb_plot %>% group_by(j) %>%
                summarise(graph_type = graph_type[1],
                          sampling = sampling[1],
                          num_tip = num_tip[1],
                          frac_resolved = mean(is_resolved),
                          frac_robust = mean(robust)) %>%
                ggboxplot(x = "num_tip",  y = "frac_resolved", facet.by = c("sampling", "graph_type")) +
                theme(axis.text.x = element_text(size = 9, angle = 65))
        ) %>%
        push_pdf(file_name = "frac_resolved", w = 7.5, h = 3., dir = "../LTModelPlots/temp_panel_eval/")

(eval_tb_plot %>% group_by(j) %>%
        summarise(graph_type = graph_type[1],
                  sampling = sampling[1],
                  num_tip = num_tip[1],
                  frac_resolved = mean(is_resolved),
                  frac_robust = mean(robust)) %>%
        ggboxplot(x = "num_tip",  y = "frac_robust", facet.by = c("sampling", "graph_type"),
                  ylim = c(0, 1)) + theme(axis.text.x = element_text(size = 9, angle = 65))
        ) %>%
        push_pdf(file_name = "frac_robust", w = 7.5, h = 3., dir = "../LTModelPlots/temp_panel_eval/")
# ggscatter(eval_tb_plot, x = "mean_cosar_error", y = "is_resolved")+ geom_smooth(method = "lm")


# updated ABOVE ----------------------------------------------------------------------------------------------------
# alternative methods below (I think only need to evaluate topology below, which does require running full function)
exp_params = e_all
big_graph_list = all_graphs
total_time = 15.0 - 0.6
# runing alternative methods
exp_params = mutate(exp_params, gr = future_map(tr, function(tr) {
        sc_celltypes = get_type_from_id(tr$tip.label)
        names(sc_celltypes) = tr$tip.label
        gr = reconstruct_graph(tr, sc_celltypes, total_time, stat_func = mean, theta = 1.0)
        gr
}, .progress = T))
exp_params = mutate(exp_params,
                    gr_eval = map2(gr, big_graph_id, function(x, i) evalute_gr(x, as.phylo(all_graphs[[i]]))))

exp_params = mutate(exp_params, gr3alt0 = future_map(tr3, function(tr) {
        sc_celltypes = get_type_from_id(tr$tip.label)
        names(sc_celltypes) = tr$tip.label
        gr = reconstruct_graph(tr, sc_celltypes, total_time, stat_func = mean, theta = 1.0)
        gr
}, .progress = T))
exp_params = mutate(exp_params,
                    gr3alt0_eval = map2(gr3alt0, big_graph_id, function(x, i) evalute_gr(x, as.phylo(all_graphs[[i]]))))

exp_params = mutate(exp_params, gr3alt1 = future_map(tr3, function(tr) {
        sc_celltypes = get_type_from_id(tr$tip.label)
        names(sc_celltypes) = tr$tip.label
        gr = reconstruct_graph(tr, sc_celltypes, total_time, stat_func = mean, theta = 0.5)
        gr
}, .progress = T))
exp_params = mutate(exp_params,
                    gr3alt1_eval = map2(gr3alt1, big_graph_id, function(x, i) evalute_gr(x, as.phylo(all_graphs[[i]]))))

eres0 = generate_exp_evals(exp_params, all_graphs, tr_col = "tr3", gr_col = "gr3", gr_eval_col = "gr3_eval")
eres1 = generate_exp_evals(exp_params, all_graphs, tr_col = "tr3", gr_col = "gr3alt0", gr_eval_col = "gr3alt0_eval")
eres2 = generate_exp_evals(exp_params, all_graphs, tr_col = "tr3", gr_col = "gr3alt1", gr_eval_col = "gr3alt1_eval")
eres3 = generate_exp_evals(exp_params, all_graphs, tr_col = "tr3", gr_col = "sps3_gr", gr_eval_col = "sps3_gr_eval")

eres_t0 = generate_exp_evals(exp_params, all_graphs, tr_col = "tr", gr_col = "gr", gr_eval_col = "gr_eval")

save(eres0, eres1, eres2, eres3, eres_t0, file = "temp_eres_all_updated.rda")

process_res <- function(eres, rev = F) {
        if (rev) {
                time_compare_xrev = bind_rows(map(eres, "gr_compare_rev"))
                time_compare_xrev = mutate(time_compare_xrev, is_resolved = as.numeric(!is.na(node)))
                time_compare_xrev = mutate(time_compare_xrev, log2_node_size_total = log2(node_size_total))
                return(time_compare_xrev)
        }
        time_compare_x = bind_rows(map(eres, "gr_compare"))
        time_compare_x = mutate(time_compare_x, log2_node_size_total = log2(node_size_total))
        time_compare_x = mutate(time_compare_x, trans_time_error = trans_time - true_time)
        time_compare_x = mutate(time_compare_x, est_time_error = est_time - true_time)
        time_compare_x = mutate(time_compare_x, is_resolved = !is.na(node))
        return(time_compare_x)
}

eres_tb0 = process_res(eres0, rev = F)
eres_tb1 = process_res(eres1, rev = F)
eres_tb2 = process_res(eres2, rev = F)
eres_tb3 = process_res(eres3, rev = F)
eres_tbt0 = process_res(eres_t0, rev = F)

summairze_node_resolve <- function(eres_tb) {
        eres_tb %>% group_by(j) %>%
                summarize(graph_type = graph_type[1],
                          sampling = sampling[1],
                          num_tip = num_tip[1],
                          ave_res = mean(is_resolved),
                          ave_trans_error = mean(trans_time_error, na.rm = T))
}

bind_rows(summairze_node_resolve(eres_tb0) %>% mutate(method = "non-disjoint"),
          summairze_node_resolve(eres_tb1) %>% mutate(method = "disijoint"),
          summairze_node_resolve(eres_tb2) %>% mutate(method = "mixed"),
          summairze_node_resolve(eres_tb3) %>% mutate(method = "sps"),
          summairze_node_resolve(eres_tbt0) %>% mutate(method = "true_disjoint")) %>%
        ggboxplot(x = "num_tip", y = "ave_res", color = "method", ylab = "Pct Resolved",
                  numeric.x.axis = T, ylim = c(0, 1.)) %>%
        facet(facet.by = c("graph_type", "sampling"))

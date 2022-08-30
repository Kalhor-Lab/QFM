library(qfm)
library(furrr)
plan(multisession, workers = 12)
# Experiments for errorous readouts
# based on simulated panels
output_dir = "./intermediate_data/panel_mod2_v1/"
plot_dir = "./plots/panel_mod2_v1_examples/"
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))

# extend commitment
# specify the fraction of the cells that commit x generation later
# gens0 = make_gens_mod2(type_graph)
# gens0$`7`$double_time
get_subgraph_nodes <- function(type_graph, node) {
        edges_nest = type_graph$edges %>%
                mutate(in_node = as.character(in_node)) %>%
                mutate(out_node = as.character(out_node)) %>%
                mutate(in_node_dup = in_node) %>% nest(out = -in_node_dup)
        edges_list = edges_nest$out
        names(edges_list) = edges_nest$in_node_dup

        tip_subgraph = character()
        node_subgraph = character()
        cur_nodes = node
        while (length(cur_nodes) > 0) {
                cur_nodes_temp = purrr::reduce(map(cur_nodes, function(x) {
                        out = edges_list[[x]]$out_node
                }), c)
                tip_subgraph = c(tip_subgraph,
                                 cur_nodes_temp[cur_nodes_temp %in% type_graph$tip_id])
                node_subgraph = c(node_subgraph,
                                  cur_nodes_temp[cur_nodes_temp %in% type_graph$node_id])
                cur_nodes = cur_nodes_temp[cur_nodes_temp %in% type_graph$node_id]
        }
        list(nodes = node_subgraph,
             tips = tip_subgraph)
}
ext_state_name <- function(state_vec) {
        map_chr(state_vec, function(state) {
                if (grepl("v", state)) {
                        sp = strsplit(state, "v")
                        state_ext = paste0(sp[[1]][1], "v", as.numeric(sp[[1]][2])+1)
                } else {
                        state_ext = paste0(state, "v1")
                }
                state_ext
        })
}
extend_state_v1 <- function(type_graph, type_mod, delay_gen, ext_prob) {
        gens0 = make_gens_mod2(type_graph)
        max_delay_gen = min(map_dbl(gens0[type_graph$merge[[type_mod]]], "num_gen"))
        assertthat::assert_that(delay_gen <= max_delay_gen)

        subgraph = get_subgraph_nodes(type_graph, type_mod)
        down_nodes = subgraph$nodes
        down_tips = subgraph$tips
        # NOTE: this version only has two downstream states, otherwise needs to copy the entire subgraph
        ext_name = ext_state_name(type_mod)
        ext_node_name = ext_state_name(down_nodes)
        ext_tip_name = ext_state_name(down_tips)

        type_graph$node_id = c(type_graph$node_id, c(ext_name, ext_node_name))
        type_graph$tip_id = c(type_graph$tip_id, ext_tip_name)
        type_graph$all_id = c(type_graph$node_id, type_graph$tip_id)

        # merges
        type_graph$merge[[ext_name]] = ext_state_name(type_graph$merge[[type_mod]])
        type_graph$merge[[type_mod]] = c(type_graph$merge[[type_mod]], ext_name)
        for (i in 1:length(down_nodes)) {
                type_graph$merge[[ext_node_name[i]]] = ext_state_name(type_graph$merge[[down_nodes[i]]])
        }
        # lifetime and difftime
        type_graph$diff_time[[ext_name]] = type_graph$diff_time[[type_mod]] + delay_gen * type_graph$lifetime[[type_mod]]
        type_graph$lifetime[[ext_name]] = type_graph$lifetime[[type_mod]]
        type_graph$prob_loss[[ext_name]] = type_graph$prob_loss[[type_mod]]
        type_graph$prob_pm[[ext_name]] = type_graph$prob_pm[[type_mod]]

        for (i in 1:length(down_nodes)) {
                type_graph$lifetime[[ext_node_name[i]]] = type_graph$lifetime[[down_nodes[i]]]
                type_graph$diff_time[[ext_node_name[i]]] = type_graph$diff_time[[down_nodes[i]]]
                type_graph$prob_loss[[ext_node_name[i]]] = type_graph$prob_loss[[down_nodes[i]]]
                type_graph$prob_pm[[ext_node_name[i]]] = type_graph$prob_pm[[down_nodes[i]]]
        }
        for (i in 1:length(down_tips)) {
                type_graph$lifetime[[ext_tip_name[i]]] = type_graph$lifetime[[down_tips[i]]]
                type_graph$diff_time[[ext_tip_name[i]]] = Inf
                type_graph$prob_loss[[ext_tip_name[i]]] = type_graph$prob_loss[[down_tips[i]]]
                type_graph$prob_pm[[ext_tip_name[i]]] = type_graph$prob_pm[[down_tips[i]]]
        }
        # probs
        type_mod_probs = type_graph$diff_mode_probs[[type_mod]]
        type_graph$diff_mode_probs[[type_mod]] = c(type_mod_probs[1:2] * (1 - ext_prob), ext_prob,
                                                   type_mod_probs[3] * (1 - ext_prob), 0, 0)
        type_graph$diff_mode_probs[[ext_name]] = c(type_mod_probs)
        for (i in 1:length(down_nodes)) {
                type_graph$diff_mode_probs[[ext_node_name[i]]] = type_graph$diff_mode_probs[[down_nodes[i]]]
        }
        type_graph = generate_edge_df_mod2(type_graph)
        type_graph
}

# commitment happens overtime
type_graph = tree_panel$type_graph[[17]]
type_graph_mod1 = extend_state_v1(type_graph, type_mod = "11", delay_gen = 1, ext_prob = 2/3)
type_graph_mod2 = extend_state_v1(type_graph_mod1, type_mod = "11v1", delay_gen = 1, ext_prob = 1/2)

gr_col_vec = gr_color_v1(type_graph)
all_subgraph_nodes = c(subgraph$nodes, subgraph$tips)

gr_col_add = c(gr_col_vec[all_subgraph_nodes],
               gr_col_vec[all_subgraph_nodes])
names(gr_col_add) = c(c(paste0(all_subgraph_nodes, "v1")),
                      c(paste0(all_subgraph_nodes, "v2")))
gr_col_vec = c(gr_col_vec, gr_col_add)

ga = plot_type_graph_clean_mod2(type_graph, node_col_mapper = function(x) gr_col_vec[x], show_node_text = T, show_node_counts = T)
gb = plot_type_graph_clean_mod2(type_graph_mod1, node_col_mapper = function(x) gr_col_vec[x], show_node_text = T, show_node_counts = T)
gc = plot_type_graph_clean_mod2(type_graph_mod2, node_col_mapper = function(x) gr_col_vec[x], show_node_text = T, show_node_counts = T)
(ga + gb + gc) %>% push_pdf("g_var_commit", w = 6.5, h = 3, dir = "./plots/panel_mod2_v1_examples/")

distribute_sample_size_version <- function(type_graph_mod, sample_size, sampling = c("fixed", "proportional")) {
        gens0 = make_gens_mod2(type_graph_mod)
        tip_type_mapped = map_chr(strsplit(type_graph_mod$tip_id, "v"), 1)
        tip_size = map_dbl(gens0[type_graph_mod$tip_id], "end_count")
        tip_size_mapped = split(tip_size, tip_type_mapped)
        tip_size_mapped_total = map_dbl(tip_size_mapped, sum)
        tip_size_mapped_total
        if (sampling == "fixed") {
                tip_sample_size = rep(sample_size, length(tip_size_mapped_total))
                names(tip_sample_size) = names(tip_size_mapped_total)
        } else {
                min_size = ceiling(sample_size / 25)
                # non stochastic sampling
                tip_sample_size = round_frac(length(tip_size_mapped_total) * sample_size,  tip_size_mapped_total / sum(tip_size_mapped_total))
                # tip_sample = extraDistr::rmvhyper(nn = 1, k = length(fm$tip_id) * sample_size, n = tip_size)[1, ]
                tip_sample_size = pmax(tip_sample_size, min_size)
                names(tip_sample_size) = names(tip_size_mapped_total)
        }
        tip_size_split = purrr::reduce(map2(tip_sample_size, tip_size_mapped[names(tip_sample_size)], function(s, m) {
                if (length(m) > 1) {
                        # out = extraDistr::rmvhyper(nn = 1, k = s, n = m)[1, ]
                        out = round_frac(s, m / sum(m))
                        names(out) = names(m)
                } else {
                        out = s
                        names(out) = names(m)
                }
                out
        }), c)
        tip_size_split
}
run_experiment_version <- function(exp_params, tree_panel, i, rep_dir, mut_p = NULL, delay = F) {
        big_graph_id = exp_params$big_graph_id[i]
        sample_size = exp_params$sample_size[i]
        sampling = exp_params$sampling[i]
        if (is.null(mut_p)) {
                mut_p = exp_params$mut_p[[i]]
        }
        ss = distribute_sample_size_version(tree_panel$type_graph[[big_graph_id]], sample_size, sampling)
        out_file = paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        if (!file.exists(out_file)) {
                try({
                        out = simulate_sc_data_mod2(tree_panel$type_graph[[big_graph_id]], mut_p, sample_size = ss, delay = delay)
                        saveRDS(out, file = out_file)
                })
        }
}
map_state_ver <- function(x) {
        map_chr(strsplit(x, "v"), 1)
}

walk_ice_fase_data_ver <- function(exp_params, input_dir, tree_panel, out_dir) {
        future_walk(1:nrow(exp_params), function(j) {
                data_out = readRDS(paste0(input_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds"))
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                if (!file.exists(out_file)) {
                        tr = data_out$tr
                        sc_celltypes_raw = get_type_from_id(tr$tip.label)
                        sc_celltypes_mapped = map_state_ver(sc_celltypes_raw)
                        graph_id = exp_params$big_graph_id[j]
                        tt = tree_panel$type_graph[[graph_id]]$target_time
                        res = ice_fase_mod(tr,
                                           sc_celltypes_mapped,
                                           total_time = tt - 0.6,
                                           root_time = 0.6,
                                           theta = 0.)
                        saveRDS(res, out_file)
                }
        }, .progress = T, .options = furrr_options(seed = T))
}
walk_ice_fase_tr3_ver <- function(exp_params, input_dir, tree_panel, out_dir, true_gr = F) {
        future_walk(1:nrow(exp_params), function(j) {
                tr3 = readRDS(paste0(input_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds"))
                out_file = paste0(out_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                if (!file.exists(out_file)) {
                        tr = tr3
                        sc_celltypes_raw = get_type_from_id(tr$tip.label)
                        sc_celltypes_mapped = map_state_ver(sc_celltypes_raw)
                        graph_id = exp_params$big_graph_id[j]
                        tt = tree_panel$type_graph[[graph_id]]$target_time
                        if (true_gr) {
                                gr_in = tree_panel$gr[[graph_id]]
                        } else {
                                gr_in = NULL
                        }
                        res = ice_fase_mod(tr,
                                           sc_celltypes_mapped,
                                           total_time = tt - 0.6,
                                           root_time = 0.6,
                                           theta = 0.,
                                           gr = gr_in)
                        saveRDS(res, out_file)
                }
        }, .progress = T, .options = furrr_options(seed = T))
}
merge_true_sampled_sizes_vec <- function(true_sampled_sizes, node_mod) {
        all_nodes = names(true_sampled_sizes)
        all_nodes_mapped = map_state_ver(all_nodes)
        all_nodes_unique = sort(unique(all_nodes_mapped))
        true_sampled_sizes_mapped = map(all_nodes_unique, function(x) {
                ss_list = true_sampled_sizes[all_nodes_mapped == x] # for debug
                if (x != node_mod) {
                        purrr::reduce(ss_list, `+`)
                } else {
                        purrr::reduce(map(ss_list, function(size_vec) {
                                out_vec = size_vec[1:3]
                                # if (length(size_vec) == 3) {
                                #         # do nothing
                                # }
                                if (length(size_vec) == 4) {
                                        out_vec[1] = size_vec[1] * (1 - size_vec[4] / sum(size_vec[2:4])) # only the commited fraction of size
                                }
                                out_vec
                        }), `+`)
                }
        })
        names(true_sampled_sizes_mapped) = all_nodes_unique
        true_sampled_sizes_mapped
}
gather_true_sampled_sizes <- function(exp_params, input_dir, node_mod) {
        exp_params$true_sampled_sizes = map(1:nrow(exp_params), function(j) {
                message(j)
                data_out = readRDS(paste0(input_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds"))
                out = merge_true_sampled_sizes_vec(data_out$true_sampled_sizes, node_mod)
                out
        })
        exp_params
}

tree_panel_mod = tibble(type_graph = list(type_graph, type_graph_mod1, type_graph_mod2))
exp_params = tibble(big_graph_id = c(rep(1, 100),
                                     rep(2, 100),
                                     rep(3, 100)),
                    sample_size = 100,
                    sampling = "fixed")
rep_dir = "./intermediate_data/panel_mod2_v1/var_commit/"
# dir.create(rep_dir)
data_dir = paste0(rep_dir, "data/")
# dir.create(data_dir)
# saveRDS(tree_panel_mod, file = paste0(rep_dir, "tree_panel_mod.rds"))
# saveRDS(exp_params, file = paste0(rep_dir, "exp_params.rds"))

tree_panel_mod = readRDS(paste0(rep_dir, "tree_panel_mod.rds"))
exp_params = readRDS(paste0(rep_dir, "exp_params.rds"))

mut_p = readRDS("./intermediate_data/panel_mod2_v1/mut_p_marc1.rds")
# future_walk(1:nrow(exp_params), function(i) {
#         run_experiment_version(exp_params, tree_panel_mod, i, data_dir, mut_p = mut_p)
# })

res_dir = paste0(rep_dir, "res/")
# dir.create(res_dir)
tree_panel_dum = tibble(type_graph = rep(list(tree_panel_mod$type_graph[[1]]), 3))
tree_panel_dum$gr = map(tree_panel_dum$type_graph, as.phylo_mod2.type_graph)

walk_ice_fase_data_ver(exp_params, input_dir = data_dir, tree_panel = tree_panel_mod, out_dir = res_dir)

tr3_dir = paste0(rep_dir, "tr3/")
res_dir = paste0(rep_dir, "res_tr3/")
dir.create(res_dir)
walk_ice_fase_tr3_ver(exp_params, input_dir = tr3_dir, tree_panel = tree_panel_mod, out_dir = res_dir)

res_dir = paste0(rep_dir, "res_gr_tr3/")
dir.create(res_dir)
walk_ice_fase_tr3_ver(exp_params, input_dir = tr3_dir, tree_panel = tree_panel_dum, out_dir = res_dir, true_gr = T)

get_true_diff_time_mod3(tree_panel_mod$type_graph[[1]])[["11"]] + tree_panel_mod$type_graph[[1]]$lifetime[["11"]]
get_true_diff_time_mod3(tree_panel_mod$type_graph[[2]])[["11v1"]] + tree_panel_mod$type_graph[[1]]$lifetime[["11"]]
get_true_diff_time_mod3(tree_panel_mod$type_graph[[3]])[["11v2"]] + tree_panel_mod$type_graph[[1]]$lifetime[["11"]]

exp_params = gather_ice_fase(exp_params, tree_panel = tree_panel_dum, out_dir = res_dir)
exp_params = gather_true_sampled_sizes(exp_params, data_dir, "11")
exp_params = run_qfm_eval(exp_params, tree_panel = tree_panel_dum)

eval_tb = extract_evals(exp_params)
eval_tb$cond = exp_params$big_graph_id[eval_tb$exp_id]
saveRDS(eval_tb, file = paste0(rep_dir, "eval_gr_tb.rds"))

eval_tb = readRDS(paste0(rep_dir, "eval_gr_tb.rds"))
# eval_tb = readRDS(paste0(rep_dir, "eval_tb.rds"))
eval_tb_mod = eval_tb[which(eval_tb$node %in% c("11", "9", "1", "13")), ]
eval_tb_mod = left_join(eval_tb_mod, eval_time_truth)
eval_tb_mod = left_join(eval_tb_mod, eval_size_truth)
eval_tb_mod = left_join(eval_tb_mod, eval_split_truth)

error_summary = eval_tb_mod %>% group_by(cond, node, exp_id) %>%
        summarise(time_rmse = sqrt(mean((true_time - gr_time_trans)^2)),
                  log2_size_rmse = sqrt(mean((log2_true_size - log2_gr_node_size_in)^2)),
                  split_rmse = sqrt(mean((true_split - gr_node_split_order)^2)))

(exp_params %>% ggline(x = "big_graph_id", y = "kc0", add = "mean_se") +
error_summary %>% ggline(x = "cond", y = "time_rmse", col = "node", add = "mean_se") +
error_summary %>% ggline(x = "cond", y = "log2_size_rmse", col = "node", add = "mean_se") +
error_summary %>% ggline(x = "cond", y = "split_rmse", col = "node", add = "mean_se") +
        plot_layout(nrow = 1, guide = "collect") &
        theme(legend.position = "bottom", text =  element_text(size = 10))) %>%
        push_pdf("var_commit_eval_gr_v1", w = 5, h= 2.5, dir = plot_dir)

# g0 = exp_params %>% ggboxplot(x = "big_graph_id", y = "kc0")
eval_time_truth = tibble(cond = c(1:3, 1:3, 1:3, 1:3),
                         node = c(rep("11", 3), rep("9", 3), rep("1", 3), rep("13", 3)),
                         true_time = c(get_true_diff_time_mod3(tree_panel_mod$type_graph[[1]])["11"],
                                       mean(get_true_diff_time_mod3(tree_panel_mod$type_graph[[2]])[c("11", "11v1")]),
                                       mean(get_true_diff_time_mod3(tree_panel_mod$type_graph[[3]])[c("11", "11v1", "11v2")]),
                                       rep(get_true_diff_time_mod3(tree_panel_mod$type_graph[[1]])["9"], 3),
                                       rep(get_true_diff_time_mod3(tree_panel_mod$type_graph[[1]])["1"], 3),
                                       rep(get_true_diff_time_mod3(tree_panel_mod$type_graph[[1]])["13"], 3)))

# g1 = (eval_tb_mod %>% ggboxplot(x = "cond", color = "node", y = "gr_time_trans")) +
#         geom_hline(aes(yintercept = time), data = eval_time_truth)

eval_size_truth = tibble(cond = c(1:3, 1:3, 1:3, 1:3),
                         node = c(rep("11", 3), rep("9", 3), rep("1", 3), rep("13", 3)),
                         log2_true_size = log2(c(get_true_size_mod2(tree_panel_mod$type_graph[[1]])[["11"]][1],
                                  get_true_size_mod2(tree_panel_mod$type_graph[[2]])[["11"]][1] * (1/3) +
                                          get_true_size_mod2(tree_panel_mod$type_graph[[2]])[["11v1"]][1],
                                  get_true_size_mod2(tree_panel_mod$type_graph[[3]])[["11"]][1] * (1/3) +
                                          get_true_size_mod2(tree_panel_mod$type_graph[[3]])[["11v1"]][1] * (1/2) +
                                          get_true_size_mod2(tree_panel_mod$type_graph[[3]])[["11v2"]][1],
                                  rep(get_true_size_mod2(tree_panel_mod$type_graph[[1]])[["9"]][1], 3),
                                  rep(get_true_size_mod2(tree_panel_mod$type_graph[[1]])[["1"]][1], 3),
                                  rep(get_true_size_mod2(tree_panel_mod$type_graph[[1]])[["13"]][1], 3))))
#
# g2 = eval_tb_mod %>% ggboxplot(x = "cond", color = "node", y = "log2_gr_node_size_in") +
#         geom_hline(aes(yintercept = log2_size),
#                    data = eval_size_truth)
eval_split_truth = tibble(node = c("9", "11", "1", "13"),
                          true_split = c(eval_tb_mod$node_split_order[eval_tb_mod$node == "9"][1],
                                         eval_tb_mod$node_split_order[eval_tb_mod$node == "11"][1],
                                         eval_tb_mod$node_split_order[eval_tb_mod$node == "1"][1],
                                         eval_tb_mod$node_split_order[eval_tb_mod$node == "13"][1]))

# g3 = eval_tb_mod %>% ggboxplot(x = "cond", color = "node", y = "gr_node_split_order") +
#         geom_hline(aes(yintercept = split),
#                    data = eval_split_truth)
#
# (g0 + g1 + g2 + g3) +
#         plot_layout(nrow = 1, guide = "collect") &
#         theme(legend.position = "bottom")


# ((exp_params %>% ggboxplot(x = "big_graph_id", y = "kc0", size = 0.4, outlier.size = 0.4)) +
#                 produce_strat_time_sum(eval_tb, c("cond", "exp_id")) %>%  ggboxplot(x = "cond", y = "rmse", size = 0.4, outlier.size = 0.4) +
#                 produce_strat_size_sum(eval_tb, c("cond", "exp_id")) %>% ggboxplot(x = "cond", y = "rmse", size= 0.4, outlier.size = 0.4) +
#                 produce_strat_split_sum(eval_tb, c("cond", "exp_id")) %>% ggboxplot(x = "cond", y = "rmse", size = 0.4, outlier.size = 0.4) +
#                 plot_layout(nrow = 2) & theme(text = element_text(size =10)))















library(qfm)
library(furrr)
plan(multisession, workers = 12)
output_dir = "./intermediate_data/panel_mod2_v1/"

exp_name = "20rep_16type"
exp_params = readRDS(paste0(output_dir, "exp_data_", exp_name, ".rds"))
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))

rep_dir = paste0(output_dir, exp_name, "/")
tr3_dir = paste0(output_dir, exp_name, "/tr3/")

# rep_indices = which(exp_params$big_graph_id == 1 & exp_params$sampling == "fixed" & exp_params$sample_size == 100)
# tr3_list = map(rep_indices, function(i) {
#         readRDS(paste0(tr3_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
# })
# sc_celltypes_list = map(tr3_list, function(tr) get_type_from_id(tr$tip.label))
# res = ice_fase_multi(tr3_list, sc_celltypes_list = sc_celltypes_list, total_time = 12.5 - 0.6, theta = 0.0)
#
#
# data_list = map(rep_indices, function(i) {
#         readRDS(paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
# })
# true_sampled_sizes_all = purrr::reduce(map(data_list, function(x) {
#         x$true_sampled_sizes
# }), merge_node_size)
# true_sampled_sizes_all = map(true_sampled_sizes_all, function(x) x/length(data_list))
# res = evaluate_gr_res(res, tree_panel$gr[[1]])
# evaluate_qfm_v1(res, type_graph = tree_panel$type_graph[[1]], true_sampled_sizes = true_sampled_sizes_all)

multi_tb = as_tibble(expand.grid(big_graph_id = which(tree_panel$num_tip == 16),
                                 sampling = c("fixed", "proportional"),
                                 sample_size = c(50, 100, 200)))
multi_tb$rep_indices = pmap(multi_tb, function(big_graph_id, sampling, sample_size) {
        which(exp_params$big_graph_id == big_graph_id & exp_params$sampling == sampling & exp_params$sample_size == sample_size)
})
# multi_tb$indices_tb = map(multi_tb$rep_indices, function(indices) {
#         tibble(n_rep = 1:4,
#                indices = list(indices[1],
#                               indices[2:3],
#                               indices[4:6],
#                               indices[7:10]))
# })
multi_tb$indices_tb = map(multi_tb$rep_indices, function(indices) {
        # make sure single sample does not repeat
        out1 = tibble(n_rep = 1, i_sim = 1:5)
        out1$indices = map(sample(indices, size = 5, replace = F), function(x) x)
        out2 = tibble(expand.grid(n_rep = 2:5,
                                  i_sim = 1:5))
        out2$indices = map(out2$n_rep, function(x) {
                sample(indices, size = x, replace = F)
                })
        out = bind_rows(list(out1, out2))
        out
})
multi_tb = select(multi_tb, -rep_indices) %>% unnest(cols = "indices_tb")

# res_dir = paste0(rep_dir, "multi_rep_res_gr/")
# dir.create(res_dir)
saveRDS(multi_tb, file = paste0(res_dir, "multi_tb.rds"))

multi_tb = readRDS(paste0(res_dir, "multi_tb.rds"))
future_walk(1:nrow(multi_tb), function(j) {
        out_file = paste0(res_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
        tr3_list = map(multi_tb$indices[[j]], function(i) {
                readRDS(paste0(tr3_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
        })
        sc_celltypes_list = map(tr3_list, function(tr) get_type_from_id(tr$tip.label))
        res = ice_fase_multi(tr3_list,
                             sc_celltypes_list = sc_celltypes_list,
                             total_time = 12.5 - 0.6,
                             theta = 0.0,
                             gr = tree_panel$gr[[multi_tb$big_graph_id[j]]])
        saveRDS(res, file = out_file)
})
multi_tb = gather_ice_fase(multi_tb, tree_panel, res_dir)

# results without topology
res_dir = paste0(rep_dir, "multi_rep_res/")
multi_tb = readRDS(paste0(res_dir, "multi_tb.rds"))
future_walk(1:nrow(multi_tb), function(j) {
        out_file = paste0(res_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
        tr3_list = map(multi_tb$indices[[j]], function(i) {
                readRDS(paste0(tr3_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
        })
        sc_celltypes_list = map(tr3_list, function(tr) get_type_from_id(tr$tip.label))
        res = ice_fase_multi(tr3_list,
                             sc_celltypes_list = sc_celltypes_list,
                             total_time = 12.5 - 0.6,
                             theta = 0.0)
        saveRDS(res, file = out_file)
})
multi_tb = gather_ice_fase(multi_tb, tree_panel, res_dir)
multi_tb$sample_size = factor(multi_tb$sample_size)

multi_tb$true_sampled_sizes = map(multi_tb$indices, function(rep_indices) {
        data_list = map(rep_indices, function(i) {
                readRDS(paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
        })
        true_sampled_sizes_all = purrr::reduce(map(data_list, function(x) {
                x$true_sampled_sizes
        }), merge_node_size)
        true_sampled_sizes_all = map(true_sampled_sizes_all, function(x) x/length(data_list))
        true_sampled_sizes_all
})
multi_tb = run_qfm_eval(multi_tb, tree_panel)
eval_tb = extract_evals(multi_tb)
eval_tb$n_rep = multi_tb$n_rep[eval_tb$exp_id]
eval_tb$sample_size = multi_tb$sample_size[eval_tb$exp_id]
saveRDS(eval_tb, file = paste0(rep_dir, "eval_tb_multi_rep_gr.rds"))

eval_tb = readRDS(paste0(rep_dir, "eval_tb_multi_rep_gr.rds"))
eval_tb$sample_size = factor(eval_tb$sample_size)

multi_tb[c(1, 2, 3, 4, 5, 8)] %>%
        write_csv(file = paste0(rep_dir, "mult_rep_stats.csv"))

multi_tb = read_csv(paste0(rep_dir, "mult_rep_stats.csv"))
multi_tb$sample_size = factor(multi_tb$sample_size)

multi_tb %>% group_by(n_rep, sample_size) %>%
        summarise(mean_kc0 = mean(kc0),
                  median_kc0 = median(kc0)) %>%
        write_csv(paste0(output_dir, "multi_average_kc0.csv"))

summarize_all_rmse(eval_tb, c("n_rep", "sample_size", "exp_id")) %>%
        write_csv(paste0(output_dir, "multi_param_eval.csv"))



g0 = ggline(multi_tb, x = "n_rep", color = "sample_size", y = "kc0", add = "mean_se")

g1 = (produce_strat_time_sum(eval_tb, c("n_rep", "sample_size", "exp_id")) %>%
              ggline(x = "n_rep", y = "rmse", add = "mean_se", color = "sample_size"))

g2 = (produce_strat_size_sum(eval_tb, c("n_rep", "sample_size", "exp_id")) %>%
              ggline(x = "n_rep", y = "rmse", add = "mean_se", color = "sample_size"))

g3 = (produce_strat_split_sum(eval_tb, c("n_rep", "sample_size", "exp_id")) %>%
              ggline(x = "n_rep", y = "rmse", add = "mean_se", color = "sample_size"))

(g0 + g1 + g2 + g3 + plot_layout(guide = "collect", nrow = 1) &
        theme(legend.position = "bottom")) %>%
        push_pdf("multi_eval", w = 6, h = 2.5, dir = "./plots/panel_mod2_v1_examples/")


g1 = (produce_strat_time_sum(eval_tb, c("n_rep", "sample_size", "exp_id", "suff_sampled")) %>%
              ggline(x = "n_rep", y = "rmse", add = "mean_se",
                     color = "sample_size", facet.by= c("suff_sampled")))
g2 = (produce_strat_size_sum(eval_tb, c("n_rep", "sample_size", "exp_id", "suff_sampled")) %>%
              ggline(x = "n_rep", y = "rmse", add = "mean_se", color = "sample_size", facet.by= c("suff_sampled")))

g3 = (produce_strat_split_sum(eval_tb, c("n_rep", "sample_size", "exp_id", "suff_sampled")) %>%
              ggline(x = "n_rep", y = "rmse", add = "mean_se", color = "sample_size", facet.by= c("suff_sampled")))

(g1 + g2 + g3 + plot_layout(guide = "collect") &
                theme(legend.position = "bottom")) %>%
        push_pdf("multi_eval_suff_sample", w = 9., h = 3.5, dir = "./plots/panel_mod2_v1_examples/")





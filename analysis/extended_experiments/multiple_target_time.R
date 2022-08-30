library(qfm)
library(furrr)
plan(multisession, workers = 12)
plot_dir = "./plots/panel_mod2_v1_examples/"

output_dir = "./intermediate_data/panel_mod2_v1/"
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))
index_sel = which(tree_panel$num_tip == 16)
tree_panel_mod = tree_panel[rep(index_sel, 5), ]
tree_panel_mod$target_time = rep(c(11.5, 12.5, 13.5, 14.5, 15.5), each = length(index_sel))
tree_panel_mod$type_graph = map2(tree_panel_mod$type_graph, tree_panel_mod$target_time, function(x, tt) {
        x$target_time = tt
        x
})
exp_params = expand.grid(big_graph_id = 1:nrow(tree_panel_mod),
                         sample_size = 100,
                         sampling = c("fixed", "proportional"),
                         i_sim = 1:2) %>% as_tibble()
exp_params$target_time = map_dbl(exp_params$big_graph_id, function(graph_id) tree_panel_mod$type_graph[[graph_id]]$target_time)

rep_dir = "./intermediate_data/panel_mod2_v1/multi_time_v1/"
dir.create(rep_dir)
data_dir = paste0(rep_dir, "data/")
dir.create(data_dir)
saveRDS(tree_panel_mod, file = paste0(rep_dir, "tree_panel_mod.rds"))
saveRDS(exp_params, file = paste0(rep_dir, "exp_params.rds"))
mut_p = readRDS("./intermediate_data/panel_mod2_v1/mut_p_marc1.rds")
future_walk(1:nrow(exp_params), function(i) {
        run_experiment(exp_params, tree_panel_mod, i, data_dir, mut_p = mut_p)
})

exp_params = readRDS(paste0(rep_dir, "exp_params.rds"))
tree_panel_mod = readRDS(paste0(rep_dir, "tree_panel_mod.rds"))

dir.create(paste0(rep_dir, "res/"))
walk_ice_fase_data(exp_params,
                   input_dir = data_dir,
                   tree_panel_mod,
                   out_dir = paste0(rep_dir, "res/"))

dir.create(paste0(rep_dir, "res_gr/"))
walk_ice_fase_data_gr(exp_params,
                      input_dir = data_dir,
                      tree_panel_mod,
                      out_dir = paste0(rep_dir, "res_gr/"))

dir.create(paste0(rep_dir, "tr3_res/"))
walk_ice_fase_data_tr3(exp_params,
                       input_dir = paste0(rep_dir, "tr3/"),
                       tree_panel_mod,
                       out_dir = paste0(rep_dir, "tr3_res/"))

dir.create(paste0(rep_dir, "tr3_res_gr/"))
walk_ice_fase_data_tr3_gr(exp_params,
                          input_dir = paste0(rep_dir, "tr3/"),
                          tree_panel_mod,
                          out_dir = paste0(rep_dir, "tr3_res_gr/"))

exp_params = gather_ice_fase(exp_params, tree_panel_mod, out_dir = paste0(rep_dir, "tr3_res/"))
exp_params = gather_ice_fase(exp_params, tree_panel_mod, out_dir = paste0(rep_dir, "res/"))

g0 = exp_params %>% ggline(x = "target_time",
                           y = "kc0",
                           add = "mean_se",
                           color = "sampling")
exp_params = gather_ice_fase(exp_params, tree_panel_mod, out_dir = paste0(rep_dir, "tr3_res_gr/"))
exp_params = gather_ice_fase(exp_params, tree_panel_mod, out_dir = paste0(rep_dir, "res_gr/"))
exp_params = load_exp_data(exp_params, data_dir)
exp_params = run_qfm_eval(exp_params, tree_panel_mod)

eval_tb = extract_evals(exp_params)
eval_tb$target_time = exp_params$target_time[eval_tb$exp_id]
saveRDS(eval_tb, file = paste0(rep_dir, "eval_true_gr_tb.rds"))

eval_tb = readRDS(paste0(rep_dir, "eval_true_gr_tb.rds")) # true phylogeny
eval_tb = readRDS(paste0(rep_dir, "eval_gr_tb.rds")) # phylotime phylogeny

g1 = produce_strat_time_sum(eval_tb, c("target_time", "sampling", "exp_id")) %>%
        ggline(x = "target_time", y = "rmse", color = "sampling", add = "mean_se")
g2 = produce_strat_size_sum(eval_tb, c("target_time", "sampling", "exp_id")) %>%
        ggline(x = "target_time", y = "rmse", color = "sampling", add = "mean_se")
g3 = produce_strat_split_sum(eval_tb, c("target_time", "sampling", "exp_id")) %>%
        ggline(x = "target_time", y = "rmse", color = "sampling", add = "mean_se")

(g0 + g1 + g2 + g3 + plot_layout(guide = "collect", nrow = 1) &
        theme(legend.position = "bottom", text = element_text(size = 10))) %>%
        # push_pdf("multi_time_eval", w = 6, h= 2.8, dir = plot_dir)
        push_pdf("multi_time_tr3_eval", w = 6, h= 2.8, dir = plot_dir)







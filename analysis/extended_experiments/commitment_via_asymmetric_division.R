library(qfm)
library(furrr)
plan(multisession, workers = 12)

output_dir = "./intermediate_data/panel_mod2_v1/"
plot_dir = "./plots/panel_mod2_v1_examples/"
# setting up
# dir.create(plot_dir)
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))
type_graph = tree_panel$type_graph[[36]]
temp_node_label_mapped = function(x) {
        out = node_label_mapper(type_graph$all_id)
        names(out) = type_graph$all_id
        for (i in 2:9) {
                out[[paste0(-i)]] = paste0("T", i+1)
        }
        out[["-9"]] = "T2"
        out[x]
}
all_id_mapped = temp_node_label_mapped(type_graph$all_id)

type_graph_mod1 = rlang::duplicate(type_graph)
type_graph_mod2 = rlang::duplicate(type_graph)
for (i in 1:16) {
        x = paste0("T", 1:16)[i]
        mapped_x = names(all_id_mapped)[all_id_mapped == x]
        type_graph_mod1$lifetime[[mapped_x]] = 0.6 - (0.6 - 0.1)/16 * i
        type_graph_mod2$lifetime[[mapped_x]] = 0.6 - (0.6 - 0.1)/16 * i
}


type_graph_mod1$lifetime
for (x in type_graph_mod1$node_id) {
        type_graph_mod1$diff_mode_probs[[x]] = c(0, 0, 1)
}
for (x in type_graph_mod2$node_id) {
        type_graph_mod2$diff_mode_probs[[x]] = c(0.5, 0.5, 0)
}
# for (x in type_graph_mod1$tip_id) {
#         type_graph_mod1$prob_pm[[x]] = 1
# }
# for (x in type_graph_mod2$tip_id) {
#         type_graph_mod2$prob_pm[[x]] = 1
# }
# type_graph_mod1$lifetime[[type_graph_mod1$root_id]] = 0.6/4
# type_graph_mod2$lifetime[[type_graph_mod1$root_id]] = 0.6/4

for (x in type_graph_mod1$node_id) {
        type_graph_mod1$prob_loss[[x]] = 0
}
for (x in type_graph_mod2$node_id) {
        type_graph_mod2$prob_loss[[x]] = 0
}
for (x in type_graph_mod1$tip_id) {
        type_graph_mod1$prob_loss[[x]] = 0.1
}
for (x in type_graph_mod2$tip_id) {
        type_graph_mod2$prob_loss[[x]] = 0.1
}
for (x in type_graph_mod1$tip_id) {
        type_graph_mod1$prob_pm[[x]] = 0.5
        type_graph_mod2$prob_pm[[x]] = 0.5
}
set.seed(73)
gr_col = gr_color_v1(type_graph_mod1)
plot_type_graph_clean_mod2(type_graph_mod1, node_col_mapper = function(x) gr_col[x], node_label_mapper = temp_node_label_mapped,
                           show_node_text = T, show_node_counts = T, ylim = c(12.5, 0),
                           add_founder = T) %>%
        push_pdf("g_assym", dir = plot_dir)
plot_type_graph_clean_mod2(type_graph_mod2, node_col_mapper = function(x) gr_col[x], node_label_mapper = temp_node_label_mapped,
                           show_node_text = T, show_node_counts = T, ylim = c(12.5, 0),
                           add_founder = T) # %>%
# push_pdf("g_assym", dir = plot_dir)

gens0 = make_gens_mod2(type_graph_mod1)
mut_p = readRDS("./intermediate_data/panel_mod2_v1/mut_p_marc1.rds")
exp_params = tibble(big_graph_id = c(rep(1, 100),
                                     rep(2, 100)),
                    sample_size = 50,
                    sampling = "fixed")
tree_panel_mod = tibble(type_graph = list(type_graph_mod1, type_graph_mod2))
tree_panel_mod$gr = map(tree_panel_mod$type_graph, as.phylo_mod2.type_graph)

data_out = simulate_sc_data_mod2(type_graph_mod1, mut_p, sample_size = 50, delay = T)

# loading results
rep_dir = "./intermediate_data/panel_mod2_v1/asymm/"
# dir.create(rep_dir)
# saveRDS(tree_panel_mod, file = paste0(rep_dir, "tree_panel_mod.rds"))
# saveRDS(exp_params, file = paste0(rep_dir, "exp_params.rds"))
# future_walk(which(exp_params$big_graph_id == 1), function(i) {
#         run_experiment(exp_params, tree_panel_mod, i, rep_dir, mut_p = mut_p, delay = T)
# })
# future_walk(which(exp_params$big_graph_id == 2), function(i) {
#         run_experiment(exp_params, tree_panel_mod, i, rep_dir, mut_p = mut_p, delay = F)
# })

tree_panel_mod = readRDS(paste0(rep_dir, "tree_panel_mod.rds"))
exp_params = readRDS(paste0(rep_dir, "exp_params.rds"))

# res_dir = paste0(rep_dir, "res/")
# dir.create(res_dir)
# # library(furrr)
# # plan(multisession, workers = 12)
# walk_ice_fase_data(exp_params, input_dir = rep_dir, tree_panel = tree_panel_mod, out_dir = res_dir)

res_dir = paste0(rep_dir, "res_tr3/")
# res_dir = paste0(rep_dir, "res_tr3_gr/")
# dir.create(res_dir)
library(furrr)
plan(multisession, workers = 12)
walk_ice_fase_data_tr3(exp_params, input_dir = paste0(rep_dir, "tr3/"), tree_panel = tree_panel_mod, out_dir = res_dir)
# walk_ice_fase_data_tr3_gr(exp_params, input_dir = paste0(rep_dir, "tr3/"), tree_panel = tree_panel_mod, out_dir = res_dir)

exp_params = gather_ice_fase(exp_params, tree_panel_mod, res_dir)
# exp_params$cond = tree_panel_mod$cond[exp_params$big_graph_id]
# exp_params$prob = tree_panel_mod$prob[exp_params$big_graph_id]

exp_params$cond = exp_params$big_graph_id
exp_params %>% group_by(cond) %>% summarise(mean_kc0 = mean(kc0),
                                            mean_rslv = mean(kc0 == 0))
exp_params %>% ggboxplot(x = "big_graph_id", y = "kc0")

exp_params = load_exp_data(exp_params, rep_dir)
exp_params_cond1 = exp_params[exp_params$big_graph_id == 1, ]
exp_params_cond1 = run_qfm_eval(exp_params_cond1, tree_panel_mod, delay = T)
exp_params_cond2 = exp_params[exp_params$big_graph_id == 2, ]
exp_params_cond2 = run_qfm_eval(exp_params_cond2, tree_panel_mod, delay = F)

eval_tb = extract_evals(bind_rows(exp_params_cond1, exp_params_cond2))
eval_tb$cond = exp_params$big_graph_id[eval_tb$exp_id]

eval_tb$cond = factor(eval_tb$cond)
eval_tb$gr_node_split_error = eval_tb$gr_node_split_order - eval_tb$node_split_order
# eval_tb %>% ggscatter(x = "log2_node_sampled", y = "gr_node_split_error", color = "cond") +
#         geom_smooth(aes(color = cond))

# eval_tb %>% ggscatter(x = "node_time", y = "gr_node_split_error", color = "cond") +
#         geom_smooth(aes(color = cond))
saveRDS(eval_tb, file = paste0(output_dir, "asymm_eval_gr_tb.rds"))

eval_tb = readRDS(paste0(output_dir, "asymm_eval_gr_tb.rds"))
eval_tb_true = extract_evals_true(bind_rows(exp_params_cond1, exp_params_cond2))
eval_tb_true$cond = exp_params$big_graph_id[eval_tb_true$exp_id]
ggviolin(eval_tb_true, y = "log2_node_sampled", x = "cond")

produce_strat_sampled_split_sum <- function(tb, group_by) {
        tb %>%
                group_by(.dots = map(group_by, as.symbol)) %>%
                filter(!is.na(gr_node_split_order) & !is.na(node_split_order)) %>%
                summarise(cc = cor(node_split_sampled_order, gr_node_split_order, method = "pearson"),
                          cs = cor(node_split_sampled_order, gr_node_split_order, method = "spearman"),
                          rmse = sqrt(mean((node_split_sampled_order - gr_node_split_order)^2)),
                          total = n())
}

temp_sampled_split = produce_strat_sampled_split_sum(eval_tb, c("cond", "exp_id"))
temp_sampled_split$cond = factor(temp_sampled_split$cond, labels = c("Asymm.", "Symm."))
temp_sampled_split %>%
        ggboxplot(x = "cond", y = "rmse", size = 0.4, outlier.size = 0.4,
                  ylab = "RMSE of sampled bias") + stat_compare_means() %>%
        push_png("asymm_sampled_bias", res = 500,
                 w = 2.5, h = 2.5,
                 dir = "./plots/panel_mod2_v1_examples/")

summarize_all_rmse(eval_tb, c("cond", "exp_id")) %>%
        write_csv(paste0(output_dir, "asymm_param_rmse.csv"))

# eval_tb$prob = exp_params$prob[eval_tb$exp_id]
((exp_params %>% ggboxplot(x = "big_graph_id", y = "kc0", size = 0.4, outlier.size = 0.4)) +
                stat_compare_means() +
produce_strat_time_sum(eval_tb, c("cond", "exp_id")) %>%  ggboxplot(x = "cond", y = "rmse", size = 0.4, outlier.size = 0.4) +
                stat_compare_means() +
        produce_strat_size_sum(eval_tb, c("cond", "exp_id")) %>% ggboxplot(x = "cond", y = "rmse", size= 0.4, outlier.size = 0.4) +
                stat_compare_means() +
        produce_strat_split_sum(eval_tb, c("cond", "exp_id")) %>% ggboxplot(x = "cond", y = "rmse", size = 0.4, outlier.size = 0.4) +
                stat_compare_means() +
        plot_layout(nrow = 2) & theme(text = element_text(size =10))) %>%
        push_pdf("asymm_eval_gr_tr3", w = 3, h= 3, dir = plot_dir)

gens0 = make_gens_mod2(type_graph_mod1)
gens1 = sample_size_gens_mod2(gens0, type_graph = type_graph_mod1, sample_size = 50)

rmse_true = map_dbl(exp_params$true_sampled_sizes, function(x_ls) {
        sqrt(mean((map_dbl(x_ls, function(x) x[2] / sum(x[2:3]))[type_graph_mod1$node_id] - 0.5)^2))
})
boxplot(cbind(rmse_true[1:100], rmse_true[101:200]))

tr = exp_params$tr[[2]]
plot_tr(tr, node_types = c(get_type_from_id(tr$node.label),
                           get_type_from_id(tr$tip.label)),
        type_col = gr_col, node_size = 3,
        ylim = c(12., -1),
        jitter_tip = T) %>%
        push_pdf("tr_asymm_jitter_v1", width = 20, h = 5, dir = plot_dir)

tr = exp_params$tr[[101]]
plot_tr(tr, node_types = c(get_type_from_id(tr$node.label),
                           get_type_from_id(tr$tip.label)),
        type_col = gr_col, node_size = 3,
        ylim = c(12., -1),
        jitter_tip = T) %>%
        push_pdf("tr_symm_jitter_v1", width = 20, h = 5, dir = plot_dir)


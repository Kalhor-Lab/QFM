library(qfm)
output_dir = "./intermediate_data/panel_mod2_v1/"
plot_dir = "./plots/panel_mod2_v1_examples/"
# dir.create(plot_dir)
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))

# condition 1: pick progenitor states to increase death rate
# condition 2: pick terminal states to increase death rate
type_graph = tree_panel$type_graph[[17]]
set.seed(73)

adjust_death_rate <- function(type_graph, prob, cond = c("tip", "node")) {
        out = rlang::duplicate(type_graph)
        true_sizes = map_dbl(get_true_size_mod2(type_graph), 1)
        node_mod = type_graph$node_id[true_sizes[type_graph$node_id] > 500]
        tip_mod = type_graph$tip_id[true_sizes[type_graph$tip_id] > 500]
        if (cond == "node") {
                out$prob_loss[node_mod] = prob
        }
        if (cond == "tip") {
                out$prob_loss[tip_mod] = prob
        }
        out
}
tree_panel_mod = as_tibble(expand.grid(cond = c("node", "tip"),
                                       prob = seq(0.1, 0.4, by = 0.05),
                                       big_graph_id = which(tree_panel$num_tip == 16)))
tree_panel_mod$type_graph = pmap(tree_panel_mod, function(cond, prob, big_graph_id) {
        adjust_death_rate(tree_panel$type_graph[[big_graph_id]],
                          prob = prob,
                          cond = cond)
})
# tree_panel_mod$type_graph = map2(tree_panel_mod$cond, tree_panel_mod$prob, function(cond, p) {
#         out = rlang::duplicate(type_graph)
#         if (cond == "node") {
#                 out$prob_loss[node_mod] = p
#         } else {
#                 out$prob_loss[tip_mod] = p
#         }
#         out
# })
tree_panel_mod$gr = map(tree_panel_mod$type_graph, as.phylo_mod2.type_graph)

node_label_mapped = node_label_mapper(tree_panel_mod$type_graph[[1]]$all_id)
names(node_label_mapped) = tree_panel_mod$type_graph[[1]]$all_id
gr_col = gr_color_v1(tree_panel_mod$type_graph[[1]])

plot_type_graph_clean_mod2(tree_panel_mod$type_graph[[1]],
                           function(x) gr_col[x],
                           node_label_mapper = node_label_mapper,
                           show_node_text = T,add_founder = T, show_node_counts = T,
                           ylim = c(12.5, -1)) %>%
        push_pdf("death_rate_g_v1", w = 4, h = 4, dir = plot_dir)

tree_panel_mod$cond[2]
node_label_mapped[which(unlist(tree_panel_mod$type_graph[[1]]$prob_loss) == tree_panel_mod$prob[1])]
node_label_mapped[which(unlist(tree_panel_mod$type_graph[[2]]$prob_loss) == tree_panel_mod$prob[2])]

# do.call(rbind, map(tree_panel_mod$type_graph, function(x) {
#         print('aa')
#         map_dbl(get_true_size_mod2(x)[tree_panel$type_graph[[81]]$all_id], 1)
# }))

exp_params = expand.grid(big_graph_id = 1:nrow(tree_panel_mod),
                         sample_size = 100,
                         num_element = c(50),
                         sampling = c("fixed"),
                         i_sim = 1) %>% as_tibble()
exp_params$target_time = map_dbl(exp_params$big_graph_id, function(graph_id) tree_panel_mod$type_graph[[graph_id]]$target_time)

# loading and running
mut_p = readRDS("./intermediate_data/panel_mod2_v1/mut_p_marc1.rds")
rep_dir = paste0("./intermediate_data/panel_mod2_v1/", "death_rate_v4/")
# dir.create(rep_dir)
# saveRDS(exp_params, file = paste0(rep_dir, "exp_params.rds"))
# saveRDS(tree_panel_mod, file = paste0(rep_dir, "tree_panel_mod.rds"))

exp_params = readRDS(paste0(rep_dir, "exp_params.rds"))
tree_panel_mod = readRDS(paste0(rep_dir, "tree_panel_mod.rds"))

# library(furrr)
# plan(multisession, workers = 12)
# future_walk(1:nrow(exp_params), function(i) {
#         run_experiment(exp_params, tree_panel_mod, i = i,
#                        mut_p = mut_p, rep_dir)
# })

library(furrr)
plan(multisession, workers = 12)
# res_dir = paste0(rep_dir, "res/")
# dir.create(res_dir)
# walk_ice_fase_data(exp_params, input_dir = rep_dir, tree_panel = tree_panel_mod, out_dir = res_dir)

tr3_dir = paste0(rep_dir, "tr3/")
# dir.create(tr3_dir)
# walk_phylotime_data(exp_params, input_dir = rep_dir, tree_panel = tree_panel_mod, mut_p = mut_p, out_dir = tr3_dir)
# res_dir = paste0(rep_dir, "res_tr3/")
# dir.create(res_dir)
res_dir = paste0(rep_dir, "res_tr3_gr/")
# dir.create(res_dir)
# walk_ice_fase_data_tr3(exp_params, input_dir = tr3_dir, tree_panel = tree_panel_mod, out_dir = res_dir)
# walk_ice_fase_data_tr3(exp_params, input_dir = tr3_dir, tree_panel = tree_panel_mod, out_dir = res_dir)
walk_ice_fase_data_tr3_gr(exp_params, input_dir = tr3_dir, tree_panel = tree_panel_mod, out_dir = res_dir)

# walk_ice_fase_gr(exp_params, tr_col = "tr", tree_panel = tree_panel_mod, out_dir = res_dir)
exp_params = gather_ice_fase(exp_params, tree_panel_mod, res_dir)
exp_params$cond = tree_panel_mod$cond[exp_params$big_graph_id]
exp_params$prob = tree_panel_mod$prob[exp_params$big_graph_id]
exp_params = load_exp_data(exp_params, rep_dir)
exp_params = run_qfm_eval(exp_params, tree_panel_mod)
# saveRDS(select(exp_params, -c(sc, tr, true_sampled_sizes)), file = paste0(rep_dir, "all_eval.rds"))
saveRDS(select(exp_params, -c(sc, tr, true_sampled_sizes)), file = paste0(rep_dir, "all_eval_gr.rds"))

exp_params %>% group_by(cond, prob) %>% summarise(mean_kc0 = mean(kc0))
exp_params %>% ggboxplot(x = "prob", color = "cond", y = "kc0")
eval_tb = extract_evals(exp_params)

eval_tb$cond = exp_params$cond[eval_tb$exp_id]
eval_tb$prob = exp_params$prob[eval_tb$exp_id]

# ((exp_params %>% ggboxplot(x = "prob", fill = "cond", y = "kc0") %>%
#           facet("cond", nrow = 2) + stat_compare_means())+
#                 (produce_strat_time_sum(eval_tb) %>%
#          ggbarplot(x = "prob", y = "rmse", fill = "cond") %>%
#                  facet(facet.by = "cond", nrow = 2)) +
# (produce_strat_size_sum(eval_tb) %>%
#          ggbarplot(x = "prob", y = "rmse", fill = "cond") %>%
#          facet(facet.by = "cond", nrow = 2)) +
# (produce_strat_split_sum(eval_tb) %>%
#          ggbarplot(x = "prob", y = "rmse", fill = "cond") %>%
#          facet(facet.by = "cond", nrow = 2)) +
# plot_layout(nrow = 1, guide = "collect") & theme(legend.position = "bottom")) %>%
#         push_pdf(file_name = "death_rate_param", dir = plot_dir, w = 6, h = 4)

# need to replace kc0 in exp_params with actual reconstructions
g0 = (exp_params %>% ggline(x = "prob", color = "cond", y = "kc0", add = "mean_se"))
g1 = (produce_strat_time_sum(eval_tb, c("cond", "prob", "exp_id")) %>% ggline(x = "prob", y = "rmse", add = "mean_se", color = "cond"))
g2 = (produce_strat_size_sum(eval_tb, c("cond", "prob", "exp_id")) %>% ggline(x = "prob", y = "rmse", add = "mean_se", color = "cond"))
g3 = (produce_strat_split_sum(eval_tb, c("cond", "prob", "exp_id")) %>% ggline(x = "prob", y = "rmse", add = "mean_se", color = "cond"))
(((g0 + g1 + g2 + g3 + g_n1) & scale_color_manual(values = c("#0b97d9", "#cc5500"))) +
        plot_layout(nrow = 1, guide = "collect") &
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))) %>%
        push_pdf(file_name = "death_rate_param_tr3_gr_new", dir = plot_dir, w = 9., h = 3)

eval_tb_true = extract_evals_true(exp_params)
eval_tb_true$prob = exp_params$prob[eval_tb_true$exp_id]
eval_tb_true$cond = exp_params$cond[eval_tb_true$exp_id]
g_n1 = (eval_tb_true %>%
        group_by(cond, prob, exp_id) %>%
        mutate(frac_undersampled = mean(log2_node_sampled < -2)) %>%
        ggline(x = "prob", color = "cond", y = "frac_undersampled", add = "mean_se",
                  ylab = "Fraction of progenitor states undersampled"))# %>%
        # push_pdf("death_rate_frac_undersampled_gr", dir = plot_dir, w = 2., h = 3.)

eval_tb_true$node_split_sampled_abs_error = abs(eval_tb_true$node_split_sampled_order - eval_tb_true$node_split_order)

ggboxplot(eval_tb_true, x = "prob", color = "cond", y = "node_split_sampled_abs_error")
group_by = c("prob", "cond")
eval_tb_true %>%
        group_by(.dots = map(group_by, as.symbol)) %>%
        filter(!is.na(gr_time_trans) & !is.na(node_time)) %>%
        summarize(sd(node_split_sampled_order))




output_dir = "./intermediate_data/panel_mod2_v1/"

# collecting KC0 results
exp_params$kc0_tr3 = map_dbl(exp_params$res3, function(x) x$gr_eval$kc0)
exp_params$kc0_control = map_dbl(exp_params$gr_control_eval, "kc0")
select(exp_params, big_graph_id, sample_size, num_element, sampling, i_sim, kc0, kc0_tr3, bsum, num_tip, kc0_control) %>%
        write_csv(paste0(output_dir, "exp_stat_1rep.csv"))

readout_exp = tibble(dir = c("./intermediate_data/panel_mod2_v1/dropout_a1/",
                             "./intermediate_data/panel_mod2_v1/dropout/",
                             "./intermediate_data/panel_mod2_v1/exist_error_allele_v1/"
),
mode = c("dropout", "dropout", "allele_error"))
exp_params = read_csv(paste0(output_dir, "exp_stat_1rep.csv"))
exp_params$exp_id = 1:nrow(exp_params)
all_stat_tb = bind_rows(c(map2(readout_exp$dir, readout_exp$mode, function(rep_dir, m) {
        read_csv(paste0(rep_dir, "all_stat.csv")) %>% mutate(mode = m)
}),
list(transmute(exp_params, frac = 0, mode = "dropout", kc0 = kc0_tr3, exp_id, big_graph_id)),
list(transmute(exp_params, frac = 0, mode = "allele_error", kc0 = kc0_tr3, exp_id, big_graph_id))
))
all_stat_summary = all_stat_tb %>% group_by(frac, mode) %>%
        summarize(mean_rslv = mean(kc0 == 0),
                  mean_kc0 = mean(kc0),
                  median_kc0 = median(kc0))
write_csv(all_stat_summary, file = paste0(output_dir, "readout_kc0_stats_v1.csv"))

g0 = all_stat_tb %>%
        # mutate(mode = factor(mode)) %>%
        ggline(x = "frac", y = "kc0", color = "mode", add = "mean_se", size = 0.25,
               ylim = c(0, NA)) + geom_hline(yintercept = mean(exp_params$kc0_control),
                                             linetype = 2)

# load processed results from readout errors
dropout_eval = readRDS(paste0("./intermediate_data/panel_mod2_v1/dropout/", "eval_tb_gr.rds"))
dropout_a1_eval = readRDS(paste0("./intermediate_data/panel_mod2_v1/dropout_a1/", "eval_tb_gr.rds"))
error_eval = readRDS(paste0("./intermediate_data/panel_mod2_v1/exist_error_allele_v1//", "eval_tb_gr.rds"))
load(paste0(output_dir, "exp_eval_true_tr_gr.rds"))
# load(paste0(output_dir, "exp_eval_true_tr_gr.rds"))
eval_tb_true$frac = -1
eval_tb_true_tr$frac = 0
dropout_a1_eval$frac = as.numeric(as.character(dropout_a1_eval$frac))

eval_tb_true_gr$frac = -1
eval_tb_true_tr_gr$frac = 0

# making exp_id unique
dropout_eval$exp_id = paste0("d_", dropout_eval$exp_id)
dropout_a1_eval$exp_id = paste0("a_", dropout_a1_eval$exp_id)
error_eval$exp_id = paste0("e_", error_eval$exp_id)
eval_tb_true_gr$exp_id = paste0("t_", eval_tb_true_gr$exp_id)
eval_tb_true_tr_gr$exp_id = paste0("tr_", eval_tb_true_tr_gr$exp_id)

tb_template_v3 <- function(sum_func, metric = "rmse") {
        group_by = c("exp_id", "frac")
        bind_rows(list(mutate(sum_func(dropout_eval, group_by), mode = "dropout"),
                        mutate(sum_func(dropout_a1_eval, group_by), mode = "dropout"),
                        mutate(sum_func(error_eval, group_by), mode = "allele_error"),
                        mutate(sum_func(eval_tb_true_tr_gr, group_by), mode = "truth"))
        )
}
tb_sum = bind_rows(list(time = tb_template_v3(produce_strat_time_sum) %>% group_by(mode, frac) %>%
             summarise(mean_rmse = mean(rmse),
                       type = "time"),
     size = tb_template_v3(produce_strat_size_sum) %>% group_by(mode, frac) %>%
             summarise(mean_rmse = mean(rmse),
                       type = "size"),
     split = tb_template_v3(produce_strat_split_sum) %>% group_by(mode, frac) %>%
             summarise(mean_rmse = mean(rmse),
                       type = "split")
     ))
write_csv(tb_sum, paste0(output_dir, "readout_summary_stat_v1.csv"))
# version without mis-labeling
plot_template_v4 <- function(sum_func, metric = "rmse") {
        group_by = c("exp_id", "frac")
        (bind_rows(list(mutate(sum_func(dropout_eval, group_by), mode = "dropout"),
                        mutate(sum_func(dropout_a1_eval, group_by), mode = "dropout"),
                        mutate(sum_func(error_eval, group_by), mode = "allele_error"),
                        mutate(sum_func(eval_tb_true_tr_gr, group_by), mode = "dropout"),
                        mutate(sum_func(eval_tb_true_tr_gr, group_by), mode = "allele_error")
        )) %>%
                        ggline(x = "frac", y = metric, color = "mode", size = 0.1, ylim = c(0, NA),
                               add = "mean_se")
        )
}

g1 = plot_template_v4(produce_strat_time_sum, metric = "rmse") + ggtitle("Time")
g2 = plot_template_v4(produce_strat_size_sum, metric ="rmse") + ggtitle("Size")
g3 = plot_template_v4(produce_strat_split_sum,  metric = "rmse") + ggtitle("Bias")
((g0 + ggtitle("Topoogy")) + g1 + g2 + g3 + plot_layout(guide = "collect", nrow = 1) & theme(legend.position = "bottom")) %>%
        push_pdf("readout_evals_v1", w = 15, h = 5, dir = "./plots/panel_mod2_v1/")





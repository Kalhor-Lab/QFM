# this script was used to process readout results

library(qfm)
library(furrr)
plan(multisession, workers = 12)
# Experiments for errorous readouts
# based on simulated panels
output_dir = "./intermediate_data/panel_mod2_v1/"
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))

# summarize processed results with true phylogeny and Phylotime phylogeny without error
exp_params = readRDS(paste0(output_dir, "exp_data_1rep_proc_v1.rds"))
exp_params$true_sampled_sizes = map(exp_params$data, "true_sampled_sizes")
exp_params$res_gr = future_map2(exp_params$tr, exp_params$big_graph_id, function(tr, graph_id) {
        res = ice_fase_mod_gr(tr,
                              get_type_from_id(tr$tip.label),
                              tree_panel$gr[[graph_id]],
                              total_time = 11.5 - 0.6,
                              root_time = 0.6,
                              theta = 0.0)
        res = evaluate_gr_res(res, tree_panel$gr[[graph_id]])
        res
}, .progress = T, .options = furrr_options(seed = T))
exp_params = run_qfm_eval(exp_params, tree_panel, "res_gr", "qfm_eval_gr")

exp_params$res_gr = future_map2(exp_params$tr, exp_params$big_graph_id, function(tr, graph_id) {
        res = ice_fase_mod_gr(tr,
                              get_type_from_id(tr$tip.label),
                              tree_panel$gr[[graph_id]],
                              total_time = 11.5 - 0.6,
                              root_time = 0.6,
                              theta = 0.0)
        res = evaluate_gr_res(res, tree_panel$gr[[graph_id]])
        res
}, .progress = T, .options = furrr_options(seed = T))
exp_params = run_qfm_eval(exp_params, tree_panel, "res_gr", "qfm_eval_gr")
exp_params$res3_gr = future_map2(exp_params$tr3, exp_params$big_graph_id, function(tr, graph_id) {
        res = ice_fase_mod_gr(tr,
                              get_type_from_id(tr$tip.label),
                              tree_panel$gr[[graph_id]],
                              total_time = 11.5 - 0.6,
                              root_time = 0.6,
                              theta = 0.0)
        res = evaluate_gr_res(res, tree_panel$gr[[graph_id]])
        res
}, .progress = T, .options = furrr_options(seed = T))
exp_params = run_qfm_eval(exp_params, tree_panel, "res3_gr", "qfm_eval3_gr")
eval_tb_true_gr = extract_evals(exp_params, eval_col = "qfm_eval_gr")
eval_tb_true_tr_gr = extract_evals(exp_params, eval_col =  "qfm_eval3_gr")
save(eval_tb_true_gr, eval_tb_true_tr_gr, file = paste0(output_dir, "exp_eval_true_tr_gr.rds"))

gr_eval_list = list(true_phylogeny = map(exp_params$res, "gr_eval"),
                    phylotime_phylogeny = map(exp_params$res3, "gr_eval"))
saveRDS(gr_eval_list, file = paste0(output_dir, "gr_eval.rds"))

exp_params$tr3 = map(1:nrow(exp_params), function(i) {
        out_file = paste0(output_dir, "tr3/", str_pad(i, width = 4, pad = "0"), ".rds")
        if (file.exists(out_file)) {
                return(readRDS(out_file))
        } else {
                return(NULL)
        }
})
exp_params$res3 = future_map(exp_params$tr3, function(tr) {
        res = ice_fase_mod(tr,
                           get_type_from_id(tr$tip.label),
                           total_time = 11.5 - 0.6,
                           root_time = 0.6,
                           theta = 0.)
        res
}, .progress = T, .options = furrr_options(seed = T))
exp_params$res3 = map2(exp_params$big_graph_id, exp_params$res3, function(i, res) {
        res = evaluate_gr_res(res, tree_panel$gr[[i]])
        res
})
exp_params = run_qfm_eval(exp_params, tree_panel, "res3", "qfm_eval3")
saveRDS(exp_params, file = paste0(output_dir, "exp_data_1rep_proc_v1.rds"))

eval_tb_true = extract_evals(exp_params)
eval_tb_true_tr = extract_evals(exp_params, eval_col =  "qfm_eval3")
save(eval_tb_true, eval_tb_true_tr, file = paste0(output_dir, "exp_eval_true_tr.rds"))
# end summarizing results

exp_params = read_csv(paste0(output_dir, "exp_stat_1rep.csv"))
exp_params$exp_id = as.numeric(1:nrow(exp_params))
exp_params$frac = 0.
exp_params$target_time = map_dbl(exp_params$big_graph_id, function(graph_id) tree_panel$type_graph[[graph_id]]$target_time)

# repeat below processing for two types of readout error
rep_dir = paste0(output_dir, "dropout/")
rep_dir = paste0(output_dir, "exist_error_allele_v1/")
# adding information to dropout_all_tb
dropout_all_tb = bind_rows(map(1:nrow(exp_params), function(i) {
        message(i)
        dropout_tb = readRDS(paste0(rep_dir, str_pad(i, width = 4, pad = "0"), ".rds"))
        dropout_tb$exp_id = i
        dropout_tb = left_join(dropout_tb,
                               select(exp_params, exp_id,
                                      big_graph_id, sample_size, num_element, sampling, i_sim,
                                      bsum, num_tip),
                               by = "exp_id")
        dropout_tb = dplyr::rename(dropout_tb, tr3 = tr)
}))
# dropout_all_tb = bind_rows(list(dropout_all_tb,
#                                 select(exp_params,
#                                        tr3, mut_p, sc,
#                                        big_graph_id, sample_size, num_element, sampling, i_sim,
#                                        exp_id, target_time,
#                                        bsum, num_tip, frac)))
# dropout_all_tb = arrange(dropout_all_tb, exp_id, frac)
saveRDS(dropout_all_tb, file = paste0(rep_dir, "all_tb.rds"))

dropout_all_tb = readRDS(paste0(rep_dir, "all_tb.rds"))
# save a base version
dropout_base_tb = dropout_all_tb[c(1, 4, 6, 7, 8, 9, 10, 11, 12, 13, 15)]
# dropout_base_tb = dropout_all_tb[c(1, 4, 6, 7, 8, 9, 10, 11, 12, 13)]
saveRDS(dropout_base_tb, file = paste0(rep_dir, "all_tb_base.rds"))
dump_tr(dropout_all_tb, tr_col = "tr3", paste0(rep_dir, "tr3/"))

dropout_all_tb = readRDS(paste0(rep_dir, "all_tb_base.rds"))
# res_dir = paste0(rep_dir, "res_gr/")
# res_dir = paste0(rep_dir, "res/")
# dir.create(res_dir)
# options(future.globals.maxSize= 5e9)
# walk_ice_fase_data_tr3(select(dropout_all_tb, big_graph_id), paste0(rep_dir, "tr3/"), tree_panel, res_dir)
# dropout_all_tb = gather_ice_fase(dropout_all_tb, tree_panel, res_dir)

dump_tr(dropout_all_tb, tr_col = "tr3", paste0(rep_dir, "tr3/"))
# walk_ice_fase(dropout_all_tb, tree_panel, "tr3", res_dir)
# saveRDS(dropout_all_tb, file = paste0(rep_dir, "all_tb.rds"))

res_dir = paste0(rep_dir, "res/")
# dir.create(res_dir)
walk_ice_fase_data_tr3(select(dropout_all_tb, big_graph_id), paste0(rep_dir, "tr3/"), tree_panel, res_dir)

res_dir = paste0(rep_dir, "res_gr/")
# dir.create(res_dir)
walk_ice_fase_data_tr3_gr(select(dropout_all_tb, big_graph_id), paste0(rep_dir, "tr3/"), tree_panel, res_dir)

# dir.create(res_dir)
error_col = RColorBrewer::brewer.pal(6, "PuRd")
dropout_all_tb$frac = factor(dropout_all_tb$frac)
walk_ice_fase_data_tr3_gr(select(dropout_all_tb, big_graph_id), paste0(rep_dir, "tr3/"), tree_panel, res_dir)
dropout_all_tb = gather_ice_fase(dropout_all_tb, tree_panel, res_dir)

ggboxplot(dropout_all_tb, x = "frac", y = "kc0",
          facet.by = c("sampling", "num_tip"),
          fill = "frac") +
        scale_fill_manual(values = error_col)

# adding true sampled sizes from simulation results for evaluation later
exp_true_sampled_sizes = map(exp_params$data, "true_sampled_sizes")
saveRDS(exp_true_sampled_sizes, file = paste0(output_dir, "temp_exp_true_sampled_sizes.rds"))
exp_true_sampled_sizes = readRDS(paste0(output_dir, "temp_exp_true_sampled_sizes.rds"))
dropout_all_tb$true_sampled_sizes = exp_true_sampled_sizes[dropout_all_tb$exp_id]

dropout_all_tb = run_qfm_eval(dropout_all_tb, tree_panel)
# dropout_all_eval = select(dropout_all_tb, -c(sc, mut_p, tr3))
# saveRDS(dropout_all_eval, file = paste0(rep_dir, "all_eval_gr.rds")) # RUN THIS
saveRDS(dropout_all_tb, file = paste0(rep_dir, "all_eval_gr.rds")) # RUN THIS
# exp_params$target_time = map_dbl(exp_params$big_graph_id, function(graph_id) tree_panel$type_graph[[graph_id]]$target_time)
# dropout_all_eval = readRDS(paste0(rep_dir, "all_eval.rds"))
eval_tb = extract_evals(dropout_all_tb)
eval_tb$frac = dropout_all_tb$frac[eval_tb$exp_id]

# saveRDS(eval_tb, file = paste0(rep_dir, "eval_tb.rds"))
saveRDS(eval_tb, file = paste0(rep_dir, "eval_tb_gr.rds"))


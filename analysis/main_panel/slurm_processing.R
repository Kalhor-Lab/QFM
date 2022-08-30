job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(qfm)
output_dir = "./output/"
# exp_name = "5rep_mut_all"
exp_name = "20rep_16type"
# panel of fate maps
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))
# set of barcoding elements
# mut_p_all = readRDS(paste0(output_dir, "mut_p_all.rds"))
mut_p_all = readRDS(paste0(output_dir, "mut_p_marc1.rds"))
exp_params = readRDS(paste0(output_dir, "exp_data_", exp_name, ".rds"))

rep_dir = paste0(output_dir, exp_name, "/")
if (!dir.exists(rep_dir)) {
        dir.create(rep_dir)
}
run_experiment(exp_params, tree_panel, job_id, rep_dir, mut_p = mut_p_all)

# example for running ice_fase below
job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
batch_assign = ceiling((1:3310)/100)
batch_indices = which(batch_assign == job_id)

tr3_dir = paste0(output_dir, "tr3_50site/")
res_dir = paste0(output_dir, "res_50site/")
if (!dir.exists(res_dir)) {
        dir.create(res_dir)
}
for (i in batch_indices) {
        message(i)
        data_file = paste0(tr3_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        tr3 = readRDS(data_file)
        out_file = paste0(res_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        if (!file.exists(out_file)) {
                tr = tr3
                graph_id = exp_params$big_graph_id[i]
                tt = tree_panel$type_graph[[graph_id]]$target_time
                # gr_in = tree_panel$gr[[graph_id]]
                gr_in = NULL
                res = ice_fase_mod(tr,
                                   get_type_from_id(tr$tip.label),
                                   total_time = tt - 0.6,
                                   root_time = 0.6,
                                   theta = 0.,
                                   gr = gr_in)
                saveRDS(res, out_file)
        }
}


job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
batch_assign = ceiling((1:3310)/100)
batch_indices = which(batch_assign == job_id)

library(qfm)
source("../R_mod/mod2_v2.R")
source("../R_mod/ice_fase_mod1.R")
source("../R_mod/plotting.R")
source("../R_mod/exp_pipeline.R")

# output_dir = "../phylotime_runs/panel_mod2_v1/5rep_mut_all/"
# tree_panel = readRDS("./output/tree_panel.rds")
# mut_flat_tb = readRDS(paste0(output_dir, "mut_flat_all.rds"))
# exp_params = readRDS("./output/exp_data_5rep_mut_all.rds")
# mut_flat_tb$big_graph_id = exp_params$big_graph_id[mut_flat_tb$exp_id]

# tr3_dir = paste0(output_dir, "tr3/")
# res_dir = paste0(output_dir, "tr3_res/")
# res_dir = paste0(output_dir, "tr3_res_gr/")
# eval_dir = paste0(output_dir, "tr3_eval_gr/")
# data_dir = "./output/5rep_mut_all/"

# for (i in batch_indices) {
# 	message(i)
# 	data_obj = readRDS(paste0(data_dir, stringr::str_pad(mut_flat_tb$exp_id[i], width = 4, pad = "0"), ".rds"))
# 	res = readRDS(paste0(res_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
# 	out_file = paste0(eval_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
# 	if (!file.exists(out_file)) {
# 	        graph_id = mut_flat_tb$big_graph_id[i]
# 	        eval_obj = evaluate_qfm_v1(res, tree_panel$type_graph[[graph_id]], data_obj$true_sampled_sizes)
# 	        saveRDS(eval_obj, out_file)
# 	}
# }

# examples running evaluations on results
job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
batch_assign = ceiling((1:3310)/100)
batch_indices = which(batch_assign == job_id)

output_dir = "../phylotime_runs/panel_mod2_v1/5rep_mut_all/"
tree_panel = readRDS("./output/tree_panel.rds")
exp_params = readRDS("./output/exp_data_5rep_mut_all.rds")
data_dir = "./output/5rep_mut_all/"

res_dir = paste0(output_dir, "res/")
eval_dir = paste0(output_dir, "eval/")
if (!dir.exists(eval_dir)) {
        dir.create(eval_dir)
}
for (i in batch_indices) {
        message(i)
        data_obj = readRDS(paste0(data_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
        res = readRDS(paste0(res_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
        out_file = paste0(eval_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        if (!file.exists(out_file)) {
                graph_id = exp_params$big_graph_id[i]
                eval_obj = evaluate_qfm_v1(res, tree_panel$type_graph[[graph_id]], data_obj$true_sampled_sizes)
                eval_obj$gr_eval = evalute_gr(res$gr, tree_panel$gr[[graph_id]])
                saveRDS(eval_obj, out_file)
        }
}

# example summarizing evaluations, different version whether breaking down by sufficiently sampled
options(dplyr.summarise.inform = FALSE)
output_dir = "../phylotime_runs/panel_mod2_v1/5rep_mut_all/"
tree_panel = readRDS("./output/tree_panel.rds")
exp_params = readRDS("./output/exp_data_5rep_mut_all.rds")
# eval_dir = paste0(output_dir, "eval_50site/")
eval_dir = paste0(output_dir, "eval/")
log2_sampling_cutoff = -2
# eval_summary = map(1:nrow(exp_params), function(i) {
eval_summary = map(1:nrow(exp_params), function(i) {
        if (i %% 100 == 0) {
                message(i)
        }
        eval_obj = readRDS(paste0(eval_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
        eval_tb = mutate(eval_obj$recon, exp_id = i)
        eval_tb$suff_sampled = eval_tb$log2_node_sampled > log2_sampling_cutoff
        eval_tb = eval_tb %>%
                mutate(log2_node_size = log2(node_size),
                       log2_gr_node_size_in = log2(gr_node_size_in))
        eval_tb$gr_time_trans_abs_error = abs(eval_tb$gr_time_trans_error)
        eval_tb$gr_node_size_abs_logfc = abs(eval_tb$gr_node_size_logfc)
        eval_tb$gr_node_split_abs_error = abs(eval_tb$gr_node_split_order - eval_tb$node_split_order)
        # list(time = produce_strat_time_sum(eval_tb, "exp_id"),
        #         size = produce_strat_size_sum(eval_tb, "exp_id"),
        #                 split = produce_strat_split_sum(eval_tb, "exp_id"))
        list(time = produce_strat_time_sum(eval_tb, c("exp_id", "suff_sampled")) %>% ungroup(),
             size = produce_strat_size_sum(eval_tb, c("exp_id", "suff_sampled")) %>% ungroup(),
             split = produce_strat_split_sum(eval_tb, c("exp_id", "suff_sampled")) %>% ungroup())
})
# saveRDS(eval_summary, file = "eval/5rep_50site_suff_eval_tb.rds")
saveRDS(eval_summary, file = "eval/5rep_suff_eval_tb.rds")

# pooling KC0 results
eval_kc0 = map_dbl(1:nrow(exp_params), function(i) {
        if (i %% 100 == 0) {
                message(i)
        }
        eval_obj = readRDS(paste0(eval_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
        eval_obj$gr_eval$kc0
})
saveRDS(eval_kc0, file = "eval/5rep_50site_eval_kc0.rds")
# end example summarizing evaluations



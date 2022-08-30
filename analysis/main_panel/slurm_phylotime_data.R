# same script used for all experiments
library(qfm)
library(furrr)
n_cores = 12
plan(multisession, workers = n_cores)

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
job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

rep_dir = "./panel_mod2_v1/death_rate/"
exp_params = readRDS(paste0(rep_dir, "exp_params.rds"))
tree_panel = readRDS(paste0(rep_dir, "tree_panel_mod.rds"))
tr3_dir = paste0(rep_dir, "tr3/")
mut_p = readRDS("./panel_mod2_v1/mut_p_marc1.rds")

if (!dir.exists(tr3_dir)) {
        dir.create(tr3_dir)
}
mut_p = readRDS("./panel_mod2_v1/mut_p_marc1.rds")
data_dir = paste0(rep_dir, "data/")
run_phylotime_data(exp_params, job_id, data_dir, tree_panel, mut_p, tr3_dir)


job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(qfm)
library(furrr)
n_cores = 12
plan(multisession, workers = n_cores)

# # runs all sc in exp_params
run_sequential_phylotime <- function(exp_tb, root_edge = 0.6) {
		exp_tb$tr = map(1:nrow(exp_tb), function(i) {
              message(i)
			        tr = phylotime(exp_tb$sc[[i]],
                       t_total = exp_tb$target_time[i] - root_edge,
                       mut_p = exp_tb$mut_p[[i]],
                       parallel = T)
			})
		exp_tb
}

# start running phylotime inference on subset of all sites
job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
input_dir = "../exp_runs/output/5rep_mut_all/"
data_out = readRDS(paste0(input_dir, stringr::str_pad(job_id, width = 4, pad = "0"), ".rds"))

output_dir = "./panel_mod2_v1/"
if (!dir.exists(output_dir)) {
 dir.create(output_dir)
}
rep_dir = paste0(output_dir, "5rep_mut_all/")
if (!dir.exists(rep_dir)) {
 dir.create(rep_dir)
}
all_mut_indices = readRDS(paste0(output_dir, "all_mut_indices.rds"))
mut_p_all = readRDS(paste0(output_dir, "mut_p_all.rds"))
all_mut_indices$sc = map(all_mut_indices$mut_indices, function(m_indices) {
 data_out$sc[, m_indices]
})
all_mut_indices$mut_p = map(all_mut_indices$mut_indices, function(m_indices) {
 subset_mut_p(mut_p_all, m_indices)
})
all_mut_indices$target_time = 11.5
all_mut_indices = run_sequential_phylotime(all_mut_indices)
out_file = paste0(rep_dir, stringr::str_pad(job_id, width = 4, pad = "0"), ".rds")
saveRDS(all_mut_indices, file = out_file)
# end running phylotime

# collecting and dumping phylotime inferred phylogenies (tr3), flattening the experiments
rep_dir = "./panel_mod2_v1/5rep_mut_all/"
tr3_dir = "./panel_mod2_v1/5rep_mut_all/tr3/"
# dir.create(tr3_dir)
exp_params = readRDS("../exp_runs/output/exp_data_5rep_mut_all.rds")
mut_tr3_all = bind_rows(map(1:nrow(exp_params), function(job_id) {
  message(job_id)
  out_file = paste0(rep_dir, stringr::str_pad(job_id, width = 4, pad = "0"), ".rds")
  mut_tb = readRDS(out_file)
  assertthat::assert_that(nrow(mut_tb) == 12)
  mut_tb$exp_id = job_id
  mut_tb$mut_id = 12*(job_id - 1) + 1:12
  walk(1:12, function(i) {
          tr = mut_tb$tr[[i]]
          out_file = paste0(tr3_dir, stringr::str_pad(mut_tb$mut_id[i], width = 4, pad = "0"), ".rds")
          saveRDS(tr, out_file)
  })
  mut_tb = mut_tb[c(1, 2, 6, 8, 9)]
  mut_tb
  }))
rep_dir = "./panel_mod2_v1/5rep_mut_all/"
saveRDS(mut_tr3_all, file = paste0(rep_dir, "mut_flat_all.rds"))

# mixed 50 hgRNAs below
input_dir = "../exp_runs/output/5rep_mut_all/"
data_out = readRDS(paste0(input_dir, stringr::str_pad(job_id, width = 4, pad = "0"), ".rds"))

output_dir = "./panel_mod2_v1/"
mut_p_all = readRDS(paste0(output_dir, "mut_p_all.rds"))
exp_params = readRDS("../exp_runs/output/exp_data_5rep_mut_all.rds")
tree_panel = readRDS("../exp_runs/output/tree_panel.rds")

rep_dir = paste0(output_dir, "5rep_mut_all/")
# tr3_dir = paste0(rep_dir, "tr3_50site/")
# tr3_dir = paste0(rep_dir, "tr3_100site/")
tr3_dir = paste0(rep_dir, "tr3_25site/")
if (!dir.exists(tr3_dir)) {
        dir.create(tr3_dir)
}
# m_indices = readRDS(paste0(output_dir, "mut_indices50.rds"))
# m_indices = readRDS(paste0(output_dir, "mut_indices100.rds"))
m_indices = readRDS(paste0(output_dir, "mut_indices25.rds"))
mut_p = subset_mut_p(mut_p_all, m_indices)

out_file = paste0(tr3_dir, stringr::str_pad(job_id, width = 4, pad = "0"), ".rds")
if (!file.exists(out_file)) {
        sc = data_out$sc[, m_indices]
        graph_id = exp_params$big_graph_id[job_id]
        tt = tree_panel$type_graph[[graph_id]]$target_time
        tr = phylotime(sc_mat = sc, t_total = tt - 0.6, mut_p = mut_p, parallel = T)
        saveRDS(tr, out_file)
}


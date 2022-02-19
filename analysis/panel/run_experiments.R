library(qfm)
#run_panel
output_dir = "./intermediate_data/panel/"
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))
exp_params = readRDS(paste0(output_dir, "exp_data_10rep.rds"))
rep_dir = paste0(output_dir, "exp_data_10rep/")
# dir.create(rep_dir)
job_id = as.numeric(Sys.getenv("SGE_TASK_ID"))
print(job_id)
set.seed(73)
job_split = split(sample(1:nrow(exp_params), replace = F), ceiling(1:nrow(exp_params)/50))
assertthat::assert_that(all(sort(purrr::reduce(job_split, c)) == 1:nrow(exp_params)))
job_indices = sort(job_split[[job_id]])
set.seed(i)
for (i in job_indices) {
        message(i)
        big_graph_id = exp_params$big_graph_id[i]
        sample_size = exp_params$sample_size[i]
        mut_p = exp_params$mut_p[[i]]
        sampling = exp_params$sampling[i]

        if (sampling == "fixed") {
                ss = sample_size
        }
        if (sampling == "proportional") {
                fm = all_graphs[[big_graph_id]]
                gens0 = make_gens(fm)
                tip_size = purrr::map_dbl(gens0[fm$tip_id], "end_count")
                tip_sample = extraDistr::rmvhyper(nn = 1, k = length(fm$tip_id) * sample_size, n = tip_size)[1, ]
                tip_sample = pmax(tip_sample, 5)
                names(tip_sample) = fm$tip_id
                ss = tip_sample
        }
        out_file = paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        if (!file.exists(out_file)) {
                try({
                        out = simulate_sc_data(all_graphs[[big_graph_id]], mut_p, sample_size = ss)
                        saveRDS(out, file = out_file)
                })
        }
}


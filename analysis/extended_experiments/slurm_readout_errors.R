library(qfm)
library(furrr)
n_cores = 12
plan(multisession, workers = n_cores)

drop_alleles <- function(x, frac = 0.05) {
        if (frac == 0.) {
                return(x)
        }
        x[sample(length(x), size = floor(length(x) * frac), replace = F)] = NA
        x
}
make_existing_error_alleles <- function(x, frac = 0.01) {
        if (frac == 0.) {
                return(x)
        }
        for (j in 1:ncol(x)) {
                y = x[, j]
                y_counts = sort(table(y))
                if (length(y_counts) == 0) {
                        return(x)
                } else {
                        y_uniq = names(y_counts)
                        y_uniq_freq = as.numeric(y_counts/sum(y_counts))
                }
                num_err = floor(length(y) * frac)
                y[sample(length(y), size = num_err, replace = F)] = sample(y_uniq, size = num_err, prob = y_uniq_freq, replace = T)
                x[, j] = y
        }
        return(x)
}
# runs all sc in exp_params
run_sequential_phylotime <- function(exp_tb, root_edge = 0.6) {
        exp_tb$tr = map(1:nrow(exp_tb), function(i) {
                tr = phylotime(exp_tb$sc[[i]],
                               t_total = exp_tb$target_time[i] - root_edge,
                               mut_p = exp_tb$mut_p[[i]],
                               parallel = T)
        })
        exp_tb
}

# start dropout
output_dir = "./panel_mod2_v1/"
exp_params = readRDS(paste0(output_dir, "phylotime_input.rds"))
# exp_params = exp_params[map_dbl(exp_params$sc, nrow) < 2000, ]

job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # 1-106
rep_dir = paste0(output_dir, "dropout/")
if (!dir.exists(rep_dir)) {
        dir.create(rep_dir)
}
sc_mat = exp_params$sc[[job_id]]
colnames(sc_mat) = paste0("el_", 1:ncol(sc_mat))
frac_vec = c(0.05, 0.1, 0.2, 0.5, 0.3, 0.4)
dropout_tb = tibble(frac = frac_vec)
dropout_tb = mutate(dropout_tb, sc = map(frac, function(val) {
        message(val)
        if (val == 0) {
                return(sc_mat)
        }
        out = im_chr(drop_alleles(sc_mat, frac = val))
        out
}))
dropout_tb$mut_p = rep(exp_params$mut_p[job_id], length(frac_vec))
dropout_tb$target_time = exp_params$target_time[job_id]
dropout_tb = run_sequential_phylotime(dropout_tb)

out_file = paste0(rep_dir, stringr::str_pad(job_id, width = 4, pad = "0"), ".rds")
saveRDS(dropout_tb, file = out_file)
# end dropout

# start allele switching
output_dir = "./panel_mod2_v1/"
exp_params = readRDS(paste0(output_dir, "phylotime_input.rds"))
# exp_params = exp_params[map_dbl(exp_params$sc, nrow) < 2000, ]

job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # 1-106
rep_dir = paste0(output_dir, "exist_error_allele_v1/")
if (!dir.exists(rep_dir)) {
        dir.create(rep_dir)
}
sc_mat = exp_params$sc[[job_id]]
colnames(sc_mat) = paste0("el_", 1:ncol(sc_mat))
# frac_vec = c(0.3, 0.4)
frac_vec = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
dropout_tb = tibble(frac = frac_vec)
dropout_tb = mutate(dropout_tb, sc = map(frac, function(val) {
        message(val)
        if (val == 0) {
                return(sc_mat)
        }
        out = make_existing_error_alleles(sc_mat, frac = val)
        out
}))
dropout_tb$mut_p = rep(exp_params$mut_p[job_id], length(frac_vec))
dropout_tb$target_time = exp_params$target_time[job_id]
dropout_tb = run_sequential_phylotime(dropout_tb)

out_file = paste0(rep_dir, stringr::str_pad(job_id, width = 4, pad = "0"), ".rds")
saveRDS(dropout_tb, file = out_file)
# end allele switching

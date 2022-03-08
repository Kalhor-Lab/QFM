# supplementary_data 1: all graphs
output_dir = "./intermediate_data/panel/"
exp_params = readRDS(paste0(output_dir, "exp_data_10rep_proc.rds"))
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))
data_dir = "./supplementary_data/supplementary_data_1_fate_maps/"
# dir.create(data_dir)
# put fate map into a more readable format and print
for (j in sort(unique(exp_params$big_graph_id))) {
        type_graph = all_graphs[[j]][c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13)]
        names(type_graph) = c("root_state",
                              "progenitor_states",
                              "terminal_states",
                              "merge",
                              "commitment_time",
                              "commitment_bias",
                              "doubling time",
                              "founder_size",
                              "target_time",
                              "commitment_side",
                              "edges",
                              "root_state_time")
        type_graph$commitment_time = unlist(type_graph$commitment_time)
        type_graph$commitment_time = type_graph$commitment_time[type_graph$progenitor_states]
        type_graph$`doubling time` = unlist(type_graph$`doubling time`)
        type_graph$commitment_bias = do.call(rbind, map(type_graph$commitment_bias, function(x) x[1:2]))
        temp = do.call(rbind, map(type_graph$commitment_side, function(x) x[1:2, 1]))
        rownames(temp) = names(type_graph$commitment_side)
        type_graph$commitment_side = temp
        sink(paste0(data_dir, str_pad(j, width = 4, pad = "0"), ".txt"))
        print(type_graph)
        sink()
}

# supplementary_data 4: true phylogram, phylotime reconstructed phylogram, lineage barcodes
dump_dir = "./supplementary_data/supplementary_data_4_simulation_data/"
dir.create(dump_dir, recursive = T)
# table for the setup for each experiment
exp_summary_tb = transmute(exp_params,
                           experiment_index = j,
                           fate_map_id = big_graph_id,
                           sample_size,
                           num_sites = num_element,
                           sampling,
                           num_tip = num_tip,
                           bsum = bsum)
write_csv(exp_summary_tb, paste0(dump_dir, "experiment_metadata.csv"))

# using indelphi probabilities to find out the ids used in each experiment
indelphi_tb = readRDS("./intermediate_data/MARC1_indelphi_predicted.rds")
indelphi_tb = indelphi_tb %>% mutate(recur_vec = purrr::map(indelphi, function(x) {
        p = x$prob
        p = sort(p/sum(p), decreasing = T)
        names(p) = paste0("R-", 1:length(p))
        p
}))
exp_id_list = map(exp_params$mut_p, function(x) {
        indelphi_tb$id[match(map_dbl(x$recur_vec_list, 1),
                             map_dbl(indelphi_tb$recur_vec, 1))]
})
exp_hgRNA = bind_cols(tibble(Experiment = 1:nrow(exp_params)),
                      as_tibble(do.call(rbind, exp_id_list)))
write_csv(exp_hgRNA, paste0(dump_dir, "experiment_MARC1_hgRNA.csv"))

dir.create(paste0(dump_dir, "/barcodes/"))
dir.create(paste0(dump_dir, "/trees_simulated/"))
dir.create(paste0(dump_dir, "/trees_phylotime/"))
dump_sc_mat(exp_params$sc,
            output_dir = paste0(dump_dir, "/barcodes/"))
dump_trees(exp_params$tr3, paste0(dump_dir, "/trees_phylotime/"))
dump_trees(exp_params$tr, paste0(dump_dir, "/trees_simulated/"))
# dump the set of MARC1 hgRNAs used for simulation

dump_mut_p_list <- function(mut_p_list, output_dir) {
        for (j in 1:length(mut_p_list)) {
                mutp_dir = paste0(output_dir, "exp_", stringr::str_pad(j, width = 4, pad = "0"), "/")
                dir.create(mutp_dir)
                dump_mut_p(mut_p_list[[j]], mutp_dir)
        }
}

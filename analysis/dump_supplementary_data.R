# supplementary_data 1: all graphs
output_dir = "./intermediate_data/panel/"
exp_params = readRDS(paste0(output_dir, "exp_data_rep1_2_proc2.rds"))
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))
data_dir = "./supplementary_data/supplementary_data_1_fate_maps/"
# dir.create(data_dir)
# put fate map into a more readable format and print
for (j in 1:length(all_graphs)) {
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
# dir.create(dump_dir)
dir.create(paste0(dump_dir, "/mutation_parameters/"))
dir.create(paste0(dump_dir, "/barcodes/"))
dir.create(paste0(dump_dir, "/trees_simulated/"))
dir.create(paste0(dump_dir, "/trees_phylotime/"))
dump_mut_p_list(exp_params$mut_p, output_dir = paste0(dump_dir, "mutation_parameters/"))
dump_sc_mat(exp_params$sc,
            output_dir = paste0(dump_dir, "/barcodes/"))
dump_trees(exp_params$tr3, paste0(dump_dir, "/trees_phylotime/"))
dump_trees(exp_params$tr, paste0(dump_dir, "/trees_simulated/"))




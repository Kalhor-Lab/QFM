library(qfm)
mut_p = readRDS("./data/example/mut_p_marc1.rds")
mut_rate_tb = tibble(site = paste0("site_", 1:length(mut_p$mut_rate)),
                     rate = mut_p$mut_rate)
allele_prob_tb = bind_rows(map(1:length(mut_p$recur_vec_list), function(i_site) {
        allele_prob = mut_p$recur_vec_list[[i_site]]
        tibble(site = paste0("Site_", i_site),
               mutation = paste0("Site", i_site, "_", str_replace_all(names(allele_prob), "R-", "Mut")),
               prob = allele_prob)
}))
write_csv(mut_rate_tb, file = "./mut_p_mut_rate.csv")
write_csv(allele_prob_tb, file = "./mut_p_recur_vec.csv")

tree_panel = readRDS("./data/example/tree_panel.rds")
type_graph = tree_panel$type_graph[[5]]
sim_data = simulate_sc_data_mod2(type_graph, mut_p = mut_p, sample_size = 100)

chr_mat = sim_data$sc
colnames(chr_mat) = paste0("Site_", 1:ncol(chr_mat))
chr_mat[] = str_replace_all(chr_mat[], "R-", "Mut")
for (j in 1:ncol(chr_mat)) {
        zero_ind = which(chr_mat[, j] == "0")
        chr_mat[, j] = paste0("Site", j, "_", chr_mat[, j])
        chr_mat[zero_ind, j] = "0"
}
cell_allele_table = as_tibble(chr_mat)
cell_allele_table = bind_cols(tibble(cell = rownames(sim_data$sc)),
                              cell_allele_table)
write_csv(cell_allele_table, file = "./mutant_allele_matrix.csv")

# generate missing
for (i in 1:50) {
        temp = cell_allele_table[[i+1]]
        temp[sample(length(temp), size = ceiling(length(temp) * 0.05), replace = F)] = NA
        cell_allele_table[[i+1]] = temp
}
write_csv(cell_allele_table, file = "./data/example/mutant_allele_matrix_with_missing.csv")

cell_type_tb = tibble(cell = cell_allele_table$cell,
                      type = get_type_from_id(cell_allele_table$cell))
cell_type_tb$type = str_replace_all(cell_type_tb$type, "-", "Type")
write_csv(cell_type_tb, file = "./data/example/cell_type.csv")

cell_type_tb = read_csv("./data/example/cell_type.csv")

# Make ground truth fate map
root_id = "5"
node_id = c("1", "2", "3", "4", "5")
tip_id = c("-1", "-2", "-3", "-4", "-5", "-6")
merge0 = do.call(rbind, list("1" = c(-1, -2),
                             "2" = c(-3, -4),
                             "3" = c(-5, -6),
                             "4" = c(1, 2),
                             "5" = c(3, 4)))
root_time = 5.
time_itv = c(4, 2, 2 ,1)
type_time0 = list("1" = root_time + cumsum(time_itv)[4],
                  "2" = root_time + cumsum(time_itv)[3],
                  "3" = root_time + cumsum(time_itv)[2],
                  "4" = root_time + cumsum(time_itv)[1],
                  "5" = root_time)
type_time1 = list("1" = root_time + cumsum(time_itv)[4],
                  "2" = root_time + cumsum(time_itv)[3],
                  "3" = root_time + cumsum(time_itv)[1],
                  "4" = root_time + cumsum(time_itv)[2],
                  "5" = root_time)
tip_time = list("-1" = Inf,
                "-2" = Inf,
                "-3" = Inf,
                "-4" = Inf,
                "-5" = Inf,
                "-6" = Inf)
type_time0 = c(type_time0, tip_time)
type_time1 = c(type_time1, tip_time)

# fixed commitment ratio
type_mode_probs0 = list("1" = c(0.5, 0.5, 0, 0, 0),
                        "2" = c(0.5, 0.5, 0, 0, 0),
                        "3" = c(0.5, 0.5, 0, 0, 0),
                        "4" = c(0.5, 0.5, 0, 0, 0),
                        "5" = c(1/3, 1 - 1/3, 0, 0, 0))
type_mode_probs1 = list("1" = c(0.5, 0.5, 0, 0, 0),
                        "2" = c(0.5, 0.5, 0, 0, 0),
                        "3" = c(0.5, 0.5, 0, 0, 0),
                        "4" = c(0.5, 0.5, 0, 0, 0),
                        "5" = c(1/3, 1 - 1/3, 0, 0, 0))

type_doubletime0 = map(c(node_id, tip_id), function(x) 0.833333)
names(type_doubletime0) = c(node_id, tip_id)
type_doubletime0[["5"]] = 0.70
type_doubletime0[["4"]] = 1.512942
type_doubletime0[["3"]] = 1.269057

type_doubletime1 = map(c(node_id, tip_id), function(x) 0.833333)
names(type_doubletime1) = c(node_id, tip_id)
type_doubletime1[["5"]] = 0.95
type_doubletime1[["4"]] = 0.85
type_doubletime1[["3"]] = 0.65

target_time = 19.
well_graph0 = make_type_graph(root_id = root_id,
                              node_id = node_id,
                              tip_id = tip_id,
                              merge_matrix = merge0,
                              differntiation_time = type_time0,
                              differntiation_mode_probs = type_mode_probs0,
                              double_time = type_doubletime0,
                              founder_size = 1,
                              target_time = target_time)
well_graph1 = make_type_graph(root_id = root_id,
                              node_id = node_id,
                              tip_id = tip_id,
                              merge_matrix = merge0,
                              differntiation_time = type_time1,
                              differntiation_mode_probs = type_mode_probs1,
                              double_time = type_doubletime1,
                              founder_size = 1,
                              target_time = target_time)
plot_type_graph(well_graph0) + plot_type_graph(well_graph1)

save(well_graph0, well_graph1, file = "./intermediate_data/iPSC/ground_truth_gr.rds")

label_map = c(paste0("P", 1:5), paste0("T", 1:6))
names(label_map) = as.character(c(1:5, (-1):(-6)))
type_col = c("3" = "#f768a1",
             "-5" = "#dd3497",
             "-6" = "#7a0177",
             "5" = "#c7e9b4",
             "1" = "#80cdc1",
             "-1" = "#35978f",
             "-2" = "#01665e",
             "2" = "#dfc27d",
             "-3" = "#d8b365",
             "-4" = "#8c510a",
             "4" = "#1d91c0")
g1 = (plot_type_graph(well_graph0,
                      node_col_mapper = function(x) type_col[x],
                      node_label_mapper = function(x) label_map[x],
                      show_node_text = F) +
              theme(text = element_text(size = 12)))
g2 = plot_type_graph(well_graph1,
                     node_col_mapper = function(x) type_col[x],
                     node_label_mapper = function(x) label_map[x],
                     show_node_text = F)
(g1 + g2) %>% push_pdf("sim_graph", w = 5, h = 3.5, dir = "./plots/iPSC/")

cell_mat_g1_fil = readRDS("./intermediate_data/iPSC/cell_mat_g1_filtered.rds")
cell_mat_g2_fil = readRDS("./intermediate_data/iPSC/cell_mat_g2_filtered.rds")
all_cell_mat_chr = readRDS("./intermediate_data/iPSC/cell_mat_all.rds")
# reading mutation rate data
mr_tb = readxl::read_excel("./data/iPSC/Lin96_hgRNA_mutation_rates.xlsx")
t_vec = c(4, 8, 11)
mr_tb$mr = map_dbl(1:nrow(mr_tb), function(i) {
        f_vec = as.numeric(mr_tb[i, c("D4", "D8", "D11")]) / 100
        mean(-log(1 - f_vec + 1e-10) / t_vec)
})
g1_mr_vec = mr_tb$mr[match(colnames(cell_mat_g1_fil), mr_tb$hgRNA)]
g1_mr_vec[colnames(cell_mat_g1_fil) == "GCCAAAAGCT"] = estimate_mut_rate(all_cell_mat_chr, t_total = 13.)[colnames(all_cell_mat_chr) == "GCCAAAAGCT"]

g2_mr_vec = mr_tb$mr[match(colnames(cell_mat_g2_fil), mr_tb$hgRNA)]
g2_mr_vec[colnames(cell_mat_g2_fil) == "GCCAAAAGCT"] = estimate_mut_rate(all_cell_mat_chr, t_total = 13.)[colnames(all_cell_mat_chr) == "GCCAAAAGCT"]

# reading indelphi predicted probabilities
Biostrings::reverseComplement(DNAString("ATCTAACCGA"))
idelphi_tb = read_tsv("./supplementary_data/supplementary_data_3_inDelphi_predictions/iPSC_indelphi_result.zip")
idelphi_tb = idelphi_tb %>% nest(data = -ID)
idelphi_tb$ID = map_chr(idelphi_tb$ID, function(x) {
        as.character(Biostrings::reverseComplement(DNAString(x)))
})
idelphi_tb = mutate(idelphi_tb, recur_vec = map(data, function(tb) {
        tb = filter(tb, Generation == "1st_Generation")
        out = sort(tb$Probability, decreasing = T)
        names(out) = paste0("R-", 1:length(out))
        out
}))
g1_recur_vec = idelphi_tb$recur_vec[match(colnames(cell_mat_g1_fil), idelphi_tb$ID)]
g2_recur_vec = idelphi_tb$recur_vec[match(colnames(cell_mat_g2_fil), idelphi_tb$ID)]

mut_p_g1 = list(mut_rate = pmax(g1_mr_vec, 1e-5),
                active_time = list(c(1, 14)),
                recur_prob = rep(1, length(g1_mr_vec)),
                recur_vec_list = g1_recur_vec)
mut_p_g2 = list(mut_rate = pmax(g2_mr_vec, 1e-5),
                active_time = list(c(1, 14)),
                recur_prob = rep(1, length(g2_mr_vec)),
                recur_vec_list = g2_recur_vec)

plan(multisession, workers = 12)
well_tb = as_tibble(expand.grid(group = c("1", "2"),
                                sample_size = c(60, 80, 100, 120, 140),
                                sim_id = 1:50))
j = 0
well_tb = well_tb %>% mutate(data = pmap(., function(group, sample_size, sim_id) {
        j <<- j + 1
        set.seed(j)
        message(paste0("running replicate: ", j, " total: ", nrow(well_tb)))
        tt = proc.time()
        if (group == "1") {
                sc_data0 = simulate_sc_data(well_graph0, mut_p_g1, sample_size)
        }
        if (group == "2") {
                sc_data0 = simulate_sc_data(well_graph1, mut_p_g2, sample_size)
        }
        plan(multisession, workers = 12)
        tr3_g0 = name_nodes(phylotime(sc_data0$sc, t_total = 13.))
        sc_celltypes = get_type_from_id(tr3_g0$tip.label)
        gr3_g0 = name_nodes(reconstruct_graph(tr3_g0, sc_celltypes = sc_celltypes, total_time = 13., theta = 0.0, stat_func = mean))

        data_obj = make_gr_tr_data(gr3_g0, tr3_g0, sc_celltypes)
        tr_node_assign = assign_node_states(data_obj)
        gr_trans_time = est_transition_time(data_obj, tr_node_assign, stat_func = mean)

        data_obj = update_edge_tb_state(data_obj, tr_node_assign)
        gr_node_sizes = get_node_size(data_obj, tr_node_assign, gr_trans_time)
        print('time lapsed: ')
        print(proc.time() - tt)
        return(list(data = sc_data0,
                    gr_tr_data = data_obj,
                    gr3_trans_time = gr_trans_time,
                    gr3_node_size = gr_node_sizes))
}))
saveRDS(well_tb, "./intermediate_data/iPSC/ice_fase_simulation.rds")

# # well_tb = as_tibble(expand.grid(group = c("2"),
# #                                 sample_size = c(60, 80, 100, 120, 140),
# #                                 sim_id = 1:50))
# # j = 0
# # well_tb = well_tb %>% mutate(data = pmap(., function(group, sample_size, sim_id) {
# #         j <<- j + 1
# #         message(paste0("running replicate: ", j, " total: ", nrow(well_tb)))
# #         tt = proc.time()
# #         if (group == "1") {
# #                 sc_data0 = simulate_sc_data(well_graph0, mut_p, sample_size)
# #         }
# #         if (group == "2") {
# #                 sc_data0 = simulate_sc_data(well_graph1, mut_p1, sample_size)
# #         }
# #         plan(multisession, workers = 12)
# #         tr3_g0 = name_nodes(reconstruct_tr1(sc_data0$sc, t_total = 13.))
# #         sc_celltypes = get_type_from_id(tr3_g0$tip.label)
# #         gr3_g0 = name_nodes(reconstruct_graph(tr3_g0, sc_celltypes = sc_celltypes, total_time = 13., theta = 0.5, stat_func = mean))
# #
# #         data_obj = make_gr_tr_data(gr3_g0, tr3_g0, sc_celltypes)
# #         tr_node_assign = assign_node_states(data_obj)
# #         gr_trans_time = est_transition_time(data_obj, tr_node_assign, stat_func = mean)
# #
# #         data_obj = update_edge_tb_state(data_obj, tr_node_assign)
# #         gr_node_sizes = get_node_size(data_obj, tr_node_assign, gr_trans_time)
# #         print('time lapsed: ')
# #         print(proc.time() - tt)
# #         return(list(data = sc_data0,
# #                     gr_tr_data = data_obj,
# #                     gr3_trans_time = gr_trans_time,
# #                     gr3_node_size = gr_node_sizes))
# # }))
# # saveRDS(well_tb, "../LTModelData/splitwell/well_tb_sample_size_g2_29hgRNA.rds")
#
# plot_barcodes(well_tb$data[[1]]$data$sc, well_tb$data[[1]]$data$tr)
#
# well_tb_g1 = readRDS("../LTModelData/splitwell/well_tb_sample_size_v2.rds")
# well_tb_g1$group = paste0(well_tb_g1$group, "_32")
# well_tb_g2 = readRDS("../LTModelData/splitwell/well_tb_sample_size_g2_29hgRNA.rds")
# well_tb_g2$group = paste0(well_tb_g2$group, "_29")
# well_tb = bind_rows(well_tb_g1, well_tb_g2)
#
# # well_tb = readRDS("../LTModelData/splitwell/well_tb_sample_size.rds")
# well_tb$node_size = map(well_tb$data, function(x) {
#         tr3 = x$gr_tr_data$tr
#         gr3 = x$gr_tr_data$gr
#         sc_celltypes = get_type_from_id(tr3$tip.label)
#         data_obj = x$gr_tr_data
#
#         node_mapped = get_node_mapping(data_obj, well_graph0)
#         tip_mapped = well_graph0$tip_id
#         names(tip_mapped) = well_graph0$tip_id
#         node_mapped = c(node_mapped, tip_mapped)
#         out = x$gr3_node_size[names(node_mapped)[match(well_graph0$node_id, node_mapped)]]
#         names(out) = well_graph0$node_id
#         out = map(out, function(x) {
#                 if (is.null(x)) return(NULL)
#                 names(x) = node_mapped[names(x)]
#                 x
#         })
#         out
# })
#
# well_tb$node_time = map(well_tb$data, function(x) {
#         tr3 = x$gr_tr_data$tr
#         gr3 = x$gr_tr_data$gr
#         sc_celltypes = get_type_from_id(tr3$tip.label)
#         data_obj = x$gr_tr_data
#
#         node_mapped = get_node_mapping(data_obj, well_graph0)
#         tip_mapped = well_graph0$tip_id
#         names(tip_mapped) = well_graph0$tip_id
#         node_mapped = c(node_mapped, tip_mapped)
#
#         out = x$gr3_trans_time[names(node_mapped)[match(well_graph0$node_id, node_mapped)]]
#         names(out) = well_graph0$node_id
#         out
# })
#
#
# well_node_size = bind_rows(pmap(well_tb, function(group, node_size, sample_size, ...) {
#         bind_rows(map(names(node_size), function(x) {
#                 if (is.null(node_size[[x]])) return(NULL)
#                 y = node_size[[x]]
#                 tibble(group = group,
#                        node = paste0(x, "_", names(y)),
#                        size = y,
#                        sample_size = sample_size)
#         }))
# }))
#
# well_node_size %>%
#         filter(node != "5_NA") %>%
#         ggboxplot(y = "size", x = "sample_size", color = "group") %>%
#                 facet(facet.by = c("node"), nrow = 2)
#
# # x = well_tb$data[[1]]
# well_tb$kc0 = map_dbl(well_tb$data, function(x) {
#         treespace::treeDist(x$gr_tr_data$gr, as.phylo(well_graph0), lambda = 0)
# })
#
# node_mapper = c("-1_-2_-3_-4_-5_-6" = "5",
#                 "-5_-6" = "3",
#                 "-1_-2" = "1",
#                 "-3_-4" = "2",
#                 "-1_-2_-3_-4" = "4")
#
# out_data = well_tb$data[[1]]
# well_tb$node_order = map_chr(well_tb$data, function(out_data) {
#         kc0 = treespace::treeDist(out_data$gr_tr_data$gr, as.phylo(well_graph0))
#         if (kc0 > 0) {
#                 return(NA)
#         }
#         node_mapped = map_chr(out_data$gr_tr_data$gr_tip_list, function(x) node_mapper[paste0(sort(x, decreasing = F), collapse = "_")])[1:5]
#         node_order = node_mapped[names(sort(out_data$gr3_trans_time[1:5]))]
#         if (all(node_order == c("5", "4", "3", "2", "1"))) {
#                 return("X")
#         }
#         if (all(node_order == c("5", "3", "4", "2", "1"))) {
#                 return("Y")
#         }
#         if (all(node_order == c("5", "4", "3", "1", "2"))) {
#                 return("XA")
#         }
#         if (all(node_order == c("5", "3", "4", "1", "2"))) {
#                 return("YA")
#         }
#         return("N")
# })
# well_tb$correct_topo = !is.na(well_tb$node_order)
#
# node_ord_col = c("X" = "#3D2C8D", "XA" = "#916BBF", "Y" = "#FF7777", "YA" = "#F5C6A5", "N" = "#C8C6C6", "NA" = "#7F7C82")
# well_tb %>% group_by(sample_size, group) %>% dplyr::count(node_order) %>%
#         mutate(node_order = replace_na(node_order, "NA")) %>%
#         mutate(node_order = factor(node_order, levels = c("X", "XA", "Y", "YA", "N", "NA"))) %>%
#         ggbarplot(x = "sample_size", y = "n", fill = "node_order", facet.by = "group") +
#         scale_fill_manual(values = node_ord_col) +
#         theme(legend.position = "right")

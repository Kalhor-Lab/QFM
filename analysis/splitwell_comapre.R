plot_dir = "./plots/iPSC/"

well_tb = readRDS("./intermediate_data/iPSC/ice_fase_simulation.rds")
rs_tb  = readRDS("./intermediate_data/iPSC/ice_fase_resampled.rds")
load("./intermediate_data/iPSC/ground_truth_gr.rds")
# cell_mat_g1 = readRDS("./intermediate_data/iPSC/cell_mat_g1.rds")
# cell_mat_g2 = readRDS("./intermediate_data/iPSC/cell_mat_g2.rds")

set.seed(37)
get_type <- function(x) {
        str_match(x, pattern = "(G\\d)T(\\d)D(.*)")[, 3]
}
get_group <- function(x) {
        str_match(x, pattern = "(G\\d)T(\\d)D(.*)")[, 2]
}
node_mapper = c("1_2_3_4_5_6" = "5",
                "5_6" = "3",
                "1_2" = "1",
                "3_4" = "2",
                "1_2_3_4" = "4")
node_dd = c("5" = "3_4",
            "4" = "1_2",
            "3" = "-5_-6",
            "2" = "-3_-4",
            "1" = "-1_-2",
            "-1" = "-1", "-2" = "-2", "-3" = "-3",
            "-4" = "-4", "-5" = "-5", "-6" = "-6")

ref_gr = as.phylo.type_graph(well_graph0)
ref_gr$tip.label = str_replace(ref_gr$tip.label, "-", "")
rs_tb$node_order = map_chr(rs_tb$data, function(out_data) {
        kc0 = treespace::treeDist(out_data$gr, ref_gr)
        if (kc0 > 0) {
                return(NA)
        }
        node_mapped = map_chr(out_data$gr_tr_data$gr_tip_list, function(x) node_mapper[paste0(sort(x), collapse = "_")])[1:5]
        out_data$gr_trans_time = correct_trans_time(out_data$gr_trans_time, out_data$gr_tr_data) # TODO: to be removed in next version
        node_order = node_mapped[names(sort(out_data$gr_trans_time[1:5]))]
        if (all(node_order == c("5", "4", "3", "2", "1"))) {
                return("X")
        }
        if (all(node_order == c("5", "3", "4", "2", "1"))) {
                return("Y")
        }
        if (all(node_order == c("5", "4", "3", "1", "2"))) {
                return("XA")
        }
        if (all(node_order == c("5", "3", "4", "1", "2"))) {
                return("YA")
        }
        # message('N')
        return("N")
})

ref_gr = as.phylo.type_graph(well_graph0)
well_tb$node_order = map_chr(well_tb$data, function(out_data) {
        kc0 = treespace::treeDist(out_data$gr_tr_data$gr, ref_gr)
        if (kc0 > 0) {
                return(NA)
        }
        node_mapped = map_chr(out_data$gr_tr_data$gr_tip_list, function(x) node_mapper[paste0(sort(stringr::str_replace_all(x, "-", ""), decreasing = F), collapse = "_")])[1:5]
        node_order = node_mapped[names(sort(out_data$gr3_trans_time[1:5]))]
        if (all(node_order == c("5", "4", "3", "2", "1"))) {
                return("X")
        }
        if (all(node_order == c("5", "3", "4", "2", "1"))) {
                return("Y")
        }
        if (all(node_order == c("5", "4", "3", "1", "2"))) {
                return("XA")
        }
        if (all(node_order == c("5", "3", "4", "1", "2"))) {
                return("YA")
        }
        return("N")
})
combine_tb = bind_rows(select(well_tb, group, sample_size, node_order) %>% mutate(type = "Sim"),
                       select(rs_tb, group, sample_size, node_order) %>% mutate(type = "Exp"))

node_ord_col = c("X" = "#3D2C8D", "XA" = "#916BBF", "Y" = "#f768a1", "YA" = "#F5C6A5", "N" = "#C8C6C6", "NA" = "#7F7C82")
(combine_tb %>% group_by(type, group, sample_size) %>% dplyr::count(node_order) %>%
        mutate(node_order = replace_na(node_order, "NA")) %>%
        mutate(node_order = factor(node_order, levels = c("X", "XA", "Y", "YA", "N", "NA"))) %>%
        ggbarplot(x = "sample_size", y = "n",
                  fill = "node_order", facet.by = c("group", "type"),
                  xlab = "Sample size per teriminal state", ylab = "# of replicates") +
        scale_fill_manual(values = node_ord_col) +
        theme(legend.position = "none",
              text = element_text(size = 10),
              strip.background = element_blank())) %>%
        push_pdf(file_name = "compare_topo_order", w = 4.5, h = 6, dir = plot_dir)

rs_tb$node_time = map(rs_tb$data, function(x) {
        gr_tr_data = x$gr_tr_data
        # getting node mapped
        node_mapped = map_chr(x$gr_tr_data$gr_tip_list,
                              function(y) node_mapper[paste0(sort(y), collapse = "_")])[1:5]
        tip_vec = as.character((-1):(-6))
        names(tip_vec) = as.character(1:6)
        node_mapped = c(node_mapped, tip_vec)
        # node_mapped

        # filter by out fate
        gr_dd_mapped = map(x$gr_tr_data$gr_dd, function(z) {
                node_mapped[z]
        })
        gr_dd_mapped = gr_dd_mapped[names(node_mapped)]
        names(gr_dd_mapped) = node_mapped
        outfate_ind = map2_lgl(gr_dd_mapped, node_dd[node_mapped], function(a, b) {
                # print(paste0(sort(a), collapse = "_"))
                paste0(sort(a), collapse = "_") == b
        })
        # outfate_ind
        node_mapped[!outfate_ind] = NA
        # end getting node mapped

        gr_trans_time = x$gr_trans_time
        names(gr_trans_time) = node_mapped[names(gr_trans_time)]
        gr_trans_time[well_graph0$node_id] + 1.
})
well_tb$node_time = map(well_tb$data, function(x) {
        tr3 = x$gr_tr_data$tr
        gr3 = x$gr_tr_data$gr
        sc_celltypes = get_type_from_id(tr3$tip.label)
        data_obj = x$gr_tr_data

        node_mapped = get_node_mapping(data_obj, well_graph0)
        tip_mapped = well_graph0$tip_id
        names(tip_mapped) = well_graph0$tip_id
        node_mapped = c(node_mapped, tip_mapped)

        out = x$gr3_trans_time[names(node_mapped)[match(well_graph0$node_id, node_mapped)]]
        names(out) = well_graph0$node_id
        out + 1.
})
combine_tb = bind_rows(select(well_tb, group, sample_size, node_time) %>% mutate(type = "Sim"),
                       select(rs_tb, group, sample_size, node_time) %>% mutate(type = "Exp"))

well_node_time = bind_rows(pmap(combine_tb, function(group, sample_size, type, node_time, ...) {
                tibble(type = type,
                       group = group,
                       sample_size = sample_size,
                       node = names(node_time),
                       time = node_time
                )
}))
compare_col = c("Exp" = "#6998AB", "Sim" = "#9C0F48")

ground_truth = tibble(node_group = c("P3 (E1)", "P3 (E2)", "P4 (E1)", "P4 (E2)"),
                      time = c(11 - 1.512942, 9 - 0.65, 9 - 1.269057, 11 - 0.85))

(well_node_time %>%
        filter(node %in% c("3", "4")) %>%
        filter(!is.na(node)) %>%
        mutate(type_group = paste0(type, "_", group)) %>%
          mutate(node_group = paste0("P", node, " (E", group, ")")) %>%
        ggboxplot(y = "time",
                  x = "sample_size",
                  fill = "type",
                  size = 0.25,
                  ylab = "Time (days)",
                  xlab = "Sample size per tip") %>%
        facet(facet.by = c("node_group"), nrow = 1) +
        scale_fill_manual(values = compare_col) +
        geom_hline(aes(yintercept = time),
                   size = 1.5,
                   color = "#F1D00A",
                   data = ground_truth) +
                theme(legend.position = "right",
                      text = element_text(size= 10),
                      axis.text.x = element_text(angle = 30))) %>%
        push_pdf(file_name = "node_time_est",
                 w = 5.5,
                 h = 2.5,
                 dir = plot_dir)

rs_tb$node_size = map(rs_tb$data, function(x) {
        gr_tr_data = x$gr_tr_data
        tr_node_assign = assign_node_states(gr_tr_data)
        gr_tr_data = update_edge_tb_state(gr_tr_data, tr_node_assign)
        gr_node_size = get_node_size(gr_tr_data, tr_node_assign, trans_time = x$gr_trans_time)

        # getting node mapped
        node_mapped = map_chr(x$gr_tr_data$gr_tip_list,
                              function(y) node_mapper[paste0(sort(y), collapse = "_")])[1:5]
        tip_vec = as.character((-1):(-6))
        names(tip_vec) = as.character(1:6)
        node_mapped = c(node_mapped, tip_vec)
        # node_mapped

        # filter by out fate
        gr_dd_mapped = map(x$gr_tr_data$gr_dd, function(z) {
                node_mapped[z]
        })
        gr_dd_mapped = gr_dd_mapped[names(node_mapped)]
        names(gr_dd_mapped) = node_mapped
        outfate_ind = map2_lgl(gr_dd_mapped, node_dd[node_mapped], function(a, b) {
                # print(paste0(sort(a), collapse = "_"))
                paste0(sort(a), collapse = "_") == b
        })
        # outfate_ind
        node_mapped[!outfate_ind] = NA
        # end getting node mapped

        names(gr_node_size) = node_mapped[names(gr_node_size)]
        gr_node_size = map(gr_node_size, function(z) {
                names(z) = node_mapped[names(z)]
                z
        })
        gr_node_size[well_graph0$node_id]
})

well_tb$node_size = map(well_tb$data, function(x) {
        tr3 = x$gr_tr_data$tr
        gr3 = x$gr_tr_data$gr
        sc_celltypes = get_type_from_id(tr3$tip.label)
        data_obj = x$gr_tr_data

        node_mapped = get_node_mapping(data_obj, well_graph0)
        tip_mapped = well_graph0$tip_id
        names(tip_mapped) = well_graph0$tip_id
        node_mapped = c(node_mapped, tip_mapped)
        out = x$gr3_node_size[names(node_mapped)[match(well_graph0$node_id, node_mapped)]]
        names(out) = well_graph0$node_id
        out = map(out, function(x) {
                if (is.null(x)) return(NULL)
                names(x) = node_mapped[names(x)]
                x
        })
        out
})

combine_tb = bind_rows(select(well_tb, group, sample_size, node_size) %>% mutate(type = "Sim"),
                       select(rs_tb, group, sample_size, node_size) %>% mutate(type = "Exp"))

well_node_size = bind_rows(pmap(combine_tb, function(group, sample_size, type, node_size, ...) {
        bind_rows(map(names(node_size), function(x) {
                if (is.null(node_size[[x]])) return(NULL)
                y = node_size[[x]]
                tibble(type = type,
                       group = group,
                       sample_size = sample_size,
                       node = paste0(x, "_", names(y)),
                       size = y
                       )
        }))
}))
ground_truth = tibble(node_group = c("P3 (E1)", "P3 (E2)", "P4 (E1)", "P4 (E2)"),
                      size = c(1060/2, 760/2, 500/2, 3250/2))
((well_node_size %>%
        filter(!node %in% c("5_NA", "4_NA")) %>%
        filter(node %in% c("4_4", "3_3")) %>%
        mutate(node_group = paste0("P", map_chr(strsplit(node, "_"), 1),
                                   " (E", group, ")")) %>%
        ggboxplot(y = "size", x = "sample_size", fill = "type", size = 0.25,
                  xlab = "Sample size per tip", ylab = "Progenitor population size") %>%
        facet(facet.by = c("node_group"), nrow = 1)) +
                scale_fill_manual(values = compare_col) +
        geom_hline(aes(yintercept = size),
                   size = 1,
                   color = "#F1D00A",
                   data = ground_truth) +
                theme(legend.position = "right",
                      text = element_text(size= 10),
                      axis.text.x = element_text(angle = 30)) +
        ylim(c(0, 250))) %>%
        push_pdf(file_name = "node_size_est",
                 w = 5.5,
                 h = 2.5,
                 dir = plot_dir)


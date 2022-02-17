plot_dir = "./plots/fate_map_figure/"

set.sed(21)
type_col = c("3" = "#f46d43",
             "-5" = "#d73027",
             "-6" = "#a50026",
             "5" = "#fee08b",
             "1" = "#66bd63",
             "-1" = "#1a9850",
             "-2" = "#006837",
             "2" = "#4393c3",
             "-3" = "#2166ac",
             "-4" = "#053061",
             "4" = "#e6ab02")
root_id = "5"
node_id = c("1", "2", "3", "4", "5")
tip_id = c("-1", "-2", "-3", "-4", "-5", "-6")
merge0 = do.call(rbind, list("1" = c(-1, -2),
                             "2" = c(-3, -4),
                             "3" = c(-5, -6),
                             "4" = c(1, 2),
                             "5" = c(3, 4)))
root_time = 3.65
time_itv = c(3.5, 2, 1 ,1)
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
type_mode_probs0 = list("1" = c(3/4, 1/4, 0, 0, 0),
                        "2" = c(2/3, 1/3, 0, 0, 0),
                        "3" = c(1/4, 3/4, 0, 0, 0),
                        "4" = c(1/3, 2/3, 0, 0, 0),
                        "5" = c(1/3, 1 - 1/3, 0, 0, 0))

type_doubletime0 = map(c(node_id, tip_id), function(x) 0.8)
names(type_doubletime0) = c(node_id, tip_id)
type_doubletime0[["5"]] = 1.2
type_doubletime0[["4"]] = 1.1
type_doubletime0[["3"]] = 1.1
type_doubletime0[["1"]] = 0.7
type_doubletime0[["2"]] = 0.5
type_doubletime0[["-1"]] = 0.65
type_doubletime0[["-2"]] = 0.75
type_doubletime0[["-3"]] = 0.45
type_doubletime0[["-4"]] = 0.55
type_doubletime0[["-5"]] = 0.9
type_doubletime0[["-6"]] = 0.8

target_time = 15.
well_graph0 = make_type_graph(root_id = root_id,
                              node_id = node_id,
                              tip_id = tip_id,
                              merge_matrix = merge0,
                              differntiation_time = type_time0,
                              differntiation_mode_probs = type_mode_probs0,
                              double_time = type_doubletime0,
                              founder_size = 1,
                              target_time = target_time)
plot_type_graph(well_graph0, node_col_mapper = function(x) type_col[x])
(plot_type_graph_clean(well_graph0, node_col_mapper = function(x) type_col[x]) +
        theme(axis.text.y = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              text = element_text(size = 8))) %>%
        push_pdf("f1_graph", width = 1.8, h = 1.4, ps = 8, dir = plot_dir)

mut_p_g1$mut_rate = 10^seq(-1.6, -0.8, length = 31)
mut_p_g1$active_time[[1]] = c(1.3, 14)
data_out = sc_data0 = simulate_sc_data(well_graph0, mut_p_g1, 50)
sc_mat = remove_uniform_id(data_out$sc, abund_thres = 1)
plot_barcodes(sc_mat,
              data_out$tr,
              tip_celltype = get_type_from_id(rownames(data_out$sc)),
              celltype_col = type_col,
              show_row_dend = F) %>%
        push_pdf("f1_barcode", width = 2.5, h = 2.5, ps = 8, dir = plot_dir)

tr = data_out$tr
tr_node_type = get_type_from_id(c(tr$node.label, tr$tip.label))

plot_tr <- function(tr, tr_node_type, node_type_col, root_len, target_time) {
        g_tr = as.igraph(tr)
        # set time
        tr_node_time = node.depth.edgelength(tr) + root_len
        names(tr_node_time) = c(tr$tip.label, tr$node.label)
        V(g_tr)$time = tr_node_time[V(g_tr)$name]
        # V(g_tr)$time[V(g_tr)$name %in% tr$tip.label] = target_time
        V(g_tr)$is_tip = factor(V(g_tr)$name %in% tr$tip.label)
        # set type
        V(g_tr)$type = tr_node_type[V(g_tr)$name]
        g_tr = set.edge.attribute(g_tr,
                                  name = "ending",
                                  index=  unlist(incident_edges(g_tr, tr$tip.label, mode = "in")),
                                  value = T)
        E(g_tr)$ending[is.na(E(g_tr)$ending)] = F
        ggraph(g_tr, layout = "dendrogram", height = time) +
                geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
                geom_node_point(aes(color = type, size = is_tip)) +
                scale_edge_alpha_manual(values = c(0.5, 0.3)) +
                scale_edge_width_manual(values = c(0.2, 0.1)) +
                scale_color_manual(values = node_type_col) +
                scale_size_continuous(range = c(0.4, 0.25)) +
                ylim(c(target_time, 0)) +
                ylab("Time (days)") +
                g_theme

}
(plot_tr(tr, tr_node_type, type_col, root_len = 0.6, target_time = 15) +
                theme(axis.text.y = element_text(size = 8),
                      axis.title.y = element_text(size = 8),
                      text = element_text(size = 8))) %>%
        push_pdf("f1_tr", width = 2.25, h = 1.4, ps = 8, dir = plot_dir)

plot_dir = "./plots/fate_map_figure/"

output_dir = "./intermediate_data/panel/"
exp_params = readRDS(paste0(output_dir, "exp_data_rep1_2_proc2.rds"))
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))

trans_col = c(#"NonTrans" = "#C8C6C6",
        "NonTrans" = NA,
        "Trans" = "#FFB319", "DisjointTrans" = "#22577A")

set.seed(73)
j = 9
type_graph = all_graphs[[exp_params$big_graph_id[j]]]
gr_col = gr_color(type_graph, jitter_amount = 2.)

# set.seed(73)
# sc_data = simulate_sc_data(type_graph, mut_p = exp_params$mut_p[[j]], sample_size = 30)
# tr = sc_data$tr
# gr = reconstruct_graph(tr, get_type_from_id(tr$tip.label), total_time = 15 - 0.6, theta = 0.5, stat_func = mean)
# gr = name_nodes(gr)
# treespace::treeDist(gr, as.phylo(type_graph), lambda = 0)
# root_len = max(node.depth.edgelength(tr)) - max(node.depth.edgelength(gr)) + 0.6

tr = exp_params$data[[j]]$tr
gr = reconstruct_graph(tr, get_type_from_id(tr$tip.label), total_time = 15 - 0.6, theta = 0.5, stat_func = mean)
gr = name_nodes(gr)
tr3 = exp_params$tr3[[j]]
gr3 = name_nodes(exp_params$gr3[[j]])
gr3_eval = exp_params$gr3_eval[[j]]
gr3_eval$kc0
gr_col = gr_color(type_graph, jitter_amount = 2.)
gr3_col = gr_color(gr3, jitter_amount = 2.)
root_len = max(node.depth.edgelength(tr3)) - max(node.depth.edgelength(gr3)) + 0.6
node_label_mapper <- function(x) {
        node_labels = toupper(letters)[1:20]
        names(node_labels) = as.character(-1:-20)
        node_label1 = as.character(1:19)
        names(node_label1) = as.character(1:19)
        node_labels = c(node_labels, node_label1)
        node_labels[x]
}

# plotting cophenetic distances
tr_out = phylotime(exp_params$data[[j]]$sc,
                   mut_p = exp_params$mut_p[[j]],
                   t_total = 15 - 0.6, return_dist = T)
dmat = dist_df2mat(tr_out)
tr3 = exp_params$tr3[[j]]
h = Heatmap(dmat[tr3$tip.label, tr3$tip.label],
        col = make_heatmap_col(max_val = 29, mid_val = 20, min_val = 11, col = rev(brewer.pal(3, "RdPu"))),
        name = "Time (days)",
        cluster_rows = as.dendrogram(tr3),
        cluster_columns = as.dendrogram(tr3),
        show_row_names = F,
        show_row_dend = F,
        show_column_names = F,
        show_column_dend = F)
push_png(h, "phylo_dist_mat", dir = plot_dir, w = 4, h = 3, ps = 10, res = 1200)

plot_type_graph_clean(type_graph,
                      node_col_mapper = function(x) gr_col[x],
                      # node_label_mapper = node_label_mapper,
                    show_node_text = F) %>%
       push_pdf(file_name = "graph_simple", dir = plot_dir, w = 3, h = 3, ps = 8)

(plot_type_graph_dendro(type_graph) + g_theme) %>%
        push_pdf(file_name = "graph_dendro", dir = plot_dir, w = 3, h = 1.5, ps = 8)

plot_barcodes(sc_mat = exp_params$sc[[j]],
                  tr = exp_params$tr3[[j]],
                  tip_celltype = get_type_from_id(rownames(exp_params$sc[[j]])),
                  celltype_col = gr_col,
                  show_row_dend = F) %>%
        push_png(file_name = "barcode_tr3", dir = plot_dir, w = 5, h = 5)

b3 = plot_barcodes(sc_mat = exp_params$sc[[j]],
                  tr = exp_params$tr3[[j]],
                  tip_celltype = get_type_from_id(rownames(exp_params$sc[[j]])),
                  celltype_col = gr_col,
                  show_row_dend = F,
                  show_heatmap_legend = F)
push_png(b, "barcode_tr", dir = plot_dir, w = 4, h = 6)

sc_data0 = simulate_sc_data(well_graph0, mut_p_g1, sample_size)

tr = tr
g_tr = as.igraph(tr)
tr_node_time = node.depth.edgelength(tr) + 0.6
names(tr_node_time) = c(tr$tip.label, tr$node.label)
V(g_tr)$time = tr_node_time[V(g_tr)$name]
V(g_tr)$type = get_type_from_id(V(g_tr)$name)
V(g_tr)$time[V(g_tr)$name %in% tr$tip.label] = 15

g_tr = set.edge.attribute(g_tr,
                   name = "ending",
                   index=  unlist(incident_edges(g_tr, tr$tip.label, mode = "in")),
                   value = T)
E(g_tr)$ending[is.na(E(g_tr)$ending)] = F

V(g_tr)$type_masked = V(g_tr)$type
V(g_tr)$type_masked[V(g_tr)$time < 15] = NA


g1 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
        geom_node_point(aes(color = type_masked), size = 0.75) +
        scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        scale_color_manual(values = gr_col) +
        ylim(c(15, 0)) + ylab("")

push_pdf(g1 + g_theme, "tr3_new_masked", width = 6.5, h = 2, dir = plot_dir)
push_png(g1 + g_theme, "tr3_new_masked", width = 6.5, h = 2, res = 1200, dir = plot_dir)

# g_tr = as.igraph(tr3)
# tr_node_time = node.depth.edgelength(tr3)
# names(tr_node_time) = c(tr3$tip.label, tr3$node.label)
# V(g_tr)$time = tr_node_time[V(g_tr)$name]
#
# g_tr = set.edge.attribute(g_tr,
#                           name = "ending",
#                           index=  unlist(incident_edges(g_tr, tr3$tip.label, mode = "in")),
#                           value = T)
# E(g_tr)$ending[is.na(E(g_tr)$ending)] = F

trans_tb = get_transitions(tr, get_type_from_id(tr$tip.label))
# trans_ind = rep("NonTrans", length(V(g_tr)$name))
# names(trans_ind) = V(g_tr)$name
# trans_ind[trans_tb$in_node] = "Trans"
# trans_ind[trans_tb$in_node[trans_tb$disjoint]] = "DisjointTrans"
# V(g_tr)$trans = trans_ind
# g2 = ggraph(g_tr, layout = "dendrogram", height = time) +
#         geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
#         geom_node_point(aes(color = trans), size = 0.75) +
#         scale_edge_alpha_manual(values = c(0.5, 0.25)) +
#         scale_edge_width_manual(values = c(0.2, 0.05)) +
#         scale_color_manual(values = trans_col) +
#         ylim(c(15, 0)) + ylab("")
# push_png(g2 + g_theme, "tr_new_fase", width = 6.5, h = 2, res = 1200, dir = plot_dir)

# plotting legend only
# plot_legend <- function(g) {
#         legend <- cowplot::get_legend(g)
#         grid.newpage()
#         grid.draw(legend)
# }

fate1 = "-1"
fate2 = "-2"
trans_tb_fase1 = trans_tb %>% mutate(fate_split_ind = map2_lgl(out_type1_vec, out_type2_vec, function(x, y) {
        # x = trans_tb$out_type1_vec[[1]]
        # y = trans_tb$out_type2_vec[[1]]
        !all(c(fate1, fate2) %in% x) & !all(c(fate1, fate2) %in% y) & all(c(fate1, fate2) %in% c(x, y))
}))

fate1 = "-7"
fate2 = "-9"
trans_tb_fase2 = trans_tb %>% mutate(fate_split_ind = map2_lgl(out_type1_vec, out_type2_vec, function(x, y) {
        # x = trans_tb$out_type1_vec[[1]]
        # y = trans_tb$out_type2_vec[[1]]
        !all(c(fate1, fate2) %in% x) & !all(c(fate1, fate2) %in% y) & all(c(fate1, fate2) %in% c(x, y))
}))
intersect(trans_tb_fase1$in_node[trans_tb_fase1$fate_split_ind],
          trans_tb_fase2$in_node[trans_tb_fase2$fate_split_ind])

trans_ind = rep("non-FASE", length(V(g_tr)$name))
names(trans_ind) = V(g_tr)$name
trans_ind[trans_tb_fase1$in_node[trans_tb_fase1$fate_split_ind]] = "FASE_-1_-2"
trans_ind[trans_tb_fase2$in_node[trans_tb_fase2$fate_split_ind]] = "FASE_-7_-9"
trans_ind[V(g_tr)$time == 15] = V(g_tr)$type[V(g_tr)$time == 15]
fase_col = c("FASE_-1_-2" = "#ea00d9",
             "FASE_-7_-9" = "#0abdc6",
             "non-FASE" = NA)
fase_col = c(fase_col, gr_col)

V(g_tr)$trans = trans_ind
V(g_tr)$is_tip = (V(g_tr)$time == 15)
g2 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
        geom_node_point(aes(color = trans, size = is_tip)) +
        # TODO: why is geom_rect not working here?
        # geom_rect(aes(xmin = subtr_xrange[1], xmax = subtr_xrange[2],
        #               ymin = subtr_yrange[1], ymax = subtr_yrange[2]),
        #           color = "red", fill = NA) +
        scale_size_manual(values = c(2.5, 0.5)) +
        scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        scale_color_manual(values = fase_col) +
        ylim(c(15, 0)) + ylab("") + g_theme
# push_png(g2 + g_theme, "tr_FASE", width = 6.5, h = 2, res = 1200, dir = plot_dir)

# subtree with fase annotation
lout = create_layout(g_tr, layout = "dendrogram", height = time)
filter(lout, y < 4 & y > 2 & x > 1500)

gr_tr_data = make_gr_tr_data(gr, tr, get_type_from_id(tr$tip.label))
tr_node_assign = assign_node_states(gr_tr_data)
V(g_tr)$assigned_state = tr_node_assign[V(g_tr)$name]

# subtree_node = "type_13_gen_-1_4"
subtree_node = "type_8_gen_0_16"
subtree_node = "type_12_gen_-1_28"
subtree_all = unique(reduce(map(gr_tr_data$tr_tip_list[[subtree_node]],
                                function(tip) {
                                        pp = shortest_paths(g_tr, subtree_node, tip)
                                        pp$vpath[[1]]$name
                                }), c))
# get ranges of subtree to draw rectangle
subtr_xrange = range(filter(lout, name %in% subtree_all)$x)
subtr_yrange = range(filter(lout, name %in% subtree_all)$y)

g_tr_sub = induced.subgraph(g_tr,
                            vids = subtree_all)
# make observed fate label
V(g_tr_sub)$observed_fate =
        map_chr(c(gr_tr_data$tr_tip_type_list,
                  get_type_from_id(gr_tr_data$tr$tip.label))[V(g_tr_sub)$name], function(x) {
                          str_c("[", str_c(str_replace(sort(x), "-", "T"), collapse = ", "),
                                "]")

                  })

# g_sub = ggraph(g_tr_sub, layout = "dendrogram", height = time) +
#         geom_edge_diagonal(color = "black", width = 0.2) +
#         geom_node_point(aes(color = trans), size = 2) +
#         geom_node_text(aes(label = observed_fate), nudge_y = 0.5, size = 2) +
#         scale_color_manual(values = fase_col) +
#         ylim(c(15, 0)) + ylab("")
#
#
#
# push_pdf(g_sub + g_theme,
#          "tr_FASE_subtree_v2",
#          ps = ,
#          width = 3., h = 3, dir = plot_dir)

V(g_tr_sub)$assigned_state_label = str_replace(V(g_tr_sub)$assigned_state, "Node-", "IP")
V(g_tr_sub)$assigned_state_label = str_replace(V(g_tr_sub)$assigned_state_label, "-", "T")
g_sub1 = ggraph(g_tr_sub, layout = "dendrogram", height = time) +
        geom_edge_diagonal(width = 0.2, color = "black") +
        geom_node_point(aes(color = assigned_state), size = 2) +
        scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        geom_node_text(aes(label = assigned_state_label), nudge_y = 0.5, size = 4) +
        scale_color_manual(values = gr3_col) +
        ylim(c(15, 0)) + ylab("")
push_pdf(g_sub1 + g_theme,
         "tr_state_subtree_v1",
         ps = 12,
         width = 4., h = 3, dir = plot_dir)


# FASE time boxplot
(bind_rows(tibble(FASE = "3",
                 Time = trans_tb_fase1$time[trans_tb_fase1$fate_split_ind]),
          tibble(FASE = "13",
                 Time = trans_tb_fase2$time[trans_tb_fase2$fate_split_ind])) %>%
        ggboxplot(x = "FASE", y = "Time", fill = "FASE", ylim = c(15, 0)) + scale_fill_manual(values = gr_col)) %>%
        push_pdf(file_name = "FASE_time_box", width = 3, height = 4., dir = plot_dir)


gr_node_mapped = get_node_mapping(gr_tr_data, type_graph = type_graph)
set.seed(73)
gr_col_alt = gr_color(type_graph, jitter_amount = 20., neg = T)
names(gr_col_alt) = type_graph$node_id

gr3_col = gr_col_alt[gr_node_mapped]
names(gr3_col) = names(gr_node_mapped)
gr3_col = c(gr3_col, gr_col[type_graph$tip_id])
g3 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
        geom_node_point(aes(color = assigned_state), size = 0.75) +
        scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        scale_color_manual(values = gr3_col) +
        ylim(c(15, 0)) + ylab("")
push_png(g3 + g_theme, "tr_node_assign_alt_col", width = 6.5, h = 2, res = 1200, dir = plot_dir)

node_assign_count = table(get_type_from_id(names(tr_node_assign)),
                          tr_node_assign)
node_assign_prob = as.matrix(node_assign_count / rowSums(node_assign_count))
class(node_assign_prob) = "matrix"

map_numeric_label <- function(lab) {
        lab_num = as.numeric(lab)
        out = character(length(lab))
        out[lab_num < 0] = str_replace(lab[lab_num < 0], "-", "T")
        out[lab_num > 0] = paste0("P", lab[lab_num > 0])
        out
}
map_gr_label <- function(lab) {
        out = str_replace(lab, "Node-", "IP")
        out = str_replace(out, "-", "T")
        out
}

rownames(node_assign_prob) = map_numeric_label(rownames(node_assign_prob))
colnames(node_assign_prob) = map_gr_label(colnames(node_assign_prob))

h_node_assign = Heatmap(t(node_assign_prob[map_numeric_label(c(gr_node_mapped, unlist(gr_tr_data$gr_tips))[lout$name[order(lout$y, decreasing = T)]]),
                                           map_gr_label(lout$name[order(lout$y, decreasing = T)])]),
        cluster_rows = F, cluster_columns = F,
        right_annotation = rowAnnotation(Type = lout$name[order(lout$y, decreasing = T)],
                                           col = list(Type = gr3_col),
                                           show_legend = F),
        top_annotation = HeatmapAnnotation(Type_True = c(gr_node_mapped, unlist(gr_tr_data$gr_tips))[lout$name[order(lout$y, decreasing = T)]],
                                            col = list(Type_True = gr_col),
                                         show_legend = F),
        column_names_side = "top",
        show_heatmap_legend = F,
        show_row_names = T,
        show_column_names = T,
        column_names_rot = 45
        )
push_pdf(h_node_assign,
         file_name = "node_assign_heatmap",
         width = 4.5,
         height = 4.5,
         dir = plot_dir)

bind_rows(tibble(`Assigned state` = "Correct",
                 Count = sum(node_assign_count[row(node_assign_count) == col(node_assign_count)])),
          tibble(`Assigned state` = "Loss potent",
                 Count = sum(node_assign_count[row(node_assign_count) != col(node_assign_count)]))) %>%
        ggpie(x = "Count", label = "Assigned state", colro = NA, fill = "Assigned state") +
        scale_fill_manual(values = c("#FBD148", "#DEDEDE")) + theme(text = element_text(size = 12),
                                                                    legend.position = "none")


gr_tip_time = rep(15., length(gr$tip.label))
names(gr_tip_time) = gr$tip.label
plot_gr(gr, root_len, gr3_col, target_time = 15.) + g_theme
# push_pdf(gg3, "gr3", dir = plot_dir, w = 3, h = 2.5, ps = 8)
gr_cat = c(rep("Node", length(gr$node.label)), rep("Tip", length(gr$tip.label)))
names(gr_cat) = c(gr$node.label, gr$tip.label)
plot_gr_clean(gr, c(gr_trans_time, gr_tip_time), gr3_col, gr_node_cat = gr_cat, target_time = 15.) %>%
        push_pdf("gr_recon_alt_col", dir = plot_dir, w = 3, h = 2.5, ps = 12)

plot_gr_dendro(gr, c(gr_trans_time, gr_tip_time), gr3_col, gr_node_cat = gr_cat, target_time = 15.) %>%
        push_pdf(file_name = "graph_dendro", dir = plot_dir, w = 3, h = 1.5, ps = 8)


total_time = 15 - 0.6
theta = 0.0
dmat_nondisjoint = transition2dist(trans_tb, mean, total_time)
dmat_disjoint = transition2dist(filter(trans_tb, disjoint), mean, total_time)
type_names = rownames(dmat_nondisjoint)
dmat = theta * dmat_disjoint[type_names, type_names] +
        (1 - theta) * dmat_nondisjoint[type_names, type_names]

ig_gr3 = as.igraph(gr)
V(ig_gr3)$Time = gr_trans_time[V(ig_gr3)$name]
V(ig_gr3)$Time[is.na(V(ig_gr3)$Time)] = total_time
lout = create_layout(ig_gr3, layout = "dendrogram")
lout_node = lout %>% filter(!leaf) %>% arrange(desc(x))
lout_leaf = lout %>% filter(leaf) %>% arrange(x)

# this is for plotting dmat


diag(dmat) = NA
library(phylogram)
rownames(dmat) = map_numeric_label(rownames(dmat))
colnames(dmat) = map_numeric_label(colnames(dmat))
h = Heatmap(dmat[map_numeric_label(lout_leaf$name), map_numeric_label(lout_leaf$name)], cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = T, show_heatmap_legend = T, name = "Time",
        column_names_side = "top", column_names_rot = 45,
        col = make_heatmap_col(max_val = 25, mid_val = 20, min_val = 15, col = rev(brewer.pal(3, "RdPu"))))
h
push_pdf(h, "dist_mat", dir = plot_dir, w = 3, h = 3, ps = 8)

trans_assign_tb = get_trans_assign_tb(gr_tr_data, tr_node_assign)
gr_trans_time = est_transition_time(gr_tr_data, tr_node_assign)
b_time = trans_assign_tb %>% mutate(in_type = factor(in_type, levels = lout_node$name)) %>%
        ggboxplot(y = "time", x = "in_type", fill = "in_type", ylab = "Time", xlab = "", width = 0.8, size = 0.25) +
        scale_fill_manual(values = gr3_col, guide = F) + ylim(c(15, 0)) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
b_time
push_pdf(b_time, "gr_trans_estimates_box", h = 3, w = 4, dir = plot_dir)

# give an example of one transition
nodes_trans = trans_assign_tb$in_node[trans_assign_tb$in_type %in% c(names(gr_node_mapped)[gr_node_mapped %in% c("3", "13")])]
V(g_tr)$assigned_state_masked = V(g_tr)$assigned_state
V(g_tr)$assigned_state_masked[!V(g_tr)$name %in% nodes_trans] = "masked"
g4 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
        geom_node_point(aes(color = assigned_state_masked), size = 1.5) +
        scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        scale_color_manual(values = c(gr3_col, "masked" = NA)) +
        ylim(c(15, 0)) + ylab("")
push_png(g4 + g_theme, "tr3_state_3_13_ice", width = 6.5, h = 2, res = 1200, dir = plot_dir)


gr_tr_data = update_edge_tb_state(gr_tr_data, tr_node_assign)
gr_node_sizes = get_node_size(gr_tr_data, tr_node_assign, gr_trans_time)
node_x = names(gr_node_mapped)[gr_node_mapped == "3"]
node_x = names(gr_node_mapped)[gr_node_mapped == "13"]
edge_diff = get_edge_diff(node_x, gr_tr_data, gr_trans_time)
edge_diff = edge_diff %>% left_join(igraph::as_data_frame(g_tr) %>% mutate(indices = 1:n()))
E(g_tr)$commit = "masked"
E(g_tr)[edge_diff$indices[edge_diff$d1_ind]]$commit = gr_tr_data$gr_dd[[node_x]][1]
E(g_tr)[edge_diff$indices[edge_diff$d2_ind]]$commit = gr_tr_data$gr_dd[[node_x]][2]
E(g_tr)[edge_diff$indices[!edge_diff$d1_ind & !edge_diff$d2_ind]]$commit = node_x
E(g_tr)$commit[is.na(E(g_tr)$commit)] = "masked"
E(g_tr)$is_masked = E(g_tr)$commit == "masked"

V(g_tr)$assigned_state_masked = V(g_tr)$assigned_state
V(g_tr)$assigned_state_masked[!V(g_tr)$name %in% c(edge_diff$from, edge_diff$to)] = "masked"

assign_commit <- function(node_y) {
        node_dd = gr_tr_data$gr_dd[[node_x]]
        gr_root_node = gr_tr_data$gr$node.label[1]
        cur = node_y
        state_path = node_y
        while (cur != gr_root_node) {
                cur_new = get_parent_node_from_edges(cur, gr_tr_data$gr_edges_tb)
                state_path = c(state_path, cur_new)
                cur = cur_new
        }
        if (node_dd[1] %in% state_path) {
                return(node_dd[1])
        }
        if (node_dd[2] %in% state_path) {
                return(node_dd[2])
        }
        return(node_x)
}
for (node_y in unique(V(g_tr)$assigned_state_masked)) {
        if (node_y != "masked") {
                V(g_tr)$assigned_state_masked[V(g_tr)$assigned_state_masked == node_y] =
                        assign_commit(node_y)
        }
}
commit_col = c("#bebebe", "#ffffbf", "#d53e4f", "#3288bd")
names(commit_col) = c("masked",
                      node_x,
                      gr_tr_data$gr_dd[[node_x]][1],
                      gr_tr_data$gr_dd[[node_x]][2])

g_size = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(color = commit,
                               alpha = factor(is_masked),
                               width = factor(is_masked))) +
        geom_node_point(aes(color = assigned_state_masked), size = 0.5) +
        scale_edge_alpha_manual(values = c(1., 0.25)) +
        scale_edge_width_manual(values = c(0.5, 0.1)) +
        scale_edge_color_manual(values = commit_col) +
        # scale_color_manual(values = c(gr3_col, "masked" = NA)) +
        scale_color_manual(values = commit_col) +
        ylim(c(15, 0)) + ylab("")

# g_size + g_theme
# plot_legend(g_size)
push_png(g_size + g_theme, "node_size_state3", width = 6.5, h = 2, res = 1200, dir = plot_dir)
push_png(g_size + g_theme, "node_size_state13", width = 6.5, h = 2, res = 1200, dir = plot_dir)

node_size_tb = bind_rows(map(names(gr_node_sizes), function(x) {
        tibble(node = x,
               type = c("Un-committed", rep("Committed", 2)),
               split = names(gr_node_sizes[[x]]),
               count = gr_node_sizes[[x]])
}))
b_size = node_size_tb %>% mutate(node = factor(node, levels = lout_node$name)) %>%
        ggbarplot(x = "node", y = "count", fill = "split", ylab = "Count", xlab = "", col = NA) %>%
        facet(facet.by = "type", ncol = 1) + scale_fill_manual(values = gr3_col, guide = F) +
        theme(text = element_text(size = 12),
              axis.text.x = element_blank(),
              # axis.text.x = element_text(angle = 30),
              axis.ticks.x = element_blank())
push_pdf(b_size, "size_barplot", dir = plot_dir, w = 3, h = 3.5, ps = 12)




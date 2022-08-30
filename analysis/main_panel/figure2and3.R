library(qfm)
source("R_mod/mod2_v2.R")
source("R_mod/ice_fase_mod1.R")
source("./analysis/MARC1/load_MARC1_data.R")
source("./R/utils.R")
source("R_mod/exp_pipeline.R")
source("R_mod/plotting.R")

output_dir = "./intermediate_data/panel_mod2_v1/"
plot_dir = "./plots/fate_map_figure/"
# dir.create(plot_dir)
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))
exp_params = readRDS(paste0(output_dir, "exp_data_1rep_proc_v1.rds"))
# exp_params = readRDS(paste0(output_dir, "exp_data_1rep_proc.rds"))

tree_indices = which(tree_panel$num_tip == 16 & tree_panel$category == "random")
tree_indices
tree_panel$bsum[tree_indices]
tree_indices[which.max(tree_panel$bsum[tree_indices])]
tree_indices

map_gr_label <- function(lab) {
        out = str_replace(lab, "Node-", "iP")
        out = str_replace(out, "-", "T")
        out
}

j = 23
type_graph = tree_panel$type_graph[[j]]
set.seed(77)
gr_col = gr_color_v1(tree_panel$type_graph[[j]], jitter_amount = 3., neg = T)
man_col_map = c("#c7ddad" = "#9ec162",
                "#d7e2d5" = "#19e2f7",
                "#eedeac" = "#f7c852",
                "#f3dbc0" = "#ce5fac",
                "#ceb3b5" = "#84ce89")
gr_col[match(toupper(names(man_col_map)), gr_col)] = man_col_map

plot_indices = rbind(c(1, 54, 162),
                     c(23, 107, 238),
                     c(36, 156, 323))

all_g = map(c(t(plot_indices)), function(i) {
        if (i == 23) {
                gr_col = gr_col
        } else {
                gr_col = gr_color_v1(tree_panel$type_graph[[i]], jitter_amount = 3., neg = T)
        }
        set.seed(77)
        plot_type_graph_clean_mod2(tree_panel$type_graph[[i]],
                                   node_col_mapper = function(x) {
                                           gr_col[x]
                                   })
})
(purrr::reduce(all_g, `+`) + plot_layout(nrow = 3, width = c(1, 2, 4))) %>%
        push_pdf(file_name = "figure1_panel_v1c3", dir = plot_dir,
                 width = 10.5, height = 6.,
                 ps = 12)

g0 = plot_type_graph_clean_mod2(type_graph,
                                node_label_mapper = node_label_mapper,
                                node_col_mapper = function(x) gr_col[x],
                                show_node_text = T)

# set.seed(88)
# data_out = simulate_sc_data_mod2(type_graph, mut_p = exp_params$mut_p[[1]], sample_size = 100)
# data_out = simulate_sc_data_mod2(type_graph, mut_p = exp_params$mut_p[[1]], sample_size = 50)
tr = exp_params$tr[[j]]
g_tr = as.igraph(tr)
tr_node_time = node.depth.edgelength(tr) + 0.6
names(tr_node_time) = c(tr$tip.label, tr$node.label)
V(g_tr)$time = tr_node_time[V(g_tr)$name]
V(g_tr)$type = get_type_from_id(V(g_tr)$name)
# V(g_tr)$time[V(g_tr)$name %in% tr$tip.label] = type_graph$target_time
V(g_tr)$time[V(g_tr)$name %in% tr$tip.label] = jitter(rep(type_graph$target_time, length(tr$tip.label)), amount = 0.5)

g_tr = set.edge.attribute(g_tr,
                          name = "ending",
                          index=  unlist(incident_edges(g_tr, tr$tip.label, mode = "in")),
                          value = T)
E(g_tr)$ending[is.na(E(g_tr)$ending)] = F

g1 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(color = ending, width = ending)) +
        geom_node_point(aes(color = type), size = 2.4) +
        # scale_size_manual(values = c(2., 0.05)) +
        # scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_color_manual(values = c("#000000", "8e8e8e")) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        scale_color_manual(values = gr_col) +
        ylim(c(type_graph$target_time, 0)) + ylab("")
push_pdf(g1 + g_theme, "tr_jitter", width = 14, h = 4.4, dir = plot_dir)

# FASE example
# T14-T15
# T8-T9
trans_tb = get_transitions(tr, get_type_from_id(tr$tip.label))
# plotting legend only
# plot_legend <- function(g) {
#         legend <- cowplot::get_legend(g)
#         grid.newpage()
#         grid.draw(legend)
# }

fate1 = "-4"
fate2 = "-5"
trans_tb_fase1 = trans_tb %>% mutate(fate_split_ind = map2_lgl(out_type1_vec, out_type2_vec, function(x, y) {
        # x = trans_tb$out_type1_vec[[1]]
        # y = trans_tb$out_type2_vec[[1]]
        !all(c(fate1, fate2) %in% x) & !all(c(fate1, fate2) %in% y) & all(c(fate1, fate2) %in% c(x, y))
}))

fate3 = "-15"
fate4 = "-16"
trans_tb_fase2 = trans_tb %>% mutate(fate_split_ind = map2_lgl(out_type1_vec, out_type2_vec, function(x, y) {
        # x = trans_tb$out_type1_vec[[1]]
        # y = trans_tb$out_type2_vec[[1]]
        !all(c(fate3, fate4) %in% x) & !all(c(fate3, fate4) %in% y) & all(c(fate3, fate4) %in% c(x, y))
}))
intersect(trans_tb_fase1$in_node[trans_tb_fase1$fate_split_ind],
          trans_tb_fase2$in_node[trans_tb_fase2$fate_split_ind])

fase1_name = paste0("FASE_", fate1, "_", fate2)
fase2_name = paste0("FASE_", fate3, "_", fate4)

trans_ind = rep("non-FASE", length(V(g_tr)$name))
names(trans_ind) = V(g_tr)$name
trans_ind[trans_tb_fase1$in_node[trans_tb_fase1$fate_split_ind]] = fase1_name
trans_ind[trans_tb_fase2$in_node[trans_tb_fase2$fate_split_ind]] = fase2_name
trans_ind[V(g_tr)$time == type_graph$target_time] = V(g_tr)$type[V(g_tr)$time == type_graph$target_time]

unmask_ind = which(V(g_tr)$name %in% tr$tip.label & V(g_tr)$type %in% c(fate1, fate2, fate3, fate4))
trans_ind[unmask_ind] = V(g_tr)$type[unmask_ind]

fase_col = c("#ea00d9",
             "#0abdc6",
             NA)
names(fase_col) = c(fase1_name,
                    fase2_name,
                    "non-FASE")
fase_col = c(fase_col, gr_col)

V(g_tr)$trans = trans_ind
V(g_tr)$is_tip = V(g_tr)$name %in% tr$tip.label
g2 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
        geom_node_point(aes(color = trans, size = is_tip)) +
        scale_size_manual(values = c(1.8, 0.5)) +
        scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        scale_color_manual(values = fase_col) +
        ylim(c(type_graph$target_time, 0)) + ylab("") + g_theme
push_pdf(g2 + g_theme, "tr_fase", width = 5.2, h = 1.5, dir = plot_dir)
# push_png(g2 + g_theme, "tr_fase", width = 6.5, h = 1.5, res = 1200, dir = plot_dir)

(bind_rows(tibble(FASE = fase1_name,
                  Time = trans_tb_fase1$time[trans_tb_fase1$fate_split_ind]),
           tibble(FASE = fase2_name,
                  Time = trans_tb_fase2$time[trans_tb_fase2$fate_split_ind])) %>%
                ggboxplot(x = "FASE", y = "Time", fill = "FASE", ylim = c(type_graph$target_time, 0)) + scale_fill_manual(values = fase_col)) %>%
        push_pdf(file_name = "FASE_time_box", width = 3, height = 4., dir = plot_dir)

res = exp_params$res[[23]]
gr = res$gr
gr_tr_data = res$gr_tr_data
tr_node_assign = res$tr_node_assign
gr_trans_time = res$gr_trans_time
V(g_tr)$assigned_state = tr_node_assign[V(g_tr)$name]

gr_node_mapped = get_node_mapping_mod2(gr_tr_data, type_graph = type_graph)
set.seed(73)
# gr_col_alt = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(15)
gr_col_alt = c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2")[c(1:5, 7)])
names(gr_col_alt) = gr$node.label
gr_col_alt = c(gr_col_alt, gr_col[type_graph$tip_id])
gr3_col = gr_col_alt

g3 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
        geom_node_point(aes(color = assigned_state), size = 0.75) +
        scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        scale_color_manual(values = gr_col_alt) +
        ylim(c(type_graph$target_time, 0)) + ylab("")
push_pdf(g3 + g_theme, "tr_node_assign_alt_col", width = 7.5, h = 1.5, dir = plot_dir)
# push_png(g3 + g_theme, "tr_node_assign_alt_col", width = 6.5, h = 1.5, res = 1200, dir = plot_dir)

(tibble(name = names(gr_col_alt),
       type = grepl("Node", names(gr_col_alt)),
       x = length(gr_col_alt):1,
       y = 1,
       size = 2) %>%
        ggscatter(x = "y", y = "x", color = "name", shape = "type") +
        scale_color_manual(values = gr_col_alt) +
        scale_shape_manual(values = c(17, 18)) +
        geom_text(aes(x = y + 0.25, label = map_gr_label(name))) +
        theme(legend.position = "none", text = element_text(size = 12))) %>%
        push_pdf(file_name = "ip_legend", width = 4, h = 4, dir = plot_dir)


total_time = type_graph$target_time - 0.6
theta = 0.0
dmat_nondisjoint = transition2dist(trans_tb, mean, total_time)
type_names = rownames(dmat_nondisjoint)
dmat = dmat_nondisjoint
dmat = dmat / 2.

ig_gr3 = as.igraph(gr)
V(ig_gr3)$Time = gr_trans_time[V(ig_gr3)$name]
V(ig_gr3)$Time[is.na(V(ig_gr3)$Time)] = total_time
lout = create_layout(ig_gr3, layout = "dendrogram")
lout_node = lout %>% filter(!leaf) %>% arrange(desc(x))
lout_leaf = lout %>% filter(leaf) %>% arrange(x)

# this is for plotting dmat
diag(dmat) = NA
library(phylogram)
rownames(dmat) = node_label_mapper(rownames(dmat))
colnames(dmat) = node_label_mapper(colnames(dmat))
h = Heatmap(dmat[node_label_mapper(lout_leaf$name), node_label_mapper(lout_leaf$name)], cluster_rows = F, cluster_columns = F,
            show_row_names = T, show_column_names = T, show_heatmap_legend = T, name = "Time",
            column_names_side = "bottom", column_names_rot = 90,
            col = make_heatmap_col(max_val = 23, mid_val = 17, min_val = 8.4, col = rev(RColorBrewer::brewer.pal(3, "RdPu"))))
h
push_pdf(h, "dist_mat", dir = plot_dir, w = 3, h = 3, ps = 8)

gr_tip_time = rep(type_graph$target_time, length(gr$tip.label))
names(gr_tip_time) = gr$tip.label
plot_gr_dendro(gr,
               c(gr_trans_time, gr_tip_time),
               gr_col,
               plot_node_point = F,
               type_graph$target_time) %>%
        push_pdf(file_name = "graph_dendro", dir = plot_dir, w = 3, h = 1.5, ps = 8)

plot_gr_clean(res$gr,
              gr_node_time = c(res$gr_trans_time, gr_tip_time),
              total_time = type_graph$target_time,
              type_col = gr_col_alt) %>%
        push_pdf("gr_recon", dir = plot_dir, w = 3, h = 2.5, ps = 12)

# node assignment
node_assign_count = table(get_type_from_id(names(tr_node_assign)),
                          tr_node_assign)
node_assign_prob = as.matrix(node_assign_count / rowSums(node_assign_count))
class(node_assign_prob) = "matrix"

rownames(node_assign_prob) = node_label_mapper(rownames(node_assign_prob))
colnames(node_assign_prob) = map_gr_label(colnames(node_assign_prob))

h_node_assign = Heatmap(t(node_assign_prob[node_label_mapper(c(gr_node_mapped, unlist(gr_tr_data$gr_tips))[lout$name[order(lout$y, decreasing = T)]]),
                                           map_gr_label(lout$name[order(lout$y, decreasing = T)])]),
                        cluster_rows = F, cluster_columns = F,
                        right_annotation = rowAnnotation(Type = lout$name[order(lout$y, decreasing = T)],
                                                         col = list(Type = gr3_col),
                                                         show_legend = F),
                        top_annotation = HeatmapAnnotation(Type_True = c(gr_node_mapped, unlist(gr_tr_data$gr_tips))[lout$name[order(lout$y, decreasing = T)]],
                                                           col = list(Type_True = gr_col),
                                                           show_legend = F),
                        name = "Fraction",
                        column_names_side = "top",
                        show_heatmap_legend = T,
                        show_row_names = T,
                        show_column_names = T,
                        column_names_rot = 90
)
h_node_assign %>% push_pdf(file_name = "node_assign_heatmap",
                           width = 4.5,
                           height = 4.5,
                           dir = plot_dir)


# ICE figure
trans_assign_tb = get_trans_assign_tb(gr_tr_data, tr_node_assign)
b_time = trans_assign_tb %>% mutate(in_type_text = factor(map_gr_label(in_type),
                                                          levels = rev(map_gr_label(lout_node$name))),
                                    in_type = factor(in_type,
                                                     levels = lout_node$name)) %>%
        ggboxplot(y = "time", x = "in_type_text", fill = "in_type", ylab = "Time", xlab = "", width = 0.8, size = 0.25) +
        scale_fill_manual(values = gr3_col, guide = F) + ylim(c(type_graph$target_time, 0)) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_text(angle = 90),
              axis.ticks.x=element_blank(),
              text = element_text(size = 8))
b_time
push_pdf(b_time, "gr_trans_estimates_box", h = 2.5, w = 3, dir = plot_dir)

ice_states = c("Node-8", "Node-5")
nodes_trans = trans_assign_tb$in_node[trans_assign_tb$in_type %in% ice_states]
node_ice_state = trans_assign_tb$in_type[trans_assign_tb$in_type %in% ice_states]
V(g_tr)$assigned_state_masked = V(g_tr)$assigned_state
V(g_tr)$assigned_state_masked[!V(g_tr)$assigned_state %in% ice_states] = "masked"
V(g_tr)$assigned_state_masked[match(nodes_trans, V(g_tr)$name)] = paste0("ICE_", node_ice_state)
gr3_ice_col = rlang::duplicate(gr3_col)
gr3_col = c(gr3_col, "masked" = NA)
gr3_ice_col_vec = gr3_col[ice_states]
names(gr3_ice_col_vec) = paste0("ICE_", names(gr3_ice_col_vec))
gr3_ice_col = c(gr3_ice_col, gr3_ice_col_vec)
gr3_ice_col[ice_states] = c("#f2b8d9", "#e2b288")
gr3_ice_col

g4 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
        geom_node_point(aes(color = assigned_state_masked), size = 1.5) +
        scale_edge_alpha_manual(values = c(0.5, 0.25), guide = F) +
        scale_edge_width_manual(values = c(0.2, 0.05), guide = F) +
        scale_color_manual(values = gr3_ice_col) +
        ylim(c(type_graph$target_time, 0)) + ylab("")
push_pdf(g4 + g_theme, "tr_ice_example", width = 7.5, h = 1.5, dir = plot_dir)
# push_png(g4 + g_theme, "tr_ice_example", width = 6.5, h = 1.5, res = 1200, dir = plot_dir)

node_assign_accuracy = bind_rows(map(1:nrow(exp_params), function(j) {
        res = exp_params$res[[j]]
        gr_node_mapped = get_node_mapping_mod2(res$gr_tr_data, type_graph = tree_panel$type_graph[[exp_params$big_graph_id[j]]])
        tibble(exp_id = j,
               accuracy = mean(get_type_from_id(names(res$tr_node_assign)) == gr_node_mapped[res$tr_node_assign], na.rm = T))
}))
node_assign_accuracy %>% gghistogram(x = "accuracy", fill = "808080", xlab = "Accuracy", ylab = "Count") %>%
        push_pdf(file_name = "node_assign_hist", w = 3, h = 2, dir = plot_dir)


# node sizes figure
E(g_tr)$commit = "masked"
V(g_tr)$assigned_state_masked = V(g_tr)$assigned_state
assign_commit <- function(node_y, node_x) {
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

# get summary of edge state pairs
plot_edge_counts <- function(edge_diff) {
        edge_pair_count = count(edge_diff, from_type, to_type, commit)
        edge_pair_count = arrange(edge_pair_count, commit, desc(n))
        edge_pair_count$pair_id = 1:nrow(edge_pair_count)
        edge_pair_count = edge_pair_count %>%
                mutate(from_y = 1., to_y = 2.)
        edge_pair_plot = gather(edge_pair_count[c(1, 2, 5)], key = "key", value = "state", -pair_id)
        edge_pair_plot$key = 3 - as.numeric(factor(edge_pair_plot$key))
        edge_pair_plot %>% ggscatter(x = "pair_id", y = "key", color = "state", size = 6) +
                scale_color_manual(values = gr3_col) +
                ggnewscale::new_scale_color() +
                geom_segment(aes(x = pair_id, xend = pair_id, y = from_y, yend = to_y, color = commit),
                             size = 2.,
                             data = edge_pair_count) +
                scale_color_manual(values = commit_col) +
                geom_text(aes(label = str_c("x", n),
                              x = pair_id + 0.5, y = 1.5),
                          size = 3,
                          data = edge_pair_count)
}

node_x = "Node-8"
edge_diff = get_edge_diff(node_x, gr_tr_data, gr_trans_time - 0.6)
edge_diff = edge_diff %>% left_join(igraph::as_data_frame(g_tr) %>% mutate(indices = 1:n()))
E(g_tr)[edge_diff$indices[edge_diff$d1_ind]]$commit = gr_tr_data$gr_dd[[node_x]][1]
E(g_tr)[edge_diff$indices[edge_diff$d2_ind]]$commit = gr_tr_data$gr_dd[[node_x]][2]
E(g_tr)[edge_diff$indices[!edge_diff$d1_ind & !edge_diff$d2_ind]]$commit = node_x
V(g_tr)$assigned_state_masked[!V(g_tr)$name %in% c(edge_diff$from, edge_diff$to)] = "masked"
edge_diff$commit[edge_diff$d1_ind] = gr_tr_data$gr_dd[[node_x]][1]
edge_diff$commit[edge_diff$d2_ind] = gr_tr_data$gr_dd[[node_x]][2]
edge_diff$commit[!edge_diff$d1_ind & !edge_diff$d2_ind] = node_x

(plot_edge_counts(edge_diff) + ylim(c(0.9, 2.1)) + theme(legend.position = "none")) %>%
        push_pdf("edge_pair_count_ip8", w = 4, h = 1.5, dir = plot_dir)


# for (node_y in unique(V(g_tr)$assigned_state_masked)) {
#         if (node_y != "masked") {
#                 V(g_tr)$assigned_state_masked[V(g_tr)$assigned_state_masked == node_y] =
#                         assign_commit(node_y, node_x)
#         }
# }
node_x = "Node-5"
edge_diff = get_edge_diff(node_x, gr_tr_data, gr_trans_time - 0.6)
edge_diff = edge_diff %>% left_join(igraph::as_data_frame(g_tr) %>% mutate(indices = 1:n()))
E(g_tr)[edge_diff$indices[edge_diff$d1_ind]]$commit = gr_tr_data$gr_dd[[node_x]][1]
E(g_tr)[edge_diff$indices[edge_diff$d2_ind]]$commit = gr_tr_data$gr_dd[[node_x]][2]
E(g_tr)[edge_diff$indices[!edge_diff$d1_ind & !edge_diff$d2_ind]]$commit = node_x
edge_diff$commit[edge_diff$d1_ind] = gr_tr_data$gr_dd[[node_x]][1]
edge_diff$commit[edge_diff$d2_ind] = gr_tr_data$gr_dd[[node_x]][2]
edge_diff$commit[!edge_diff$d1_ind & !edge_diff$d2_ind] = node_x

(plot_edge_counts(edge_diff) + ylim(c(0.9, 2.1)) + theme(legend.position = "none")) %>%
        push_pdf("edge_pair_count_ip5", w = 2.5, h = 1.5, dir = plot_dir)

# unmask "Node-5" related states
node_state_ind = V(g_tr)$name %in% c(edge_diff$from, edge_diff$to)
V(g_tr)$assigned_state_masked[node_state_ind] = V(g_tr)$assigned_state[node_state_ind]
# for (node_y in unique(V(g_tr)$assigned_state_masked)) {
#         if (node_y != "masked") {
#                 V(g_tr)$assigned_state_masked[V(g_tr)$assigned_state_masked == node_y] =
#                         assign_commit(node_y, node_x)
#         }
# }

V(g_tr)$is_masked = V(g_tr)$assigned_state_masked == "masked"
E(g_tr)$commit[is.na(E(g_tr)$commit)] = "masked"
E(g_tr)$is_masked = E(g_tr)$commit == "masked"
gr3_col = c(gr3_col, "masked" = "#d3d3d3")

# for edges
# commit_col = c("#d3d3d3",
#                "#bebebe",
#                "#ffffbf",
#                "#d53e4f",
#                "#3288bd")
# names(commit_col) = c("masked",
#                       node_x,
#                       gr_tr_data$gr_dd[[node_x]][1],
#                       gr_tr_data$gr_dd[[node_x]][2])

node_x = "Node-8"
commit_col1 = c("#a6611a", "#b2abd2", "#5e3c99")
names(commit_col1) = c(node_x,
                       gr_tr_data$gr_dd[[node_x]][1],
                       gr_tr_data$gr_dd[[node_x]][2])
node_x = "Node-5"
commit_col2 = c("#018571", "#fdb863", "#e66101")
names(commit_col2) = c(node_x,
                       gr_tr_data$gr_dd[[node_x]][1],
                       gr_tr_data$gr_dd[[node_x]][2])
commit_col = c("masked" = "#7e7e7e", commit_col1, commit_col2)

g_size = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_node_point(aes(color = assigned_state_masked,
                            size = is_masked,
                            alpha = is_masked)) +
        geom_edge_diagonal(aes(color = commit,
                               alpha = factor(is_masked),
                               width = factor(is_masked))) +
        scale_size_manual(values = c(0.5, 0.1)) +
        scale_alpha_manual(values = c(1., 0.3)) +
        scale_edge_alpha_manual(values = c(1., 0.2)) +
        scale_edge_width_manual(values = c(0.3, 0.1)) +
        scale_edge_color_manual(values = commit_col) +
        # scale_color_manual(values = c(gr3_col, "masked" = NA)) +
        scale_color_manual(values = gr3_col) +
        ylim(c(type_graph$target_time, 0)) + ylab("")
push_pdf(g_size + g_theme, "node_size_example", width = 6.2, h = 1.8, dir = plot_dir)

# g_size + g_theme
# plot_legend(g_size)
# push_png(g_size + g_theme, "node_size_state3", width = 6.5, h = 2, res = 1200, dir = plot_dir)
# push_png(g_size + g_theme, "node_size_state13", width = 6.5, h = 2, res = 1200, dir = plot_dir)

gr_col_alt_mask = rep("#d6d6d6", length(gr_col_alt))
names(gr_col_alt_mask) = names(gr_col_alt)
node_x = "Node-5"
node_vec = c(node_x, gr_tr_data$gr_dd[[node_x]])
gr_col_alt_mask[node_vec] = gr_col_alt[node_vec]
node_x = "Node-8"
node_vec = c(node_x, gr_tr_data$gr_dd[[node_x]])
gr_col_alt_mask[node_vec] = gr_col_alt[node_vec]
# also color downstream states
for (node_x in c("Node-6", "Node-11", "Node-12", "Node-15", "Node-14")) {
        node_vec = c(node_x, gr_tr_data$gr_dd[[node_x]])
        gr_col_alt_mask[node_vec] = gr_col_alt[node_vec]
}
# edge colors
# adding downstreams
node_x = "Node-5"
commit_col2 = c("#fdb863", rep("#e66101", 10))
names(commit_col2) = c(gr_tr_data$gr_dd[[node_x]][1],
                       gr_tr_data$gr_dd[[node_x]][2],
                       "Node-6", "Node-11", "Node-12", "Node-15",
                       "-11", "-12", "-13", "-14", "-15")
node_x = "Node-8"
commit_col1 = c("#b2abd2", rep("#5e3c99", 3))
names(commit_col1) = c(gr_tr_data$gr_dd[[node_x]][1],
                       gr_tr_data$gr_dd[[node_x]][2],
                       "-3", "-4")
commit_col = c("masked" = "#7e7e7e", commit_col1, commit_col2)

gr_col_alt_commit = rlang::duplicate(gr3_col)
gr_col_alt_commit[!names(gr_col_alt_commit) %in% names(commit_col)] = commit_col["masked"]
gr_col_alt_commit[names(commit_col)] = commit_col

plot_gr_clean(res$gr,
              gr_node_time = c(res$gr_trans_time, gr_tip_time),
              total_time = type_graph$target_time,
              type_col = gr_col_alt_mask,
              edge_col = gr_col_alt_commit,
              edge_width = 0.5) %>%
        push_pdf(file_name = "gr_commit_edge", dir = plot_dir,
                 width = 3., height = 3.,
                 ps = 12)

gr_node_sizes = res$gr_node_sizes
node_size_tb = bind_rows(map(names(gr_node_sizes), function(x) {
        tibble(node = x,
               type = c("Prior to commitment", rep("Immediately post commitment", 2)),
               split = names(gr_node_sizes[[x]]),
               count = gr_node_sizes[[x]])
}))

node_size_tb$node_text = map_gr_label(node_size_tb$node)
node_size_tb$node_text = factor(node_size_tb$node_text,
                                 levels = node_size_tb$node_text[match(lout_node$name, node_size_tb$node)])
node_size_tb$split_text = paste0(map_chr(gr_tr_data$gr_dd[node_size_tb$node], function(x) {
        stringr::str_pad(paste0(map_gr_label(x), collapse = " | "), width = 11, " ", side = "both")
}))
node_size_tb$split_text = factor(node_size_tb$split_text,
                                 levels = node_size_tb$split_text[match(lout_node$name, node_size_tb$node)])

node_size_tb$split = factor(node_size_tb$split, levels = c("Node-1",
                                                           rev(unlist(gr_tr_data$gr_dd[gr_tr_data$gr$node.label]))))

b_size1 = node_size_tb %>%
        filter(type == "Prior to commitment") %>%
        ggbarplot(x = "node_text", y = "count", fill = "split", ylab = "Count", xlab = "", col = NA) +
        # facet(facet.by = "type", ncol = 1) +
        scale_fill_manual(values = gr3_col, guide = F)
b_size2 = node_size_tb %>%
        filter(type == "Immediately post commitment") %>%
        ggbarplot(x = "split_text", y = "count", fill = "split", ylab = "Count", xlab = "", col = NA,
                  position = position_stack()) +
        # facet(facet.by = "type", ncol = 1) +
        scale_fill_manual(values = gr3_col, guide = F)

((b_size1 / b_size2) &
        theme(text = element_text(size = 10),
              # axis.text.x = element_blank(),
              axis.text.x = element_text(angle = 90),
              axis.ticks.x = element_blank())) %>%
        push_pdf("size_barplot", dir = plot_dir, w = 3, h = 4, ps = 10)
# end gr_node_sizes

# generating subtree
lout = create_layout(g_tr, layout = "dendrogram", height = time)
filter(lout, y < 4 & y > 2 & x > 1500)

gr_tr_data = make_gr_tr_data(gr, tr, get_type_from_id(tr$tip.label))
tr_node_assign = assign_node_states(gr_tr_data)
V(g_tr)$assigned_state = tr_node_assign[V(g_tr)$name]

# subtree_node = "type_13_gen_-1_4"
subtree_node = "type_8_gen_0_16"
subtree_node = "type_13_gen_-1_19"
subtree_all = unique(purrr::reduce(map(gr_tr_data$tr_tip_list[[subtree_node]],
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

V(g_tr_sub)$assigned_state_label = str_replace(V(g_tr_sub)$assigned_state, "Node-", "IP")
V(g_tr_sub)$assigned_state_label = str_replace(V(g_tr_sub)$assigned_state_label, "-", "T")
V(g_tr_sub)$time[V(g_tr_sub)$name %in% tr$tip.label] = 11.5
g_sub1 = ggraph(g_tr_sub, layout = "dendrogram", height = time) +
        geom_edge_diagonal(width = 0.2, color = "black") +
        geom_node_point(aes(color = assigned_state), size = 2) +
        scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        geom_node_text(aes(label = assigned_state_label), nudge_y = 0.5, size = 4) +
        scale_color_manual(values = gr3_col) +
        ylim(c(type_graph$target_time, 0)) + ylab("")
push_pdf(g_sub1 + g_theme,
         "tr_state_subtree_v1",
         ps = 12,
         width = 4., h = 3, dir = plot_dir)

# Figure 8 plots, MRCA, correlation
tr_out = phylotime(exp_params$data[[j]]$sc,
                   mut_p = exp_params$mut_p[[j]],
                   t_total = 11.5 - 0.6,
                   return_dist = T)
dmat = dist_df2mat(tr_out)
tr3 = exp_params$tr3[[j]]
dmat_col = make_heatmap_col(max_val = 22., mid_val = 18.5, min_val = 15., col = rev(RColorBrewer::brewer.pal(3, "RdPu")))
h = Heatmap(dmat[tr3$tip.label, tr3$tip.label],
            col = dmat_col,
            show_heatmap_legend = F,
            cluster_rows = as.dendrogram(tr3),
            cluster_columns = as.dendrogram(tr3),
            show_row_names = F,
            show_row_dend = F,
            show_column_names = F,
            show_column_dend = F)

h_legend = recordPlot({
        plot.new()
        draw(Legend(col_fun = dmat_col, title = "Time (days)"))
})
push_pdf(h_legend, "phylo_dist_mat_legend", dir = plot_dir)
push_png(h, "phylo_dist_mat", dir = plot_dir, w = 3, h = 3, ps = 10, res = 1200)

plot_barcodes(exp_params$sc[[j]],
              tr = tr3,
              tip_celltype = get_type_from_id(tr3$tip.label),
              celltype_col = gr_col,
              show_row_dend = F) %>%
        push_pdf("barcodes", w = 3.2, h = 3, dir = plot_dir)

# Figure 8D
# define dmat and dmat3, then
dmat = cophenetic(exp_params$tr[[j]])
dmat3 = cophenetic(exp_params$tr3[[j]])
dmat3 = dmat3[rownames(dmat), colnames(dmat)]
dmat_tb = tibble(dist3 = c(dmat3),
                 dist = c(dmat))
dmat_plot_tb = filter(dmat_tb, dist > 0)
dmat_plot_tb_p1 = dmat_plot_tb[dmat_plot_tb$dist <= 15, ]
dmat_plot_tb_p2 = dmat_plot_tb[dmat_plot_tb$dist > 15, ]
dmat_plot_tb_sel = bind_rows(sample_n(dmat_plot_tb_p1, 1e4),
                             sample_n(dmat_plot_tb_p2, 2e4))

(dmat_plot_tb %>%
                ggscatter(x = "dist",
                          y = "dist3",
                          size = 0.25) + # ,
                #xlim = c(15, 30),
                #ylim = c(15, 30)) +
                xlab("True Cophenetic Dist") +
                ylab("Reconstructed Cophenetic Dist") +
                stat_cor() +
                geom_abline(size = 0.5) +
                theme(text = element_text(size = 10))) %>%
        push_pdf("tr_cophenetic_scatter", w = 3, h = 2.5, ps = 10, dir = plot_dir)

# plotting tr3
tr = exp_params$tr3[[j]]
g_tr = as.igraph(tr)
tr_node_time = node.depth.edgelength(tr) + 0.6
names(tr_node_time) = c(tr$tip.label, tr$node.label)
V(g_tr)$time = tr_node_time[V(g_tr)$name]
V(g_tr)$type = get_type_from_id(V(g_tr)$name)
# V(g_tr)$time[V(g_tr)$name %in% tr$tip.label] = type_graph$target_time
V(g_tr)$time[V(g_tr)$name %in% tr$tip.label] = jitter(rep(type_graph$target_time, length(tr$tip.label)), amount = 0.8)

g_tr = set.edge.attribute(g_tr,
                          name = "ending",
                          index=  unlist(incident_edges(g_tr, tr$tip.label, mode = "in")),
                          value = T)
E(g_tr)$ending[is.na(E(g_tr)$ending)] = F

g1 = ggraph(g_tr, layout = "dendrogram", height = time) +
        geom_edge_diagonal(aes(color = ending, width = ending)) +
        geom_node_point(aes(color = type), size = 2.) +
        # scale_size_manual(values = c(2., 0.05)) +
        # scale_edge_alpha_manual(values = c(0.5, 0.25)) +
        scale_edge_color_manual(values = c("#000000", "8e8e8e")) +
        scale_edge_width_manual(values = c(0.2, 0.05)) +
        scale_color_manual(values = c(gr_col, "Unknown" = NA)) +
        ylim(c(type_graph$target_time, 0)) + ylab("")
push_pdf(g1 + g_theme, "tr3_jitter", width = 14, h = 4.4, dir = plot_dir)

# plotting dropout
drop_alleles <- function(x, frac = 0.05) {
        if (frac == 0.) {
                return(x)
        }
        x[sample(length(x), size = floor(length(x) * frac), replace = F)] = NA
        x
}
sc_err = drop_alleles(
        make_existing_error_alleles(exp_params$sc[[j]],
                                    frac = 0.05),
        frac = 0.2)
sc_err_im = im_chr(sc_err)
tr = phylotime(mat_im, t_total = 11.5 - 0.6, mut_p = exp_params$mut_p[[j]], parallel = T)
plot_barcodes(sc_err, tr = tr, show_row_dend = F)
plot_barcodes(exp_params$sc[[j]], tr = exp_params$tr[[j]], show_row_dend = F)





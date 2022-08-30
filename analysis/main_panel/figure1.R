plot_dir = "./plots/fate_map_figure/"

# rnd_indices = which(tree_panel$category == "random" & tree_panel$num_tip == 16)
# all_g_temp = map(rnd_indices, function(i) {
#         set.seed(77)
#         type_graph = tree_panel$type_graph[[i]]
#         gr_col = gr_color_v1(tree_panel$type_graph[[i]], jitter_amount = 3.0, neg = T)
#         g0 = plot_type_graph_clean_mod2(type_graph,
#                                         node_label_mapper = node_label_mapper,
#                                         node_col_mapper = function(x) gr_col[x],
#                                         show_node_text = F,
#                                         add_founder = T)
#         g0
# })
# purrr::reduce(all_g_temp, `+`)

node_label_mapper <- function(node_names) {
        out = node_names
        node_ind = !grepl("-", node_names)
        node_names[node_ind] = sum(node_ind) + 1 - as.numeric(node_names[node_ind])
        for(i in 1:length(out)) {
                if (grepl("-", node_names[i])) {
                        out[i] = str_replace(node_names[i], "-", "T")
                } else {
                        out[i] = paste0("P", node_names[i])
                }
        }
        out
}

set.seed(77)
i = 23
type_graph = tree_panel$type_graph[[i]]
gr_col = gr_color_v1(tree_panel$type_graph[[i]], jitter_amount = 3.0, neg = T)
man_col_map = c("#c7ddad" = "#9ec162",
                "#d7e2d5" = "#19e2f7",
                "#eedeac" = "#f7c852",
                "#f3dbc0" = "#ce5fac",
                "#ceb3b5" = "#84ce89")
gr_col[match(toupper(names(man_col_map)), gr_col)] = man_col_map
g0 = plot_type_graph_clean_mod2(type_graph,
                                node_label_mapper = node_label_mapper,
                                node_col_mapper = function(x) gr_col[x],
                                show_node_text = T,
                                add_founder = T) + theme(legend.position = "none")
g0
theme0 = theme(
        panel.background = element_rect(fill = NA, color = NA, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=9),
        text = element_text(size = 12)
)
theme1 = theme(
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = NA, size = 10),
        legend.text=element_text(size=9),
        text = element_text(size = 12),
        axis.title.x = element_blank()
)

gens0 = make_gens_mod2(type_graph)
gens0 = sample_size_gens_mod2(gens0, type_graph, 100)
mut_p = list(mut_rate = 1.,
             recur_prob = 1.,
             recur_vec = c(0.5, 0.2, 0.1, 0.1, 0.05, 0.05),
             active_time = list(c(0, 3.)))
mut_param = do.call(make_mut_param_by_rate, mut_p)
gens0 = simulate_muts(gens0, type_graph, mut_param = mut_param)
# mut_counts_type = lapply(gens0, extract_mut_counts_history, type_graph = type_graph)
# mut_counts_history = merge_mut_counts(mut_counts_type)

# plot_cell_count_graph <- function(type_graph, gr_col) {
#         gens0 = make_gens_mod2(type_graph)
#         gens0 = sample_size_gens_mod2(gens0, type_graph, 100)
#         g = make_empty_graph()
#         for(type_gens in gens0[c(type_graph$tip_id, backward_merge_sequence_mod2(type_graph))]) {
#                 message(type_gens$cell_type)
#                 # ignore cases of unsampled branches for now
#                 # when there's no unsampled branching
#                 # vertices
#                 if (type_gens$num_gen <= 1) {
#                         i = 0
#                         g = g %>% add_vertices(1,
#                                                name = paste0("type_", type_gens$cell_type, "_gen_", i),
#                                                count = sum(type_gens$cell_mode_counts[i+1, c(1, 3)]),
#                                                time = type_gens$start_time + type_gens$double_time * i,
#                                                cell_type =type_gens$cell_type,
#                                                vertex_type = "cells")
#                         if (!is.null(type_gens$sample_size)) {
#                                 V(g)[paste0("type_", type_gens$cell_type, "_gen_", i)]$sample_size = type_gens$sample_size[i+1]
#                         }
#                 } else {
#                         for (i in c(0, type_gens$num_gen)) {
#                                 g = g %>% add_vertices(1,
#                                                        name = paste0("type_", type_gens$cell_type, "_gen_", i),
#                                                        count = sum(type_gens$cell_mode_counts[i+1, c(1, 3)]),
#                                                        time = type_gens$start_time + type_gens$double_time * i,
#                                                        cell_type =type_gens$cell_type,
#                                                        vertex_type = "cells")
#                                 if (!is.null(type_gens$sample_size)) {
#                                         V(g)[paste0("type_", type_gens$cell_type, "_gen_", i)]$sample_size = type_gens$sample_size[type_gens$num_gen - i + 1]
#                                 }
#                         }
#                 }
#                 if (type_gens$num_gen > 0) {
#                         # one edge connecting first and last gen
#                         g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", 0),
#                                               paste0("type_", type_gens$cell_type, "_gen_", c(type_gens$num_gen - 1))),
#                                             edge_type = "cell division")
#                 }
#                 # cell type transition vertives and edges
#                 if (type_gens$active) {
#                         biases = type_gens$end_mode_counts[1:2] / sum(type_gens$end_mode_counts[1:2])
#                         out_types = type_graph$merge[[type_gens$cell_type]]
#                         out_gen_names = paste0("type_", out_types, "_gen_0")
#                         g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", pmax(type_gens$num_gen -1, 0)),
#                                               out_gen_names[1]),
#                                             edge_type = "cell fate commit",
#                                             bias = biases[1]) %>%
#                                 add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", pmax(type_gens$num_gen -1, 0)),
#                                             out_gen_names[2]),
#                                           edge_type = "cell fate commit",
#                                           bias = biases[2])
#                 }
#         }
#         V(g)$count_chr = map_chr(V(g)$count, function(x) {
#                 if (x <= 9999) {
#                         return(as.character(x))
#                 } else {
#                         return(format(x, digits = 2, justify = "right", scientific = T))
#                 }
#         })
#         g1 = ggraph(g, layout = "dendrogram", height = time) +
#                 geom_edge_bend(aes(linetype = edge_type, edge_width = bias)) +
#                 geom_node_point(aes(shape = vertex_type, color = cell_type, size = log2(count+1))) +
#                 geom_node_text(aes(label = count_chr), nudge_y = -0.15, size = 2.) +
#                 ylim(c(type_graph$target_time, 0)) +
#                 scale_size_continuous(range = c(4, 9), guide = F) +
#                 scale_color_manual(values = gr_col, na.value = "grey") +
#                 scale_edge_width(range = c(0., 1.), na.value = 0.5) +
#                 theme0
#         g1
# }

make_cell_count_graph_v1 <- function(type_graph, gr_col) {
        gens0 = make_gens_mod2(type_graph)
        gens0 = sample_size_gens_mod2(gens0, type_graph, 100)
        g = make_empty_graph()
        for(type_gens in gens0[c(type_graph$tip_id, backward_merge_sequence_mod2(type_graph))]) {
                message(type_gens$cell_type)
                # ignore cases of unsampled branches for now
                # when there's no unsampled branching
                # vertices
                # if (type_gens$num_gen <= 1) {
                #         i = 0
                #         g = g %>% add_vertices(1,
                #                                name = paste0("type_", type_gens$cell_type, "_gen_", i),
                #                                count = sum(type_gens$cell_mode_counts[i+1, c(1, 3)]),
                #                                time = type_gens$start_time + type_gens$double_time * i,
                #                                cell_type =type_gens$cell_type,
                #                                vertex_type = "cells")
                #         if (!is.null(type_gens$sample_size)) {
                #                 V(g)[paste0("type_", type_gens$cell_type, "_gen_", i)]$sample_size = type_gens$sample_size[i+1]
                #         }
                # } else {
                        for (i in unique(c(0, type_gens$num_gen))) {
                                if (i == type_gens$num_gen & !type_gens$active) {
                                        # only time is different
                                        g = g %>% add_vertices(1,
                                                               name = paste0("type_", type_gens$cell_type, "_gen_", i),
                                                               count = sum(type_gens$cell_mode_counts[i+1, c(1, 3)]),
                                                               time = type_graph$target_time, # using target time here
                                                               cell_type =type_gens$cell_type,
                                                               vertex_type = "cells",
                                                               class = "terminal")
                                } else {
                                        if (type_gens$cell_type == type_graph$root_id & i == 0) {
                                                class_chr = "root"
                                        } else {
                                                class_chr = ifelse(type_gens$num_gen == i, "pre", "post")
                                        }
                                        g = g %>% add_vertices(1,
                                                               name = paste0("type_", type_gens$cell_type, "_gen_", i),
                                                               count = sum(type_gens$cell_mode_counts[i+1, c(1, 3)]),
                                                               time = type_gens$start_time + type_gens$double_time * i,
                                                               cell_type =type_gens$cell_type,
                                                               vertex_type = "cells",
                                                               class = class_chr)
                                }
                                if (!is.null(type_gens$sample_size)) {
                                        V(g)[paste0("type_", type_gens$cell_type, "_gen_", i)]$sample_size = type_gens$sample_size[type_gens$num_gen - i + 1]
                                }
                        }
                # }
                if (type_gens$num_gen > 0) {
                        # one edge connecting first and last gen
                        g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", 0),
                                              paste0("type_", type_gens$cell_type, "_gen_", c(type_gens$num_gen))),
                                            edge_type = "cell division")
                }
                # cell type transition vertives and edges
                if (type_gens$active) {
                        biases = type_gens$end_mode_counts[1:2] / sum(type_gens$end_mode_counts[1:2])
                        out_types = type_graph$merge[[type_gens$cell_type]]
                        out_gen_names = paste0("type_", out_types, "_gen_0")
                        g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", pmax(type_gens$num_gen, 0)),
                                              out_gen_names[1]),
                                            edge_type = "cell fate commit",
                                            bias = biases[1]) %>%
                                add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", pmax(type_gens$num_gen, 0)),
                                            out_gen_names[2]),
                                          edge_type = "cell fate commit",
                                          bias = biases[2])
                }
        }
        V(g)$count_chr = map_chr(V(g)$count, function(x) {
                if (x <= 9999) {
                        return(as.character(x))
                } else {
                        return(format(x, digits = 2, justify = "right", scientific = T))
                }
        })
        g
}

g = make_cell_count_graph_v1(type_graph, gr_col)
g1 = ggraph(g, layout = "dendrogram", height = time) +
        geom_node_point(aes(shape = vertex_type,
                            color = cell_type,
                            size = log2(count+1))) +
        geom_edge_bend(aes(color = edge_type, edge_width = bias)) +
        geom_node_text(aes(label = count_chr), nudge_y = +0.15, size = 2.) +
        ylim(c(type_graph$target_time, 0)) +
        scale_size_continuous(range = c(2, 5), guide = F) +
        scale_color_manual(values = gr_col) +
        scale_edge_width(range = c(0.5, 2.), na.value = 0.5) +
        scale_edge_color_manual(values = c("#fb9a99", "#8dd3c7")) +
        theme0

# g1
# library(ggraph)
# ggraph(g, layout = "dendrogram", height = time) + geom_node_label(aes(label = name, color = cell_type)) + geom_edge_bend() + ylim(c(type_graph$target_time, 0))
# g1 = ggraph(g, layout = "dendrogram", height = time) +
#         geom_edge_bend(aes(linetype = edge_type, edge_width = bias)) +
#         geom_node_point(aes(shape = vertex_type, color = cell_type, size = log2(count+1))) +
#         geom_node_text(aes(label = count_chr), nudge_y = -0.15, size = 2.) +
#         ylim(c(type_graph$target_time, 0)) +
#         scale_size_continuous(range = c(6, 9), guide = F) +
#         scale_color_manual(values = gr_color(type_graph), na.value = "grey") +
#         scale_edge_width(range = c(0., 1.), na.value = 0.5) +
#         theme0
# push_pdf(g1, file_name = "cell_count", width = 8.5, h = 3.5, ps = 12, dir = plot_dir)

V(g)$unsampled_count = V(g)$count - V(g)$sample_size
xy <- ggraph::create_layout(g, "dendrogram", height = time)
V(g)$x <- xy[match(V(g)$name, xy$name), 1]
V(g)$y <- xy[match(V(g)$name, xy$name), 2]

frac_vec = V(g)$sample_size / V(g)$count
frac_chr = map_chr(frac_vec, function(x) {
        if (x < 0.01) {
                return("<0.01")
        } else {
                if (x == 1) {
                        return("1.00")
                }
                return(format(x, digits = 2, justify = "right", scientific = F))
        }
})

V(g)$sampled_fraction = frac_chr

# sample size
# library(scatterpie)
# dd = igraph::as_data_frame(g, "vertices")
# dd$sampled_color = toy_col_mapper(dd$cell_type)
# dd$sampled_color = toy_col_mapper(dd$cell_type)

aspect_ratio = 0.5
g2 = ggraph(g, layout = "manual", x = V(g)$x / aspect_ratio, y = type_graph$target_time - V(g)$y) +
        geom_edge_bend(aes(color = edge_type, edge_width = bias)) +
        geom_scatterpie(aes(x = x / aspect_ratio,
                            y =  type_graph$target_time - y,
                            r = log2(count+1)/50+0.2),
                        cols = c("unsampled_count", "sample_size"),
                        color = NA,
                        data = filter(igraph::as_data_frame(g, "vertices"), class != "post")) +
        scale_fill_manual(values = c("#d1d1d1", "#3182bd")) +
        scale_edge_width(range = c(0.5, 2.), na.value = 0.5) +
        scale_edge_color_manual(values = c("#fb9a99", "#8dd3c7")) +
        # scale_edge_linetype(guide=F) +
        geom_node_text(aes(label = count_chr),
                       nudge_y = 0.4,
                       size = 2.) +
        ylim(c(-1.0, type_graph$target_time + 0.5)) +
        ylab("Time") +
        coord_fixed() + theme0 + theme(legend.position = "right")
# push_pdf(g2, file_name = "sampled_fraction", width = 8.5, h = 3.5, ps = 12, dir = plot_dir)
((g0 + g2) + plot_layout(guide = "collect", width = c(1, 2.5))) %>%
        push_pdf(file_name = "figure1_example_v1", dir = plot_dir,
                 width = 11.5, height = 3.0,
                 ps = 10)

g_lout = create_layout(g, layout = "manual", x = V(g)$x, y = 11.5 - V(g)$y)
g_lout = g_lout[g_lout$time > 11., ]
g_lout$name[order(g_lout$x)]

# proportional sampling
tip_size = map_dbl(gens0[type_graph$tip_id], "end_count")
tip_sample = rmvhyper(nn = 1, k = length(type_graph$tip_id) * 100, n = tip_size)[1, ]
tip_sample = pmax(tip_sample, 5)
names(tip_sample) = type_graph$tip_id
tip_sample[map_chr(strsplit(g_lout$name[order(g_lout$x)], "_"), 2)]










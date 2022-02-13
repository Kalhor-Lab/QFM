# toy example in figure 2
plot_dir = "./plots/toy/"
g_theme <-                 theme(
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = NA, size = 10),
        legend.position = 'none',
        axis.title.x = element_blank(),
        text = element_text(size = 12)
)

output_dir = "./intermediate_data/panel/"
exp_params = readRDS(paste0(output_dir, "exp_data_rep1_2_proc2.rds"))
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))
# experiment 9 used as an example
j = 9
type_graph = all_graphs[[exp_params$big_graph_id[j]]]

set.seed(36)
gens0 = make_gens(type_graph)
gens0 = sample_size_gens(gens0, type_graph, 100)
mut_p = list(mut_rate = 1.,
             recur_prob = 1.,
             recur_vec = c(0.5, 0.2, 0.1, 0.1, 0.05, 0.05),
             active_time = list(c(0, 3.)))
mut_param = do.call(make_mut_param_by_rate, mut_p)
gens0 = simulate_muts(gens0, type_graph, mut_param = mut_param)
# mut_counts_type = lapply(gens0, extract_mut_counts_history, type_graph = type_graph)
# mut_counts_history = merge_mut_counts(mut_counts_type)

g = make_empty_graph()
for(type_gens in gens0[c(tip_id, backward_merge_sequence(type_graph))]) {
        message(type_gens$cell_type)
        # ignore cases of unsampled branches for now
        # when there's no unsampled branching
        # vertices
        if (type_gens$num_gen <= 1) {
                i = 0
                g = g %>% add_vertices(1,
                                       name = paste0("type_", type_gens$cell_type, "_gen_", i),
                                       count = type_gens$cell_counts[i+1],
                                       time = type_gens$start_time + type_gens$double_time * i,
                                       cell_type =type_gens$cell_type,
                                       vertex_type = "cells")
                if (!is.null(type_gens$sample_size)) {
                        V(g)[paste0("type_", type_gens$cell_type, "_gen_", i)]$sample_size = type_gens$sample_size[i+1]
                }
        } else {
                for (i in c(0, type_gens$num_gen - 1)) {
                        g = g %>% add_vertices(1,
                                               name = paste0("type_", type_gens$cell_type, "_gen_", i),
                                               count = type_gens$cell_counts[i+1],
                                               time = type_gens$start_time + type_gens$double_time * i,
                                               cell_type =type_gens$cell_type,
                                               vertex_type = "cells")
                        if (!is.null(type_gens$sample_size)) {
                                V(g)[paste0("type_", type_gens$cell_type, "_gen_", i)]$sample_size = type_gens$sample_size[i+1]
                        }
                }
        }
        if (type_gens$num_gen > 0) {
                # one edge connecting first and last gen
                g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", 0),
                                      paste0("type_", type_gens$cell_type, "_gen_", c(type_gens$num_gen - 1))),
                                    edge_type = "cell division")
        }
        # cell type transition vertives and edges
        if (type_gens$active) {
                biases = type_gens$end_mode_counts[1:2] / sum(type_gens$end_mode_counts[1:2])
                out_types = type_graph$merge[type_gens$cell_type, ]
                out_gen_names = paste0("type_", out_types, "_gen_0")
                g = g %>% add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", pmax(type_gens$num_gen -1, 0)),
                                      out_gen_names[1]),
                                    edge_type = "cell fate commit",
                                    bias = biases[1]) %>%
                        add_edges(c(paste0("type_", type_gens$cell_type, "_gen_", pmax(type_gens$num_gen -1, 0)),
                                    out_gen_names[2]),
                                  edge_type = "cell fate commit",
                                  bias = biases[2])
        }
}

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

V(g)$count_chr = map_chr(V(g)$count, function(x) {
        if (x <= 9999) {
                return(as.character(x))
        } else {
                return(format(x, digits = 2, justify = "right", scientific = T))
        }
})

library(ggraph)
# ggraph(g, layout = "dendrogram", height = time) + geom_node_label(aes(label = name, color = cell_type)) + geom_edge_bend() + ylim(c(type_graph$target_time, 0))
g1 = ggraph(g, layout = "dendrogram", height = time) +
        geom_edge_bend(aes(linetype = edge_type, edge_width = bias)) +
        geom_node_point(aes(shape = vertex_type, color = cell_type, size = log2(count+1))) +
        geom_node_text(aes(label = count_chr), nudge_y = -0.15, size = 2.) +
        ylim(c(type_graph$target_time, 0)) +
        scale_size_continuous(range = c(6, 9), guide = F) +
        scale_color_manual(values = gr_color(type_graph), na.value = "grey") +
        scale_edge_width(range = c(0., 1.), na.value = 0.5) +
        theme0
push_pdf(g1, file_name = "cell_count", width = 8.5, h = 3.5, ps = 12, dir = plot_dir)

V(g)$unsampled_count = V(g)$count - V(g)$sample_size
xy <- ggraph::create_layout(g, "dendrogram", height = time)
V(g)$x <- xy[match(V(g)$name, xy$name), 1]
V(g)$y <- xy[match(V(g)$name, xy$name), 2]

frac_vec = V(g)$sample_size / V(g)$count
frac_chr = map_chr(frac_vec, function(x) {
        if (x < 0.01) {
                return(format(x, digits = 2, scientific = T))
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

create_layout(g, layout = "manual", x = V(g)$x, y = 15 - V(g)$y)
igraph::as_data_frame(g, "vertices")

g2 = ggraph(g, layout = "manual", x = V(g)$x, y = 15. - V(g)$y) +
        geom_edge_bend(aes(linetype = edge_type)) +
        geom_scatterpie(aes(x = x,
                            y =  15. - y,
                            r = log2(count+1)/160+0.2),
                        cols = c("unsampled_count", "sample_size"),
                        color = NA,
                        data = igraph::as_data_frame(g, "vertices")) +
        scale_fill_manual(values = c("#d1d1d1", "black")) +
        geom_node_text(aes(label = sampled_fraction),
                       nudge_y = -0.4,
                       size = 2.) +
        ylab("Time") +
        scale_edge_linetype(guide=F) + theme1
push_pdf(g2, file_name = "sampled_fraction", width = 8.5, h = 3.5, ps = 12, dir = plot_dir)

# proportional sampling
tip_size = map_dbl(gens0[type_graph$tip_id], "end_count")
tip_sample = rmvhyper(nn = 1, k = length(type_graph$tip_id) * 100, n = tip_size)[1, ]
tip_sample = pmax(tip_sample, 5)
names(tip_sample) = fm$tip_id
tip_sample

# push_pdf(g_all_simple[[9]],
#          file_name = "gr_simple",
#          width = 6.5, h = 3.5, ps = 12, dir = plot_dir)





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
        names(out) = node_names
        out
}
map_gr_label <- function(lab) {
        out = str_replace(lab, "Node-", "iP")
        out = str_replace(out, "-", "T")
        names(out) = lab
        out
}
g_theme = theme(
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        text = element_text(size = 12)
)
generate_mut_color <- function(mat) {
        unique_val = sort(unique(c(mat)))
        unique_val = unique_val[unique_val!="0"]
        set.seed(76)
        mut_colors = ComplexHeatmap:::default_col(factor(1:(length(unique_val))))
        stopifnot(mut_colors != "#000000")
        colors = structure(c("#000000", mut_colors),
                           names = c("0", unique_val))
        colors
}
plot_barcodes <- function(sc_mat,
                          tr,
                          tip_celltype = NULL,
                          celltype_col = NULL,
                          tip_batch = NULL,
                          mut_col = NULL,
                          id_ord = NULL,
                          ...) {
        # make sure each element has unique allele names
        if (is.null(mut_col)) {
                sc_mat_mod = rlang::duplicate(sc_mat)
                for (j in 1:ncol(sc_mat_mod)) {
                        sc_mat_mod[, j] = paste0("el_", j, "_", sc_mat_mod[, j])
                        sc_mat_mod[is.na(sc_mat[, j]), j] = NA
                        sc_mat_mod[sc_mat[, j] == "0", j] = "0"
                }
                colors = generate_mut_color(sc_mat_mod)
        }
        else {
                sc_mat_mod = sc_mat
                colors = mut_col
        }
        if (is.null(id_ord)) {
                id_ord = order(colMeans(sc_mat_mod == "0", na.rm = T))
        }
        heat = Heatmap(sc_mat_mod[tr$tip.label, id_ord],
                       col = colors,
                       cluster_rows = phylogram::as.dendrogram(tr),
                       cluster_columns = F,
                       show_row_names = F,
                       show_heatmap_legend = F,
                       ...)
        if (!is.null(tip_celltype)) {
                if (is.null(names(tip_celltype))) {
                        names(tip_celltype) = rownames(sc_mat)
                }
                if (!is.null(celltype_col)) {
                        if (!is.null(tip_batch)) {
                                ra = rowAnnotation(celltype = tip_celltype[tr$tip.label],
                                                   batch = tip_batch[tr$tip.label],
                                                   col = list(celltype = celltype_col))
                        } else {
                                ra = rowAnnotation(celltype = tip_celltype[tr$tip.label],
                                                   col = list(celltype = celltype_col))
                        }
                } else {
                        ra = rowAnnotation(celltype = tip_celltype[tr$tip.label])
                }
                return(ra + heat)
        } else {
                return(heat)
        }
}
make_ig_from_tr <- function(tr, add_root = F) {
        g_tr = as.igraph(tr)
        V(g_tr)$type = get_type_from_id(V(g_tr)$name)
        V(g_tr)$time = get_node_time(tr, V(g_tr)$name)
        if (add_root) {
                g_tr = g_tr %>% add_vertices(1, name = "root", type = "root", time = 0) %>% add_edges(c("root",
                                                                                                        "Node-1"))
        }
        g_tr
}
scale_min_max <- function(x, min, max) {
        x_norm  = (x- min(x)) /(max(x)-min(x))
        x_norm * (max - min) + min
}
gr_color <- function(gr, jitter_amount = 2., neg = F) {
        if (class(gr) == "type_graph") {
                edges_tb = gr$edges
        }
        if (class(gr) == "phylo") {
                edges_tb = make_edge_tb(gr)[1:3]
        }
        # make an embeding  and color scheme
        ig_gr = graph_from_data_frame(edges_tb, directed = F)
        E(ig_gr)$weight = E(ig_gr)$length

        # scale coordinates
        gr_dist = distances(ig_gr)
        # gr_embed = uwot::umap(gr_dist, n_components = 3, repulsion_strength = 10.)
        # rownames(gr_embed) = rownames(gr_dist)
        gr_embed = cmdscale(distances(ig_gr), k = 3)
        if (neg) {
                gr_embed[, 1] = -gr_embed[, 2]
                gr_embed[, 2] = -gr_embed[, 3]
                gr_embed[, 3] = -gr_embed[, 1]
        }
        gr_embed[, 1] = scale_min_max(jitter(gr_embed[, 1], amount = jitter_amount), 0.1, 0.9)
        gr_embed[, 2] = scale_min_max(jitter(gr_embed[, 2], amount = jitter_amount), 0.1, 0.9)
        gr_embed[, 3] = scale_min_max(jitter(gr_embed[, 3], amount = jitter_amount), 0.1, 0.9)
        col = colorspace::hex(colorspace::RGB(gr_embed))
        # qplot(gr_embed[, 1], gr_embed[, 2], col = names(col)) + scale_color_manual(values = col)
        col
}
gr_color_v1 <- function (gr, jitter_amount = 2, neg = F)
{
        if (class(gr) == "type_graph") {
                edges_tb = gr$edges
        }
        if (class(gr) == "phylo") {
                edges_tb = make_edge_tb(gr)[1:3]
        }
        edges_tb$length = 1.
        ig_gr = graph_from_data_frame(edges_tb, directed = F)

        n_all = length(V(ig_gr))
        if (n_all <= 4) {
                col = RColorBrewer::brewer.pal(n_all, name = "Set2")
                names(col) = c(gr$node.label, gr$tip.label)
                return(col)
        }

        E(ig_gr)$weight = E(ig_gr)$length
        gr_dist = distances(ig_gr)
        gr_embed = cmdscale(distances(ig_gr), k = 3)
        if (neg) {
                gr_embed[, 1] = -gr_embed[, 2]
                gr_embed[, 2] = -gr_embed[, 3]
                gr_embed[, 3] = -gr_embed[, 1]
        }
        gr_embed[, 1] = scale_min_max(jitter(gr_embed[, 1], amount = jitter_amount),
                                      0.1, 0.9)
        gr_embed[, 2] = scale_min_max(jitter(gr_embed[, 2], amount = jitter_amount),
                                      0.1, 0.9)
        gr_embed[, 3] = scale_min_max(jitter(gr_embed[, 3], amount = jitter_amount),
                                      0.1, 0.9)
        col = colorspace::hex(colorspace::RGB(gr_embed))
        col
}
plot_topology <- function(gr, gr_node_time = NULL, total_time = NULL,
                          type_col = NULL, node_label = NULL,
                          edge_col = NULL,
                          node_size = 3, edge_width = 0.5,
                          gr_node_cat = NULL,
                          show_node_label = F) {
        ig_gr = as.igraph(gr)
        if (is.null(gr_node_time)) {
                gr_node_time = node.depth.edgelength(gr)
                names(gr_node_time) = c(gr$tip.label, gr$node.label)
                V(ig_gr)$Time = gr_node_time[V(ig_gr)$name]
        }
        V(ig_gr)$Time = gr_node_time[V(ig_gr)$name]
        V(ig_gr)$type = V(ig_gr)$name

        if (is.null(total_time)) {
                total_time = max(gr_node_time)
        }
        # V(ig_gr)$Type = V(ig_gr)$name
        if (!is.null(node_label)) {
                V(ig_gr)$node_label = node_label[V(ig_gr)$name]
        } else {
                V(ig_gr)$node_label = V(ig_gr)$name
        }
        if (!is.null(gr_node_cat)) {
                V(ig_gr)$Cat = gr_node_cat[V(ig_gr)$name]
        } else {
                V(ig_gr)$Cat = "Tip"
                V(ig_gr)$Cat[V(ig_gr)$name %in% gr$node.label] = "Node"
        }
        if (!is.null(edge_col)) {
                E(ig_gr)$type = igraph::as_data_frame(ig_gr)$to
        }
        # V(ig_gr)$node_label = rep("", length(V(ig_gr)$name))
        if (is.null(type_col)) {
                type_col = gr_color(gr)
        }
        g_out = ggraph(ig_gr, layout = 'dendrogram', height = Time)
        if (!is.null(edge_col)) {
                g_out = g_out +
                        geom_edge_bend(
                                aes(color = type,
                                    width = factor(1),
                                    fontface = 'plain'),
                                arrow = arrow(
                                        type = "closed",
                                        length = unit(2, "pt"),
                                        angle = 45
                                ),
                                start_cap = circle(5, 'pt'),
                                end_cap = circle(5, 'pt')
                        )
                g_out = g_out + scale_edge_color_manual(values = edge_col)
        } else {
                g_out = g_out +
                        geom_edge_bend(
                                aes(width = factor(1),
                                    fontface = 'plain'),
                                arrow = arrow(
                                        type = "closed",
                                        length = unit(2, "pt"),
                                        angle = 45
                                ),
                                start_cap = circle(5, 'pt'),
                                end_cap = circle(5, 'pt')
                        )
        }
        g_out = g_out +
                geom_node_point(aes(color = type, shape = Cat),
                                size = node_size)
        if (show_node_label) {
                g_out = g_out + geom_node_text(aes(label = node_label), nudge_x = - 0.2, nudge_y = 0.3)
        }
        if (!is.null(gr_node_cat)) {
                assertthat::assert_that(length(unique(gr_node_cat)) == 2)
                g_out  = g_out + scale_shape_manual(values = c(18, 17))
        } else {
                g_out  = g_out + scale_shape_manual(values = c(18, 17))
        }
        g_out = g_out +
                # scale_fill_manual(values = type_col) +
                scale_color_manual(values = type_col) +
                ylim(c(total_time, -0.5)) + ylab('Time') +
                scale_edge_width_manual(values = edge_width, guide = "none") +
                scale_size_manual(values = 4, guide = "none") +
                theme(
                        axis.line.y = element_line(),
                        axis.ticks.y = element_line(),
                        axis.text.y = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        panel.background = element_rect(fill = NA, color = NA, size = 10),
                        legend.position = 'none',
                        axis.title.x = element_blank(),
                        text = element_text(size = 12)
                )
        g_out
}
# alias
plot_gr_clean <- plot_topology
plot_gr_dendro <- function(gr, gr_node_time, type_col, total_time, node_size = 3, gr_node_cat = NULL, plot_node_point = T) {
        ig_gr = as.igraph(gr)
        # gr_node_time = node.depth.edgelength(gr)
        # names(gr_node_time) = c(gr$tip.label, gr$node.label)
        # V(ig_gr)$Time = gr_node_time[V(ig_gr)$name] + gr_root_time
        V(ig_gr)$Time = gr_node_time[V(ig_gr)$name]
        V(ig_gr)$type = V(ig_gr)$name
        V(ig_gr)$node_label = V(ig_gr)$name
        if (!is.null(gr_node_cat)) {
                V(ig_gr)$Cat = gr_node_cat[V(ig_gr)$name]
        } else {
                V(ig_gr)$Cat = "None"
        }
        # V(ig_gr)$node_label = rep("", length(V(ig_gr)$name))

        g_out = ggraph(ig_gr, layout = 'dendrogram', height = Time) +
                geom_edge_elbow()
        if (plot_node_point) {
                g_out = g_out +
                        geom_node_point(aes(color = type, shape = Cat),
                                        size = node_size)
        }
        # geom_node_label(aes(
        #         label = node_label,
        #         fill = name,
        #         size = factor(1)
        # ),
        # label.padding = unit(6, 'pt'))\
        if (!is.null(gr_node_cat)) {
                assertthat::assert_that(length(unique(gr_node_cat)) == 2)
                g_out  = g_out + scale_shape_manual(values = c(18, 15))
        } else {
                g_out  = g_out + scale_shape_manual(values = c(15))
        }
        g_out = g_out +
                # scale_fill_manual(values = type_col) +
                scale_color_manual(values = type_col) +
                ylim(c(total_time, 0)) + ylab('Time') +
                # scale_edge_width_manual(values = 0.5, guide = "none") +
                scale_size_manual(values = 4, guide = "none") +
                theme(
                        axis.line.y = element_line(),
                        axis.ticks.y = element_line(),
                        axis.text.y = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        panel.background = element_rect(fill = NA, color = NA, size = 10),
                        legend.position = 'none',
                        axis.title.x = element_blank(),
                        text = element_text(size = 12)
                )
        g_out
}
plot_tr <- function(tr,
                    node_types = NULL,
                    root_time = 0.0,
                    type_col = NULL,
                    total_time = NULL,
                    node_size = 0.75,
                    edge_alpha = 0.5,
                    end_alpha_terminal = 0.05,
                    edge_width = 0.2,
                    edge_width_terminal = 0.02,
                    ylim = c(total_time, 0),
                    jitter_tip = F
                    ) {
        g_tr = as.igraph(tr)
        tr_node_time = ape::node.depth.edgelength(tr) + root_time
        names(tr_node_time) = c(tr$tip.label, tr$node.label)
        V(g_tr)$time = tr_node_time[V(g_tr)$name]
        V(g_tr)$type = "Unknown"
        # setting node types
        if (!is.null(node_types)) {
                assertthat::assert_that(all(names(node_types) %in% V(g_tr)$name))
                V(g_tr)$type[match(names(node_types), V(g_tr)$name)] = node_types
        }
        if (!is.null(total_time)) {
                V(g_tr)$time[V(g_tr)$name %in% tr$tip.label] = total_time
        } else {
                total_time = max(V(g_tr)$time)
        }
        g_tr = set.edge.attribute(g_tr,
                                  name = "ending",
                                  index=  unlist(incident_edges(g_tr, tr$tip.label, mode = "in")),
                                  value = T)
        E(g_tr)$ending[is.na(E(g_tr)$ending)] = F
        unknown_col = "#808080"
        names(unknown_col) = "Unknown"
        if (!is.null(type_col)) {
                type_col = c(type_col, unknown_col)
        } else {
                type_col = unknown_col
        }
        if (jitter_tip) {
                V(g_tr)$time[V(g_tr)$name %in% tr$tip.label] = jitter(rep(total_time,
                                                                          length(tr$tip.label)), amount = 0.5)
        }
        g_out = ggraph(g_tr, layout = "dendrogram", height = time) +
                geom_edge_diagonal(aes(alpha = ending, width = ending), color = "black") +
                geom_node_point(aes(color = type), size = node_size) +
                scale_edge_alpha_manual(values = c(edge_alpha, end_alpha_terminal)) +
                scale_edge_width_manual(values = c(edge_width, edge_width_terminal)) +
                scale_color_manual(values = type_col) +
                ylim(ylim) + ylab("")
        g_out + g_theme
}


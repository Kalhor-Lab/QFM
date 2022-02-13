g_theme <- theme(
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
plot_gr_clean <- function(gr3, gr_node_time, node_col, target_time, node_size = 3, gr_node_cat = NULL) {
        ig_gr3 = as.igraph(gr3)
        # gr_node_time = node.depth.edgelength(gr3)
        # names(gr_node_time) = c(gr3$tip.label, gr3$node.label)
        # V(ig_gr3)$Time = gr_node_time[V(ig_gr3)$name] + gr_root_time
        V(ig_gr3)$Time = gr_node_time[V(ig_gr3)$name]
        V(ig_gr3)$type = V(ig_gr3)$name
        V(ig_gr3)$node_label = V(ig_gr3)$name
        if (!is.null(gr_node_cat)) {
                V(ig_gr3)$Cat = gr_node_cat[V(ig_gr3)$name]
        } else {
                V(ig_gr3)$Cat = "None"
        }
        # V(ig_gr3)$node_label = rep("", length(V(ig_gr3)$name))

        g_out = ggraph(ig_gr3, layout = 'dendrogram', height = Time) +
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
                ) +
                geom_node_point(aes(color = type, shape = Cat),
                                size = node_size)
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
                # scale_fill_manual(values = node_col) +
                scale_color_manual(values = node_col) +
                ylim(c(target_time, 0)) + ylab('Time') +
                scale_edge_width_manual(values = 0.5, guide = "none") +
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

plot_gr <- function(gr3, gr_node_time, node_col, target_time, node_label = NULL) {
        ig_gr3 = as.igraph(gr3)
        # gr_node_time = node.depth.edgelength(gr3)
        # names(gr_node_time) = c(gr3$tip.label, gr3$node.label)
        # V(ig_gr3)$Time = gr_node_time[V(ig_gr3)$name] + gr_root_time
        V(ig_gr3)$Time = gr_node_time[V(ig_gr3)$name]
        V(ig_gr3)$type = V(ig_gr3)$name
        if (!is.null(node_label)) {
                V(ig_gr3)$node_label = node_label[V(ig_gr3)$name]
        } else {
                V(ig_gr3)$node_label = V(ig_gr3)$name
        }

        # V(ig_gr3)$node_label = rep("", length(V(ig_gr3)$name))
        g_out = ggraph(ig_gr3, layout = 'dendrogram', height = Time) +
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
                ) +
                geom_node_label(aes(
                        label = node_label,
                        fill = name,
                        size = factor(1)
                ),
                label.padding = unit(6, 'pt'))
        g_out = g_out +
                scale_fill_manual(values = node_col) +
                ylim(c(target_time, 0)) + ylab('Time') +
                scale_edge_width_manual(values = 0.5, guide = "none") +
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

plot_gr_dendro <- function(gr3, gr_node_time, node_col, target_time, node_size = 3, gr_node_cat = NULL) {
        ig_gr3 = as.igraph(gr3)
        # gr_node_time = node.depth.edgelength(gr3)
        # names(gr_node_time) = c(gr3$tip.label, gr3$node.label)
        # V(ig_gr3)$Time = gr_node_time[V(ig_gr3)$name] + gr_root_time
        V(ig_gr3)$Time = gr_node_time[V(ig_gr3)$name]
        V(ig_gr3)$type = V(ig_gr3)$name
        V(ig_gr3)$node_label = V(ig_gr3)$name
        if (!is.null(gr_node_cat)) {
                V(ig_gr3)$Cat = gr_node_cat[V(ig_gr3)$name]
        } else {
                V(ig_gr3)$Cat = "None"
        }
        # V(ig_gr3)$node_label = rep("", length(V(ig_gr3)$name))

        g_out = ggraph(ig_gr3, layout = 'dendrogram', height = Time) +
                geom_edge_elbow() +
                geom_node_point(aes(color = type, shape = Cat),
                                size = node_size)
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
                # scale_fill_manual(values = node_col) +
                scale_color_manual(values = node_col) +
                ylim(c(target_time, 0)) + ylab('Time') +
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

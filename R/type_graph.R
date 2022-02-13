mode_offspring = do.call(rbind, list(c(0, 2, 0),
                                     c(0, 0, 2),
                                     c(1, 1, 0),
                                     c(1, 0, 1),
                                     c(0, 1, 1)))
#' Constructor for type graph
#' @export
#' @param root_id name of root, which is the integer num_tip - 1
#' @param node_id vector of names of the graph nodes, i.e. progenitor states
#' @param tip_id vector of names of the graph tip, i.e. terminal states
#' @param merge_matrix character matrix specifying the sequence of merges
#' @param differentiation_time commitment time of states
#' @param differntiation_mode_probs commitment bias, need to be vector of length 5, only first values two are non-zero
#' @param double_time list of doubling times for each state
#' @param founder_size founder size
#' @param target_time time of sample collection
make_type_graph <- function(root_id,
                            node_id,
                            tip_id,
                            merge_matrix,
                            differntiation_time,
                            differntiation_mode_probs,
                            double_time,
                            founder_size,
                            target_time,
                            transition_selection = NULL) {
        assertthat::assert_that(length(node_id) == (length(tip_id) - 1))
        assertthat::assert_that(all(root_id %in% node_id))
        all_id = c(node_id, tip_id)
        non_root_id = all_id[all_id != root_id]
        assertthat::assert_that(all(non_root_id %in% as.character(merge_matrix)))
        assertthat::assert_that(all(all_id %in% names(differntiation_time)))
        assertthat::assert_that(all(node_id %in% names(differntiation_mode_probs)))
        assertthat::assert_that(all(all_id %in% names(double_time)))
        out = list(
                root_id = root_id,
                node_id = node_id,
                tip_id = tip_id,
                all_id = c(node_id, tip_id),
                merge = merge_matrix,
                diff_time = differntiation_time,
                diff_mode_probs = differntiation_mode_probs,
                lifetime = double_time,
                founder_size = founder_size,
                trans_sel = transition_selection,
                target_time = target_time
        )
        out = generate_diff_outcome(out)
        out = generate_edge_df(out)
        class(out) = "type_graph"
        return(out)
}
#' Generate differentiation outcome matrix for type graph
#' @export
generate_diff_outcome = function(type_graph) {
        # allowed differentiation outcomes
        merge_mat = type_graph$merge
        type_graph$diff_outcome = lapply(1:nrow(merge_mat), function(parent) {
                out = rbind(
                        c(parent, parent),
                        rep(merge_mat[parent, 1], 2),
                        rep(merge_mat[parent, 2], 2),
                        merge_mat[parent,],
                        c(merge_mat[parent, 1], parent),
                        c(merge_mat[parent, 2], parent)
                )
                out = matrix(as.character(out), nrow(out))
        })
        names(type_graph$diff_outcome) = rownames(merge_mat)
        type_graph
}
#' Generate edge data.frame for type graph, getting edge lengths from laid out graph
#' @export
generate_edge_df <- function(type_graph) {
        diff_edge = do.call(rbind, lapply(1:nrow(type_graph$merge), function(i) {
                data.frame(in_node = i,
                           out_node = type_graph$merge[i,])
        }))
        diff_time_edge = list()
        gens0 = make_gens(type_graph)
        # gens0_time = map_dbl(gens0, function(x) {
        #         ifelse(x$active,
        #                x$end_time - x$double_time,
        #                type_graph$target_time)
        # })
        type_graph$edges = diff_edge
        gens0_time = get_true_diff_time(type_graph)
        diff_edge$length = gens0_time[as.character(diff_edge[, "out_node"])] - gens0_time[as.character(diff_edge[, "in_node"])]
        assertthat::assert_that(all(diff_edge$length > 0))
        # diff_time_edge[[type_graph$root_id]] = gens0_time[[type_graph$root_id]]
        # for (x in c(type_graph$node_id, type_graph$tip_id)) {
        #         if (x != root_id) {
        #                 in_node = diff_edge[diff_edge[, "out_node"] == x, "in_node"]
        #                 diff_time_edge[[x]] = gens0_time[[x]] - type_graph$diff_time[[as.character(in_node)]]
        #                 assertthat::assert_that(diff_time_edge[[x]] > 0)
        #         }
        # }
        # diff_edge$length = unlist(diff_time_edge)[match(diff_edge$out_node, names(diff_time_edge))]
        type_graph$edges = diff_edge
        type_graph$root_time = gens0_time[[type_graph$root_id]]
        type_graph
}
get_true_diff_time <- function(type_graph) {
        gens0 = make_gens(type_graph = type_graph)
        diff_time = map_dbl(c(type_graph$node_id, type_graph$tip_id), function(x) {
                message(x)
                x_gens = gens0[[x]]
                if (x_gens$active) {
                        if(x_gens$num_gen >= 1) {
                                return(x_gens$start_time + (x_gens$num_gen - 1) * x_gens$double_time)
                        } else {
                                assertthat::assert_that(x_gens$num_gen == 0)
                                parent_node = get_parent_node(x, type_graph)
                                x_parent_gens = gens0[[parent_node]]
                                out = x_parent_gens$start_time + x_parent_gens$num_gen * x_parent_gens$double_time
                                return(out)
                        }
                } else {
                        out = x_gens$start_time + x_gens$num_gen * x_gens$double_time
                        return(out)
                }
        })
        names(diff_time) = c(type_graph$node_id, type_graph$tip_id)
        # make diff_time target at tip
        diff_time[type_graph$tip_id] = type_graph$target_time
        diff_time
}

as.phylo <- function(tg) UseMethod("as.phylo")
#' Convert type graph to phylo object
#' @export
as.phylo.type_graph <- function(type_graph) {
        root_id = type_graph$root_id
        diff_edge = type_graph$edges
        # need to redefine the edge length, by subtracting one doubling time
        leaf_newick = lapply(type_graph$tip_id, function(id) {
                paste0(id, ":", diff_edge[diff_edge$out_node == id, "length"])
        })
        eligible_nodes = type_graph$tip_id
        eligible_newick = leaf_newick
        while (length(eligible_nodes) > 1) {
                p0_df = diff_edge[diff_edge$out_node %in% eligible_nodes,]
                p0_count = table(p0_df$in_node)
                p0_merging = as.numeric(names(p0_count)[p0_count == 2])
                p0_df_merge = p0_df[p0_df$in_node %in% p0_merging,]
                d0_newick = eligible_newick[match(p0_df_merge$out_node, eligible_nodes)]
                p0_newick0 = lapply(split(d0_newick, p0_df_merge$in_node), function(x) {
                        paste0("(", x[[1]],  ",", x[[2]] , ")")
                })
                p0_length = sapply(as.numeric(names(p0_newick0)), function(x)
                        diff_edge[diff_edge$out_node == x, "length"])
                p0_newick = mapply(
                        function(z, x, y)
                                paste0(x, z, ":", y),
                        names(p0_newick0),
                        p0_newick0,
                        p0_length
                )

                eligible_nodes1 = c(eligible_nodes[!eligible_nodes %in% p0_df_merge$out_node],
                                    as.numeric(names(p0_newick)))
                eligible_newick1 = c(eligible_newick[!eligible_nodes %in% p0_df_merge$out_node],
                                     p0_newick)
                eligible_nodes = eligible_nodes1
                eligible_newick = eligible_newick1
        }
        tr = ape::read.tree(text = paste0(eligible_newick[[1]],
                                          type_graph$diff_time[[root_id]], ";"))
        tr$root.edge = type_graph$root_time
        tr
}

#' Plot type graph using ggraph
#' @import igraph
#' @import ggraph
#' @export
plot_type_graph <- function(type_graph, node_col_mapper = white_col_mapper,
                            node_label_mapper = NULL, effective_time = T,
                            show_node_text = T) {
        gens0 = make_gens(type_graph)
        type_ig = ape::as.igraph.phylo(as.phylo(type_graph))
        V(type_ig)$name = str_replace_all(V(type_ig)$name, ":", "")
        V(type_ig)$Type = V(type_ig)$name
        V(type_ig)$Time = get_true_diff_time(type_graph)[V(type_ig)$name]
        # node_frac = do.call(rbind,
        #                     lapply(type_graph$node_id, function(id)
        #                             data.frame(
        #                                     id = as.character(type_graph$merge[id,]),
        #                                     fraction = type_graph$diff_mode_probs[[id]][1:2]
        #                             )))
        # V(type_ig)$fraction = node_frac[match(V(type_ig)$name, node_frac$id), "fraction"]
        V(type_ig)$counts = map_dbl(gens0[V(type_ig)$name], function(x) {
                ifelse(x$active,
                       ifelse(effective_time, x$end_count/2, x$end_count),
                       x$end_count)
        })
        double_time = unlist(type_graph$lifetime)
        V(type_ig)$double_time = double_time[match(V(type_ig)$name, names(double_time))]
        # V(type_ig)$label = paste0(
        #         "double time: ",
        #         V(type_ig)$double_time,
        #         "\n",
        #         ifelse(is.na(V(type_ig)$fraction),
        #                "",
        #                paste0(                "fraction: ",
        #                                       V(type_ig)$fraction,
        #                                       "\n")),
        #         "counts: ",
        #         ifelse(V(type_ig)$counts > 1000,
        #                format(V(type_ig)$counts, digits = 2),
        #                V(type_ig)$counts),
        #         "\n")
        V(type_ig)$label = paste0(ifelse(V(type_ig)$counts > 1000,
                       format(V(type_ig)$counts, digits = 2),
                       V(type_ig)$counts))
        if (!is.null(node_label_mapper)) {
                V(type_ig)$node_label = node_label_mapper(V(type_ig)$name)
        } else {
                V(type_ig)$node_label =V(type_ig)$name
        }
        node_col = node_col_mapper(sort(V(type_ig)$name))
        type_ig = type_ig %>% add_vertices(1, name = "Founder",
                                           node_label = "F",
                                           Time = 0) %>%
                add_edges(c("Founder", as.character(type_graph$root_id)))
        V(type_ig)[["Founder"]]$label = paste0("count: ", type_graph$founder_size)
        node_col = c(node_col, "#FFFFFF")
        names(node_col)[length(node_col)] = "Founder"

        for (x_n in names(type_graph$trans_sel)) {
                x = type_graph$trans_sel[[x_n]]
                x_gens = gens0[[x_n]]
                node_out = x_n
                if (x_n == type_graph$root_id) {
                        node_in = "Founder"
                } else {
                        node_in = type_graph$edges$in_node[type_graph$edges$out_node == x_n]
                }
                for (j in 1:length(x)) {
                        y = x[[j]]
                        trans_name = paste0(x_n, " unsampled ", j)
                        type_ig = type_ig - edge(paste0(node_in, "|", node_out))
                        # need to get the exact time, can be tricky
                        trans_time = ifelse(effective_time,
                                            x_gens$start_time + x_gens$double_time * (which(x_gens$cell_mode_counts[, 2] > 0)[j] - 3),
                                            x_gens$start_time + x_gens$double_time * (which(x_gens$cell_mode_counts[, 2] > 0)[j] - 2)
                        )
                        type_ig = type_ig %>% add_vertices(1, name = trans_name,
                                                           Time = trans_time) %>%
                                add_edges(c(node_in, trans_name)) %>%
                                add_edges(c(trans_name, node_out))

                        # V(type_ig)[[trans_name]]$label = paste0("frac: ", y$fraction[1])
                        V(type_ig)[[trans_name]]$label = NA
                        node_col = c(node_col, "#FFFFFF")
                        names(node_col)[length(node_col)] = trans_name
                        node_in  = trans_name
                }
        }
        # things not mapped by label mapper
        V(type_ig)$node_label[is.na(V(type_ig)$node_label)] = V(type_ig)$name[is.na(V(type_ig)$node_label)]
        g_out = ggraph(type_ig, layout = 'dendrogram', height = Time) +
                geom_edge_bend(
                        aes(width = factor(1),
                            fontface = 'plain'),
                        arrow = arrow(
                                type = "closed",
                                length = unit(3, "pt"),
                                angle = 45
                        ),
                        start_cap = ggraph::circle(5, 'pt'),
                        end_cap = ggraph::circle(5, 'pt')
                ) +
                geom_node_label(aes(
                        label = node_label,
                        fill = name,
                        size = factor(1)
                ),
                label.padding = unit(6, 'pt'))
        if (show_node_text) {
                g_out = g_out +
                        geom_node_text(aes(label = label), nudge_x = - 0.2, nudge_y = 0.3)
        }
        g_out = g_out +
                scale_fill_manual(values = node_col) +
                ylim(c(type_graph$target_time, 0)) + ylab('Time') +
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
#' Calculate the counts before next division
#' @export
calculate_type_counts <- function(type_graph, init_count = type_graph$founder_size) {
        states_list_end = list()
        cur_states_list = initial_states_list(type_graph, init_count)
        while (any(sapply(cur_states_list, "[[", "active"))) {
                for (cell_type in names(cur_states_list)) {
                        cells_state = cur_states_list[[cell_type]]
                        if (cells_state$active) {
                                end_states = count_single_type(
                                        cell_type,
                                        cells_state,
                                        type_graph,
                                        differentiate = F
                                )
                                assertthat::assert_that(length(end_states) == 1)
                                if (end_states[[1]]$active) {
                                        states_list_end[[names(end_states)]] = end_states[[1]]
                                }
                        }
                }
                new_states_list = advance_cells_state(cur_states_list, type_graph = type_graph)
                cur_states_list = new_states_list
        }
        type_end_counts = sapply(states_list_end, "[[", "count")
        names(type_end_counts) = names(states_list_end)
        type_end_counts
}
calculate_type_counts_new <- function(type_graph) {
        gens0 = make_gens(type_graph = type_graph)
        return(sapply(gens0, "[[", "end_count"))
}
#' Make initial states list
#' @export
initial_states_list <- function(type_graph, count = 1) {
        states_list = list()
        states_list[[type_graph$root_id]] =  make_cells_state(time = 0,
                                                              count = count,
                                                              active = T)
        states_list
}
#' Change target time
#' @export
change_target_time <- function(type_graph, target_time) {
        type_graph$target_time = target_time
        type_graph
}
check_type_count_seq <- function(count_vec) {
        all(diff(count_vec) >= 0 | count_vec[-1] == 0)
}
backward_merge_sequence <- function(type_graph) {
        node_seq = character()
        cur_nodes = type_graph$tip_id
        while(length(node_seq) < length(type_graph$node_id)) {
                for(node_id in type_graph$node_id) {
                        node_dau = as.character(type_graph$merge[node_id, ])
                        if(all(node_dau  %in% cur_nodes)) {
                                cur_nodes = cur_nodes[!cur_nodes %in% node_dau]
                                cur_nodes = c(cur_nodes, node_id)
                                node_seq = append(node_seq, node_id)
                        }
                }
        }
        return(node_seq)
}
forward_merge_sequence <- function(type_graph) {
        node_seq = character()
        cur_nodes = type_graph$root_id
        while(length(node_seq) < length(type_graph$node_id)) {
                for(node_id in cur_nodes) {
                        if (!node_id %in% type_graph$tip_id) {
                                node_dau = as.character(type_graph$merge[node_id, ])
                                cur_nodes = cur_nodes[cur_nodes != node_id]
                                cur_nodes = c(cur_nodes, node_dau)
                                node_seq = append(node_seq, node_id)
                        }
                }
        }
        return(node_seq)
}
merge2phylo <- function(merge_df) {
        edge_nodes = c(merge_df$merge1, merge_df$merge2)
        tip_id = unique(edge_nodes[as.numeric(edge_nodes) < 0])
        node_height = c(map(merge_df$height, function(x) x),
                        map(tip_id, function(x) 0))
        names(node_height) = c(merge_df$node, tip_id)
        edge_mat = bind_rows(map(1:nrow(merge_df), function(i) {
                tibble(in_node = merge_df$node[i],
                       out_node = c(merge_df$merge1[i], merge_df$merge2[i]),
                       length = c(node_height[[merge_df$node[i]]] - node_height[[merge_df$merge1[i]]],
                                  node_height[[merge_df$node[i]]] - node_height[[merge_df$merge2[i]]]))
        }))
        leaf_newick = lapply(tip_id, function(id) {
                paste0(id, ":", edge_mat[edge_mat$out_node == id, "length"])
        })
        eligible_nodes = tip_id
        eligible_newick = leaf_newick
        while (length(eligible_nodes) > 1) {
                p0_df = edge_mat[edge_mat$out_node %in% eligible_nodes,]
                p0_count = table(p0_df$in_node)
                p0_merging = as.numeric(names(p0_count)[p0_count == 2])
                p0_df_merge = p0_df[p0_df$in_node %in% p0_merging,]
                d0_newick = eligible_newick[match(p0_df_merge$out_node, eligible_nodes)]
                p0_newick0 = lapply(split(d0_newick, p0_df_merge$in_node), function(x) {
                        paste0("(", x[[1]],  ",", x[[2]] , ")")
                })
                p0_length = sapply(as.numeric(names(p0_newick0)), function(x)
                        edge_mat[edge_mat$out_node == x, "length"])
                p0_newick = mapply(
                        function(z, x, y)
                                paste0(x, z, ":", y),
                        names(p0_newick0),
                        p0_newick0,
                        p0_length
                )
                eligible_nodes1 = c(eligible_nodes[!eligible_nodes %in% p0_df_merge$out_node],
                                    as.numeric(names(p0_newick)))
                eligible_newick1 = c(eligible_newick[!eligible_nodes %in% p0_df_merge$out_node],
                                     p0_newick)
                eligible_nodes = eligible_nodes1
                eligible_newick = eligible_newick1
        }
        tr = ape::read.tree(text = paste0(eligible_newick[[1]],
                                          0, ";"))
        tr
}

plot_type_graph_clean <- function(type_graph, node_col_mapper = white_col_mapper,
                            node_label_mapper = NULL, effective_time = T,
                            show_node_text = F, node_size = 3) {
        gens0 = make_gens(type_graph)
        type_ig = ape::as.igraph.phylo(as.phylo(type_graph))
        V(type_ig)$name = str_replace_all(V(type_ig)$name, ":", "")
        V(type_ig)$Type = V(type_ig)$name
        V(type_ig)$Time = get_true_diff_time(type_graph)[V(type_ig)$name]
        # node_frac = do.call(rbind,
        #                     lapply(type_graph$node_id, function(id)
        #                             data.frame(
        #                                     id = as.character(type_graph$merge[id,]),
        #                                     fraction = type_graph$diff_mode_probs[[id]][1:2]
        #                             )))
        # V(type_ig)$fraction = node_frac[match(V(type_ig)$name, node_frac$id), "fraction"]
        V(type_ig)$counts = map_dbl(gens0[V(type_ig)$name], function(x) {
                ifelse(x$active,
                       ifelse(effective_time, x$end_count/2, x$end_count),
                       x$end_count)
        })
        double_time = unlist(type_graph$lifetime)
        V(type_ig)$double_time = double_time[match(V(type_ig)$name, names(double_time))]
        # V(type_ig)$label = paste0(
        #         "double time: ",
        #         V(type_ig)$double_time,
        #         "\n",
        #         ifelse(is.na(V(type_ig)$fraction),
        #                "",
        #                paste0(                "fraction: ",
        #                                       V(type_ig)$fraction,
        #                                       "\n")),
        #         "counts: ",
        #         ifelse(V(type_ig)$counts > 1000,
        #                format(V(type_ig)$counts, digits = 2),
        #                V(type_ig)$counts),
        #         "\n")
        V(type_ig)$label = paste0(ifelse(V(type_ig)$counts > 1000,
                                         format(V(type_ig)$counts, digits = 2),
                                         V(type_ig)$counts))
        if (!is.null(node_label_mapper)) {
                V(type_ig)$node_label = node_label_mapper(V(type_ig)$name)
        } else {
                V(type_ig)$node_label =V(type_ig)$name
        }
        node_col = node_col_mapper(sort(V(type_ig)$name))
        # type_ig = type_ig %>% add_vertices(1, name = "Founder",
        #                                    node_label = "F",
        #                                    Time = 0) %>%
        #         add_edges(c("Founder", as.character(type_graph$root_id)))
        # V(type_ig)[["Founder"]]$label = paste0("count: ", type_graph$founder_size)
        # node_col = c(node_col, "#FFFFFF")
        # names(node_col)[length(node_col)] = "Founder"

        for (x_n in names(type_graph$trans_sel)) {
                x = type_graph$trans_sel[[x_n]]
                x_gens = gens0[[x_n]]
                node_out = x_n
                if (x_n == type_graph$root_id) {
                        node_in = "Founder"
                } else {
                        node_in = type_graph$edges$in_node[type_graph$edges$out_node == x_n]
                }
                for (j in 1:length(x)) {
                        y = x[[j]]
                        trans_name = paste0(x_n, " unsampled ", j)
                        type_ig = type_ig - edge(paste0(node_in, "|", node_out))
                        # need to get the exact time, can be tricky
                        trans_time = ifelse(effective_time,
                                            x_gens$start_time + x_gens$double_time * (which(x_gens$cell_mode_counts[, 2] > 0)[j] - 3),
                                            x_gens$start_time + x_gens$double_time * (which(x_gens$cell_mode_counts[, 2] > 0)[j] - 2)
                        )
                        type_ig = type_ig %>% add_vertices(1, name = trans_name,
                                                           Time = trans_time) %>%
                                add_edges(c(node_in, trans_name)) %>%
                                add_edges(c(trans_name, node_out))

                        # V(type_ig)[[trans_name]]$label = paste0("frac: ", y$fraction[1])
                        V(type_ig)[[trans_name]]$label = NA
                        node_col = c(node_col, "#FFFFFF")
                        names(node_col)[length(node_col)] = trans_name
                        node_in  = trans_name
                }
        }
        # things not mapped by label mapper
        V(type_ig)$node_label[is.na(V(type_ig)$node_label)] = V(type_ig)$name[is.na(V(type_ig)$node_label)]
        g_out = ggraph(type_ig, layout = 'dendrogram', height = Time) +
                geom_edge_bend(
                        aes(width = factor(1),
                            fontface = 'plain'),
                        arrow = arrow(
                                type = "closed",
                                length = unit(2, "pt"),
                                angle = 45
                        ),
                        start_cap = ggraph::circle(5, 'pt'),
                        end_cap = ggraph::circle(5, 'pt')
                ) +
                geom_node_point(aes(color = Type), shape = 15, size = node_size)
        if (show_node_text) {
                g_out = g_out +
                        geom_node_text(aes(label = label), nudge_x = - 0.2, nudge_y = 0.3)
        }
        g_out = g_out +
                scale_color_manual(values = node_col) +
                ylim(c(type_graph$target_time, 0)) + ylab('Time') +
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
plot_type_graph_dendro <- function(type_graph, node_col_mapper = white_col_mapper,
                                  node_label_mapper = NULL, effective_time = T,
                                  show_node_text = F, node_size = 3) {
        gens0 = make_gens(type_graph)
        type_ig = ape::as.igraph.phylo(as.phylo(type_graph))
        V(type_ig)$name = str_replace_all(V(type_ig)$name, ":", "")
        V(type_ig)$Type = V(type_ig)$name
        V(type_ig)$Time = get_true_diff_time(type_graph)[V(type_ig)$name]
        # node_frac = do.call(rbind,
        #                     lapply(type_graph$node_id, function(id)
        #                             data.frame(
        #                                     id = as.character(type_graph$merge[id,]),
        #                                     fraction = type_graph$diff_mode_probs[[id]][1:2]
        #                             )))
        # V(type_ig)$fraction = node_frac[match(V(type_ig)$name, node_frac$id), "fraction"]
        V(type_ig)$counts = map_dbl(gens0[V(type_ig)$name], function(x) {
                ifelse(x$active,
                       ifelse(effective_time, x$end_count/2, x$end_count),
                       x$end_count)
        })
        double_time = unlist(type_graph$lifetime)
        V(type_ig)$double_time = double_time[match(V(type_ig)$name, names(double_time))]
        # V(type_ig)$label = paste0(
        #         "double time: ",
        #         V(type_ig)$double_time,
        #         "\n",
        #         ifelse(is.na(V(type_ig)$fraction),
        #                "",
        #                paste0(                "fraction: ",
        #                                       V(type_ig)$fraction,
        #                                       "\n")),
        #         "counts: ",
        #         ifelse(V(type_ig)$counts > 1000,
        #                format(V(type_ig)$counts, digits = 2),
        #                V(type_ig)$counts),
        #         "\n")
        V(type_ig)$label = paste0(ifelse(V(type_ig)$counts > 1000,
                                         format(V(type_ig)$counts, digits = 2),
                                         V(type_ig)$counts))
        if (!is.null(node_label_mapper)) {
                V(type_ig)$node_label = node_label_mapper(V(type_ig)$name)
        } else {
                V(type_ig)$node_label =V(type_ig)$name
        }
        node_col = node_col_mapper(sort(V(type_ig)$name))
        # type_ig = type_ig %>% add_vertices(1, name = "Founder",
        #                                    node_label = "F",
        #                                    Time = 0) %>%
        #         add_edges(c("Founder", as.character(type_graph$root_id)))
        # V(type_ig)[["Founder"]]$label = paste0("count: ", type_graph$founder_size)
        # node_col = c(node_col, "#FFFFFF")
        # names(node_col)[length(node_col)] = "Founder"

        for (x_n in names(type_graph$trans_sel)) {
                x = type_graph$trans_sel[[x_n]]
                x_gens = gens0[[x_n]]
                node_out = x_n
                if (x_n == type_graph$root_id) {
                        node_in = "Founder"
                } else {
                        node_in = type_graph$edges$in_node[type_graph$edges$out_node == x_n]
                }
                for (j in 1:length(x)) {
                        y = x[[j]]
                        trans_name = paste0(x_n, " unsampled ", j)
                        type_ig = type_ig - edge(paste0(node_in, "|", node_out))
                        # need to get the exact time, can be tricky
                        trans_time = ifelse(effective_time,
                                            x_gens$start_time + x_gens$double_time * (which(x_gens$cell_mode_counts[, 2] > 0)[j] - 3),
                                            x_gens$start_time + x_gens$double_time * (which(x_gens$cell_mode_counts[, 2] > 0)[j] - 2)
                        )
                        type_ig = type_ig %>% add_vertices(1, name = trans_name,
                                                           Time = trans_time) %>%
                                add_edges(c(node_in, trans_name)) %>%
                                add_edges(c(trans_name, node_out))

                        # V(type_ig)[[trans_name]]$label = paste0("frac: ", y$fraction[1])
                        V(type_ig)[[trans_name]]$label = NA
                        node_col = c(node_col, "#FFFFFF")
                        names(node_col)[length(node_col)] = trans_name
                        node_in  = trans_name
                }
        }
        # things not mapped by label mapper
        V(type_ig)$node_label[is.na(V(type_ig)$node_label)] = V(type_ig)$name[is.na(V(type_ig)$node_label)]
        g_out = ggraph(type_ig, layout = 'dendrogram', height = Time) +
                geom_edge_elbow()
        if (show_node_text) {
                g_out = g_out +
                        geom_node_text(aes(label = label), nudge_x = - 0.2, nudge_y = 0.3)
        }
        g_out = g_out +
                scale_color_manual(values = node_col) +
                ylim(c(type_graph$target_time, 0)) + ylab('Time')
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

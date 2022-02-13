#' Convert cells to edge data.frame
cells_to_edges <- function(cells) {
        if (length(cells) == 1) {
                stopifnot(is.na(cells))
                return(NULL)
        } else {
                return(data.frame(in_node = cells$parent,
                                  out_node = cells$id))
        }
}
#' Convert cells to node data.frame
cells_to_nodes <- function(cells) {
        if (length(cells) == 1) {
                stopifnot(is.na(cells))
                return(NULL)
        } else {
                return(data.frame(id = cells$id,
                                  birth_time = cells$birth_time,
                                  type = cells$type_state,
                                  double_time = cells$life_duration,
                                  sample_size = cells$sample_size))
        }
}
#' Convert sim_history to edges data.frame
sim_history_to_edges <- function(sim_history) {
        edges_df = do.call(rbind, lapply(sim_history[2:length(sim_history)], function(x) {
                rbind(cells_to_edges(x$active),
                      cells_to_edges(x$inactive))
        }))
        edges_df = edges_df[!duplicated(edges_df), ]
        edges_df
}
#' Convert sim_history to nodes data.frame
sim_history_to_nodes <- function(sim_history) {
        nodes_df = do.call(rbind, lapply(sim_history[1:length(sim_history)], function(x) {
                rbind(cells_to_nodes(x$active),
                      cells_to_nodes(x$inactive))
        }))
        nodes_df = nodes_df[!duplicated(nodes_df), ]
        edges_df = sim_history_to_edges(sim_history)

        # Identify leaf and internal nodes
        all_nodes = nodes_df$id
        leaf_nodes = sim_history[[length(sim_history)]]$inactive$id
        internal_nodes = all_nodes[!all_nodes %in% leaf_nodes]

        nodes_df$leaf = nodes_df$id %in% leaf_nodes
        internal_nodes_deg = table(edges_df[, 1])
        nodes_df$degree = 0
        nodes_df$degree[nodes_df$leaf] = 1
        nodes_df$degree[!nodes_df$leaf] = internal_nodes_deg[match(nodes_df$id[!nodes_df$leaf],
                                                                   as.numeric(names(internal_nodes_deg)))]
        nodes_df
}
#' Add edge length to edge df
add_edge_length <- function(edges_df, nodes_df) {
        edges_df$length = nodes_df[match(edges_df$out_node, nodes_df$id), "double_time"]
        edges_df
}
#' Remove internal nodes with single daughter
simplify_edges_df <- function(edges_df, nodes_df) {
        edges_df = add_edge_length(edges_df, nodes_df)
        single_internal_node = nodes_df$id[(!nodes_df$leaf) & (nodes_df$degree == 2)]
        # degree 1 or 2
        for (x in single_internal_node) {
                in_edge = which(edges_df$out_node == x)
                out_edge = which(edges_df$in_node == x)
                len_in = edges_df[in_edge, "length"]
                len_out = edges_df[out_edge, "length"]
                # modify
                edges_df[in_edge, "length"] = len_in + len_out
                edges_df[in_edge, "out_node"] = edges_df[out_edge, "out_node"]
                edges_df = edges_df[-out_edge, ]
        }
        return(edges_df)
}
#' Convert simplified edges df to phylo
edges_df_to_newick <- function(edges_df_simple, nodes_df) {
        assertthat::assert_that(!is.null(edges_df_simple$length))
        nodes_leaf = nodes_df[nodes_df$leaf, ]
        leaf_newick = paste0(nodes_leaf[, "id"], ":", edges_df_simple[match(nodes_leaf$id, edges_df_simple$out_node), "length"])
        eligible_nodes = nodes_leaf$id
        eligible_newick = leaf_newick
        eligible_total_old = length(eligible_nodes)
        while(length(eligible_nodes)> 1){
                p0_df = edges_df_simple[edges_df_simple$out_node %in% eligible_nodes, ]
                p0_count = table(p0_df$in_node)
                p0_merging = names(p0_count)[p0_count == 2]
                p0_df_merge = p0_df[p0_df$in_node %in% p0_merging, ]
                d0_newick = eligible_newick[match(p0_df_merge$out_node, eligible_nodes)]
                p0_newick0 = lapply(split(d0_newick, factor(p0_df_merge$in_node)), function(x) {
                        paste0("(", x[[1]],  ",", x[[2]] , ")")
                })
                # names(p0_newick0) = levels(p0_df_merge$in_node)[sort(unique(p0_df_merge$in_node))]
                p0_length = sapply(names(p0_newick0), function(x) edges_df_simple[edges_df_simple$out_node == x, "length"])
                p0_newick = mapply(function(z, x, y) paste0(x, z, ":", y), names(p0_newick0), p0_newick0, p0_length)
                eligible_nodes1 = c(as.character(eligible_nodes[!eligible_nodes %in% p0_df_merge$out_node]),
                                    names(p0_newick))
                eligible_newick1 = c(eligible_newick[!eligible_nodes %in% p0_df_merge$out_node],
                                     p0_newick)
                eligible_nodes = eligible_nodes1
                eligible_newick = eligible_newick1
                eligible_total = length(eligible_nodes)
                assertthat::assert_that(eligible_total < eligible_total_old)
                eligible_total_old = eligible_total
                print(eligible_total)
        }
        assertthat::assert_that(length(eligible_newick) == 1)
        eligible_newick[[1]]
}
#' Get sampled cells of at target time
#' @export
get_sampled_cells <- function(sim_history) {
        sim_history[[length(sim_history)]]$inactive
}
#' Get sample size of sim_history
#' @export
get_sample_size <- function(sim_history) {
        num_cells(get_sampled_cells(sim_history))
}
#' Get target time from sim_history
get_target_time <- function(sim_history) {
        sampled_cells = get_sampled_cells(sim_history)
        tt = sampled_cells$birth_time + sampled_cells$life_duration
        assertthat::assert_that(is_zero_range(tt))
        tt[1]
}
#' Construct lineage tree from sim_history
#' @export
# construct_lineage <- function(sim_history) {
#         edges_df = sim_history_to_edges(sim_history)
#         nodes_df = sim_history_to_nodes(sim_history)
#         assertthat::assert_that(nrow(nodes_df) == (nrow(edges_df)+1))
#         assertthat::assert_that(sum(nodes_df$degree) == (nrow(edges_df) + get_sample_size(sim_history)))
#
#         edges_df_simple = simplify_edges_df(edges_df, nodes_df)
#         assertthat::assert_that(nrow(edges_df_simple) == 2 * (get_sample_size(sim_history) - 1))
#         tr = ape::read.tree(text=paste0(edges_df_to_newick(edges_df_simple), nodes_df[nodes_df$id == 0, 'double_time'], ";"))
#         list(phylo = tr,
#              target_time = get_target_time(sim_history),
#              sample_size = get_sample_size(sim_history),
#              nodes = nodes_df,
#              edges = edges_df,
#              edges_simple = edges_df_simple)
# }
#' Plot lineage
#' @export
# plot_lineage <- function(lin_obj) {
#         g = as.igraph(lin_obj$phylo)
#         nodes_match = lin_obj$nodes[match(V(g)$name, as.character(lin_obj$nodes$id)), ]
#         vertex_attr(g, "Time") <- nodes_match$birth_time + nodes_match$double_time
#         vertex_attr(g, "Type") <- as.character(nodes_match$type)
#         ggraph(g, layout='dendrogram', height = Time) +
#                 geom_edge_diagonal() +
#                 geom_node_point(aes(col=Type)) +
#                 theme_void() + scale_color_manual(values = type_col_mapper(sort(levels(nodes_match$type)))) +
#                 ylim(c(lin_obj$target_time, 0)) + ylab('Time') +
#                 theme(axis.line.y = element_line(),
#                       axis.ticks.y = element_line(),
#                       axis.text.y = element_text(),
#                       axis.title = element_text(size = 12),
#                       panel.background = element_blank(),
#                       axis.title.x = element_blank())
# }
normalize_tree <- function(tree, check.ultrametric=TRUE){

        if(check.ultrametric){
                if(!is.ultrametric(tree))
                        stop("the input tree is not ultrametric")
        }

        nTips  = length(tree$tip.label)
        rNode  = nTips + 1
        nEdges = Nedge(tree)

        g        = graph.edgelist(tree$edge, directed = TRUE)
        root2tip = get.shortest.paths(g, rNode, to=1:nTips, mode="out", output="epath")$epath

        root.edge <- ifelse(is.null(tree$root.edge), 0, tree$root.edge)
        Tval     = root.edge + sum(tree$edge.length[root2tip[[1]] ])
        #Tval = mean ( sapply( 1:nTips, FUN=function(x) sum(tree$edge.length[root2tip[[x]]])   )  )

        tree$edge.length = tree$edge.length / Tval
        if(!is.null(tree$root.edge)){
                tree$root.edge <- tree$root.edge / Tval
        }
        return(tree)
}



get_gr_progenitors <- function(gr, tr_r, total_time) {
        gr_tips = list_dd_and_tips(gr)$tips
        gr_depth = node.depth.edgelength(gr)
        names(gr_depth) = c(gr$tip.label, gr$node.label)
        gr_root_edgelen = total_time - max(gr_depth)
        gr_depth = gr_depth + gr_root_edgelen
        gr_node_field = map(gr$node.label, function(x) {
                get_progenitor_nodes(tr_r, gr_depth[x], gr_tips[[x]])
        })
        names(gr_node_field) = gr$node.label
        gr_node_field
}
get_parent_node <- function(node, type_graph) {
        as.character(type_graph$edges[type_graph$edges$out_node == node, "in_node"])
}
get_parent_split_index <- function(node, parent_node, type_graph) {
        which(as.character(type_graph$merge[parent_node, ]) == node)
}
mean_top_half <- function(vec) {
        if (length(vec) < 2) {
                return(vec)
        } else {
                vec_top = vec[order(vec, decreasing = T)[1:ceiling(length(vec)/2)]]
                return(mean(vec_top))
        }
}
evalute_gr_topo <- function(gr, type_graph) {
        list(rf = phangorn::RF.dist(gr, as.phylo(type_graph)),
             kc0 = treespace::treeDist(gr, as.phylo(type_graph), lambda = 0))
}
evalute_gr_temp <- function(gr, type_graph, total_time) {
        gr_vec_func <- function(gr, total_time) {
                gr_vec = treespace::treeVec(gr, lambda = 1.)
                gr_root = total_time - max(node.depth.edgelength(gr))
                gr_vec + gr_root
        }
        kc1_func <- function(gr, type_graph, total_time) {
                ref_vec = treespace::treeVec(as.phylo(type_graph), lambda = 1.) + type_graph$root_time
                eval_indices = 1:(length(ref_vec)-length(type_graph$tip_id))
                sqrt(mean((gr_vec_func(gr, total_time)[eval_indices] - ref_vec[eval_indices])^2))
        }
        list(rf = phangorn::RF.dist(gr, as.phylo(type_graph)),
             kc0 = treespace::treeDist(gr, as.phylo(type_graph), lambda = 0),
             kc1 = kc1_func(gr, type_graph, total_time))
}
get_progenitor_nodes <- function(tr_r, gr_node_time, gr_tip_set, type_func = get_type_from_id) {
        tr_r = name_nodes(tr_r)
        out_list = list_dd_and_tips(tr_r)
        tr_dd = out_list$dd
        tr_tip_list = out_list$tips

        edge_df = make_edge_df(tr_r)
        node_depth = node.depth.edgelength(tr_r)
        names(node_depth) = c(tr_r$tip.label, tr_r$node.label)
        edge_df$from_depth = node_depth[edge_df$from]
        edge_df$to_depth = node_depth[edge_df$to]

        edge_sel = filter(edge_df, from_depth <= gr_node_time, to_depth > gr_node_time)
        progenitor_nodes = unique(edge_sel$from)
        progenitor_nodes = progenitor_nodes[map_lgl(tr_tip_list[progenitor_nodes], function(x) any(gr_tip_set %in% type_func(x)))]
        progenitor_nodes
}

make_edge_df <- function(tr_r) {
        # get edge_df with node names
        edge_mapper = character(max(tr_r$edge))
        edge_mapper[(1:Nnode(tr_r))+Ntip(tr_r)] = tr_r$node.label
        edge_mapper[1:Ntip(tr_r)] = tr_r$tip.label
        edge_df = tibble(from = edge_mapper[tr_r$edge[, 1]],
                         to = edge_mapper[tr_r$edge[, 2]])
        edge_df$length = tr_r$edge.length
        edge_df
}

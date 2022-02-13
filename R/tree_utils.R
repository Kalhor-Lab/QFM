name_nodes <- function(tr) {
        if(!is.null(tr$node.label)) {
                return(tr)
        }
        tr$node.label = paste0("Node-", 1:tr$Nnode)
        tr
}
evalute_gr <- function(gr, ref) {
        list(rf = phangorn::RF.dist(gr, ref),
             kc0 = treespace::treeDist(gr, ref, lambda = 0),
             kc1 = treespace::treeDist(gr, ref, lambda = 1),
             treevec_cor = cor(treespace::treeVec(gr),
                               treespace::treeVec(ref), method = "spearman"))
}
get_node_time <- function(tr, node_labels) {
        node_depth = node.depth.edgelength(tr)
        names(node_depth) = c(tr$tip.label, tr$node.label)
        node_depth[node_labels]
}
get_mrca_graph <- function(type_graph, types_vec) {
        gr = as.phylo(type_graph)
        gr$node.label[getMRCA(gr, types_vec) - Ntip(gr)]
}
list_dd_and_tips <- function(tr_r) {
        # get edge_df with node names
        edge_mapper = character(max(tr_r$edge))
        edge_mapper[(1:Nnode(tr_r))+Ntip(tr_r)] = tr_r$node.label
        edge_mapper[1:Ntip(tr_r)] = tr_r$tip.label
        edge_df = tibble(from = edge_mapper[tr_r$edge[, 1]],
                         to = edge_mapper[tr_r$edge[, 2]])
        edge_df$length = tr_r$edge.length
        tr_dd_new = split(edge_df$to, edge_df$from)
        tr_dd_new = tr_dd_new[tr_r$node.label]
        # identical(tr_dd_new, tr_dd)

        m_seq = make_merge_sequences(dplyr::rename(edge_df, out_node = to, in_node = from),
                                     tip_id = tr_r$tip.label)
        # build node_tip_lists with merges
        tr_tip_list_new = map(tr_r$tip.label, function(x) x)
        names(tr_tip_list_new) = tr_r$tip.label
        for (x in m_seq) {
                y = map(x, function(dd) {
                        purrr::reduce(tr_tip_list_new[dd], c)
                })
                tr_tip_list_new = append(tr_tip_list_new, y)
        }
        tr_tip_list_new = tr_tip_list_new[tr_r$node.label]
        # identical(tr_tip_list_new, tr_tip_list)
        list(tips = tr_tip_list_new,
             dd = tr_dd_new)
}

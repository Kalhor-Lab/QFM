get_barcodes <- function(gens, tip_id) {
        do.call(rbind, lapply(tip_id, function(tip) {
                gens[[tip]]$end_barcodes
        }))
}
remove_allele_ver <- function(x) {
        x_rm = matrix(sapply(strsplit(x, "\\("), "[[", 1), nrow(x))
        rownames(x_rm) = rownames(x)
        x_rm
}
barcode_allele_onehot <- function(x) {
        # NA essentially treated as zero
        onehot_list = lapply(1:ncol(x), function(i) {
                y = x[, i]
                y[is.na(y)] = "0"
                levels = sort(unique(y[y != "0"]))
                if (length(levels) == 0) {
                        return(NULL)
                }
                1 * do.call(cbind, lapply(levels, function(z) y == z))
        })
        out = do.call(cbind, onehot_list)
        if (!is.null(out)) {
                rownames(out) = rownames(x)
        }
        out
}
barcode_allele_onehot_new <- function(x, include_unmutated = T) {
        # NA essentially treated as zero
        onehot_list = lapply(1:ncol(x), function(i) {
                y = x[, i]
                y[is.na(y)] = "0"
                if (all(y == "0")) {
                        return(NULL)
                }
                levels = sort(unique(y))
                out_mat = 1 * do.call(cbind, lapply(levels, function(z) y == z))
                colnames(out_mat) = levels
                if (!include_unmutated) {
                        out_mat = out_mat[, colnames(out_mat) != "0", drop = F]
                }
                out_mat
        })
        out = do.call(cbind, onehot_list)
        if (!is.null(out)) {
                rownames(out) = rownames(x)
        }
        out
}
dropout_sc_random <- function(sc_mat, prob_vec) {
        # prob_vec provides element specific dropout prob
        out = purrr::reduce(purrr::map(1:ncol(sc_mat), function(j) {
                val = sc_mat[, j]
                val[sample(c(T, F), size = length(val), prob = c(prob_vec[j], 1 - prob_vec[j]), replace = T)] = NA
                val
        }), cbind)
        out
}
filter_allele_onehot <- function(x, lo_cut = 0.0, hi_cut = 1.) {
        x = x[, colSums(x) > 1]
        x = x[, colSums(x) > floor(nrow(x)*lo_cut), drop = F]
        x = x[, colSums(x) < ceiling(nrow(x)*hi_cut), drop = F]
        x
}
get_type_from_id <- function(id_name) {
        out = map_chr(strsplit(id_name, "_"), function(x) {
                ifelse(length(x) == 5, x[2], "Unknown")
        })
        names(out) = id_name
        out
}
get_time_from_id <- function(id_name, type_graph) {
        gens = make_gens(type_graph)
        id_name_sp = strsplit(id_name, "_")
        map_dbl(id_name_sp, function(x)
                gens[[x[2]]]$start_time + as.numeric(x[[4]]) * gens[[x[2]]]$double_time)
}
sps_ver1 <- function(lin_tr, tip_id, type_func = get_type_from_id) {
        if (is.null(lin_tr)) {
                return(NULL)
        }
        lin_tr$edge.length[lin_tr$edge.length < 0] = 0
        node_tips_list = lapply(adephylo::listTips(lin_tr), names)
        if (lin_tr$Nnode == 1) {
                return(NULL)
        }
        node_type_list = lapply(node_tips_list[2:length(node_tips_list)], function(x) table(type_func(x))[tip_id])
        sim_mat = matrix(0,
                         nrow = length(tip_id),
                         ncol = length(tip_id))
        rownames(sim_mat) = colnames(sim_mat) = tip_id
        for (x in node_type_list) {
                x = x[!is.na(x)]
                p_score = 1/2^(sum(x > 0)-1)
                # p_score = 1/sum(x > 0)
                if (sum(x>0) > 1 & sum(x>0) < length(tip_id)) {
                        all_comb = combn(names(x)[x>0], m = 2)
                        for (i in 1:ncol(all_comb)) {
                                sim_mat[all_comb[1, i], all_comb[2, i]] = sim_mat[all_comb[1, i], all_comb[2, i]] + p_score
                                sim_mat[all_comb[2, i], all_comb[1, i]] = sim_mat[all_comb[2, i], all_comb[1, i]] + p_score
                        }
                }
                if (sum(x>0) == 1) {
                        rr = names(x)[x>0]
                        sim_mat[rr, rr] = sim_mat[rr, rr] + 1
                }
                # else {
                #         y = names(x)[x>0]
                #         sim_mat[y, y] = sim_mat[y, y] +  2 * p_score
                # }
        }
        return(sim_mat)
}
# reconstruction
shared_progenitor_score <- function(lin_tr, tip_id, type_func = get_type_from_id) {
        if (is.null(lin_tr)) {
                return(NULL)
        }
        lin_tr$edge.length[lin_tr$edge.length < 0] = 0
        node_tips_list = list_dd_and_tips(lin_tr)$tips
        # node_tips_list = lapply(adephylo::listTips(lin_tr), names)
        if (lin_tr$Nnode == 1) {
                return(NULL)
        }
        node_type_list = lapply(node_tips_list[2:length(node_tips_list)], function(x) table(type_func(x))[tip_id])
        sim_mat = matrix(0,
                         nrow = length(tip_id),
                         ncol = length(tip_id))
        rownames(sim_mat) = colnames(sim_mat) = tip_id
        for (x in node_type_list) {
                x = x[!is.na(x)]
                p_score = 1/2^(sum(x > 0)-1)
                # p_score = 1/sum(x > 0)
                if (sum(x>0) > 1 & sum(x>0) < length(tip_id)) {
                        all_comb = combn(names(x)[x>0], m = 2)
                        for (i in 1:ncol(all_comb)) {
                                sim_mat[all_comb[1, i], all_comb[2, i]] = sim_mat[all_comb[1, i], all_comb[2, i]] + p_score
                                sim_mat[all_comb[2, i], all_comb[1, i]] = sim_mat[all_comb[2, i], all_comb[1, i]] + p_score
                        }
                }
                # else {
                #         y = names(x)[x>0]
                #         sim_mat[y, y] = sim_mat[y, y] +  2 * p_score
                # }
        }
        return(sim_mat)
}
reconstruct_lineage <- function(x, filter_lo_cut = 0.0, tr_func = phangorn::upgma) {
        x_onehot = barcode_allele_onehot_new(x)
        if (is.null(x_onehot)) {
                return(NULL)
        }
        x_onehot = filter_allele_onehot(x_onehot, lo_cut = filter_lo_cut)
        # x_onehot = x_onehot[rowSums(x_onehot) > 0, ]
        if (ncol(x_onehot) == 0) {
                return(NULL)
        }
        shuffle_index = sample(nrow(x_onehot), replace = F)
        lin_tr = tr_func(dist(x_onehot[shuffle_index, ], method = "manhattan"))
        lin_tr = di2multi(lin_tr)
        return(lin_tr)
}
produce_distance_from_lin <- function(lin_tr, tip_id) {
        if (is.null(lin_tr)) {
                return(NULL)
        }
        if (Nnode(lin_tr) == 1) {
                return(NULL)
        }
        sim_mat = shared_progenitor_score(lin_tr, tip_id)
        if (max(sim_mat) == 0 ) {
                dmat = NULL
        } else {
                dmat = 1 - sim_mat/max(sim_mat)
        }
        return(dmat)
}
sc_pipeline <- function(x, tip_id) {
        lin_tr = reconstruct_lineage(x)
        sim_mat = shared_progenitor_score(lin_tr, tip_id)
        dmat = 1 - sim_mat/max(sim_mat)
        return(dmat)
}
#' extend out node to the next bifurcation in the lineage tree
get_next <- function(y, edges_list) {
        l = y$length
        n = as.character(y$out_node)
        # previous sampled internal node encounters leaf
        if (is.null(edges_list[[n]])) {
                return(list(node = as.character(n),
                            length = l,
                            tip = T))
        }
        while(nrow(edges_list[[n]]) < 2) {
                old_n = n
                n = as.character(edges_list[[n]]$out_node)
                # singleton encounters leaf
                if (is.null(edges_list[[n]])) {
                        return(list(node = as.character(n),
                                    length = l,
                                    tip = TRUE))
                } else {
                        l = l + edges_list[[old_n]]$length
                }
        }
        return(list(node = as.character(n),
                    length = l,
                    tip = FALSE))
}
construct_true_lineage <- function(gens1, tip_id) {
        nodes_df = collect_nodes(gens1, tip_id)
        edges_df = collect_edges(gens1)

        edges_df = add_edge_length(edges_df, nodes_df)
        edges_list = split(as.data.table(edges_df), edges_df$in_node)
        # each element is a progenitor
        x_cell = edges_df[1, "in_node"]
        while (nrow(edges_list[[x_cell]]) == 1) {
                x_cell = edges_list[[x_cell]]$out_node
        }
        sp0 = list(edges_list[[x_cell]][1, ],
                   edges_list[[x_cell]][2, ])
        n_tip = 0
        node_count = 0
        edges_simple_collect = list()
        while(TRUE) {
                # p0 = bplapply(sp0, get_next, edges_list = edges_list, BPPARAM = bp_param)
                p0 = lapply(sp0, get_next, edges_list = edges_list)
                # for (x in sp0) {
                #         print(x)
                #         get_next(x, edges_list = edges_list)
                # }
                # table(sapply(p0, "[[", "tip"))
                # as.character(sapply(p0, "[[", "node"))
                # table(sapply(p0, "[[", "length"))

                # p0 = clusterApply(cl, sp0, get_next, edges_list = edges_list)
                edges_simple_collect = append(edges_simple_collect,
                                              list(do.call(rbind, mapply(function(x, y) c(in_node = as.character(x$in_node),
                                                                                          out_node = y$node,
                                                                                          length = y$length),
                                                                         sp0, p0, SIMPLIFY = F))))
                p0_node = p0[!sapply(p0, "[[", "tip")]
                tip_count = sum(sapply(p0, "[[", "tip"))
                node_count = length(p0_node)
                n_tip = n_tip + tip_count
                message(paste0("merging_node: ", node_count))
                if (length(p0_node) > 0) {
                        sp1 = do.call(c, lapply(edges_list[as.character(sapply(p0_node, "[[", "node"))],
                                                function(x) list(x[1, ], x[2, ])))
                        sp0 = sp1
                } else {
                        assertthat::assert_that(sum(nodes_df$leaf) == n_tip)
                        break()
                }
        }
        edges_df_simple = data.frame(do.call(rbind, edges_simple_collect))
        edges_df_simple$length = as.numeric(as.character(edges_df_simple$length))
        edges_df_simple$in_node = as.character(edges_df_simple$in_node)
        edges_df_simple$out_node = as.character(edges_df_simple$out_node)

        tr_nw = edges_df_to_newick(edges_df_simple, nodes_df)
        tr_nw = paste0(tr_nw, nodes_df[1, 'double_time'], ";")
        tr = ape::read.tree(text=tr_nw)
        tr
}

construct_true_lineage_multi <- function(gens1, tip_id) {
        nodes_df = collect_nodes(gens1, tip_id)
        edges_df = collect_edges(gens1)
        edges_df = add_edge_length(edges_df, nodes_df)
        edges_list = split(as.data.table(edges_df), edges_df$in_node)

        founder_cells = rownames(gens1[[1]]$start_barcodes)
        tr_list = purrr::map(founder_cells, function(x_cell) {
                # message(x_cell)
                # skipping to first cell division
                while (nrow(edges_list[[x_cell]]) == 1) {
                        x_cell = edges_list[[x_cell]]$out_node
                        if (!x_cell %in% names(edges_list)) {
                                return(x_cell)
                        }
                }
                sp0 = list(edges_list[[x_cell]][1, ],
                           edges_list[[x_cell]][2, ])
                n_tip = 0
                node_count = 0
                edges_simple_collect = list()
                while(TRUE) {
                        # p0 = bplapply(sp0, get_next, edges_list = edges_list, BPPARAM = bp_param)
                        p0 = lapply(sp0, get_next, edges_list = edges_list)
                        # for (x in sp0) {
                        #         print(x)
                        #         get_next(x, edges_list = edges_list)
                        # }
                        # table(sapply(p0, "[[", "tip"))
                        # as.character(sapply(p0, "[[", "node"))
                        # table(sapply(p0, "[[", "length"))

                        # p0 = clusterApply(cl, sp0, get_next, edges_list = edges_list)
                        edges_simple_collect = append(edges_simple_collect,
                                                      list(do.call(rbind, mapply(function(x, y) c(in_node = as.character(x$in_node),
                                                                                                  out_node = y$node,
                                                                                                  length = y$length),
                                                                                 sp0, p0, SIMPLIFY = F))))
                        p0_node = p0[!sapply(p0, "[[", "tip")]
                        tip_count = sum(sapply(p0, "[[", "tip"))
                        node_count = length(p0_node)
                        n_tip = n_tip + tip_count
                        # message(paste0("merging_node: ", node_count))
                        if (length(p0_node) > 0) {
                                sp1 = do.call(c, lapply(edges_list[as.character(sapply(p0_node, "[[", "node"))],
                                                        function(x) list(x[1, ], x[2, ])))
                                sp0 = sp1
                        } else {
                                # assertthat::assert_that(sum(nodes_df$leaf) == n_tip)
                                break()
                        }
                }
                edges_df_simple = data.frame(do.call(rbind, edges_simple_collect))
                edges_df_simple$length = as.numeric(as.character(edges_df_simple$length))
                edges_df_simple$in_node = as.character(edges_df_simple$in_node)
                edges_df_simple$out_node = as.character(edges_df_simple$out_node)

                tr_nw = edges_df_to_newick(edges_df_simple, nodes_df[nodes_df$id %in% c(edges_df_simple$in_node, edges_df_simple$out_node), ])
                tr_nw = paste0(tr_nw, nodes_df[nodes_df$id == x_cell, "time"], ";")
                tr = ape::read.tree(text=tr_nw)
                tr
        })
        tr_list
}

# appended additional functions below
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




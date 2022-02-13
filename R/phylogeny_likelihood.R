# get the seuqnce of nodes from a tip to root
get_node_seq <- function(tip_id, edge_df) {
        node_cur = tip_id
        node_seq = character()
        while(node_cur != "N-1") {
                node_cur = edge_df$from[edge_df$to == node_cur]
                node_seq = c(node_seq, node_cur)
        }
        node_seq
}
# generate a sequence of merges based on edge mat
make_merge_sequences <- function(edge_mat, tip_id) {
        merge_sequence = list()
        eligible_nodes = tip_id
        while (length(eligible_nodes) > 1) {
                p0_df = edge_mat[edge_mat$out_node %in% eligible_nodes,]
                p0_count = table(p0_df$in_node)
                # nodes to be merged this round
                p0_merging = names(p0_count)[p0_count == 2]
                p0_df_merge = p0_df[p0_df$in_node %in% p0_merging, ]
                cur_merge_list = split(p0_df_merge$out_node, p0_df_merge$in_node)
                merge_sequence = append(merge_sequence, list(cur_merge_list))
                eligible_nodes1 = c(eligible_nodes[!eligible_nodes %in% p0_df_merge$out_node],
                                    unique(p0_df_merge$in_node))
                eligible_nodes = eligible_nodes1
        }
        merge_sequence
}
# compute allele prob at each node given mutation parameters
#' @param recur_vec_list a list of named allele probability vectors
compute_allele_prob <- function(tree_data, recur_vec_list) {
        sc_mat = tree_data$sc_mat
        node_tip_list = tree_data$tr_proc$node_tip_list
        assertthat::assert_that(length(recur_vec_list) == ncol(sc_mat))
        node_allele_prob_mat = reduce(map(1:ncol(sc_mat), function(k) {
                node_allele_prob = purrr::map_dbl(node_tip_list, function(x) {
                        y = unique(sc_mat[x, k])
                        if (length(y) == 1 & !all(y == "0")) {
                                return(recur_vec_list[[k]][y])
                        } else {
                                return(0)
                        }
                })
                node_allele_prob
        }), rbind)
        colnames(node_allele_prob_mat) = names(node_tip_list)
        tip_allele_prob_mat = reduce(map(1:ncol(sc_mat), function(k) {
                out = recur_vec_list[[k]][sc_mat[, k]]
                out[is.na(out)] = 0
                out
        }), rbind)
        colnames(tip_allele_prob_mat) = rownames(sc_mat)
        allele_prob_mat = cbind(node_allele_prob_mat, tip_allele_prob_mat)
        tree_data$allele_prob_mat = allele_prob_mat
        tree_data
}
# preprocess tree
preproccess_tr <- function(tree_data) {
        tr_r = tree_data$tr
        # rename tr nodes
        tr_r$node.label = paste0("N-", 1:(tr_r$Nnode))
        # list tips for each node
        node_tip_list = list_dd_and_tips(tr_r)$tips
        # node_tip_list = map(adephylo::listTips(tr_r), names)
        # get edge_df with node names
        edge_mapper = character(max(tr_r$edge))
        edge_mapper[(1:Nnode(tr_r))+Ntip(tr_r)] = tr_r$node.label
        edge_mapper[1:Ntip(tr_r)] = tr_r$tip.label
        edge_df = tibble(from = edge_mapper[tr_r$edge[, 1]],
                         to = edge_mapper[tr_r$edge[, 2]])
        edge_df$length = tr_r$edge.length
        edge_df = bind_rows(tibble(from = "root", to = "N-1", length = tree_data$root_edge_len),
                            edge_df)
        m_seq = make_merge_sequences(dplyr::rename(edge_df, out_node = to, in_node = from),
                                     tip_id = tree_data$tip_id)
        tip_node_seq = map(tree_data$tip_id, get_node_seq, edge_df = edge_df)
        edge_df_node = edge_df[!edge_df$to %in% tree_data$tip_id, ]

        tree_data$tr_proc = list(edge_df = edge_df,
                                 edge_df_node = edge_df_node,
                                 m_seq = m_seq,
                                 tip_node_seq = tip_node_seq,
                                 node_tip_list = node_tip_list)
        tree_data
}
#' @param node_edge_len edge lengths vector for nodes, in the order of edge_df_node
compute_loglike <- function(tree_data, mut_rate, node_edge_len = NULL) {
        edge_df = tree_data$tr_proc$edge_df
        edge_df_node = tree_data$tr_proc$edge_df_node
        tip_node_seq = tree_data$tr_proc$tip_node_seq
        sc_mat = tree_data$sc_mat
        assertthat::assert_that(length(mut_rate) == ncol(sc_mat))
        if (is.null(node_edge_len)) {
                node_edge_len = edge_df_node$length
        }
        assertthat::assert_that(length(node_edge_len) == nrow(edge_df_node))
        tip_edge_len = map_dbl(tip_node_seq, function(x) {
                tree_data$total_time - sum(node_edge_len[edge_df_node$to %in% x])
        })
        # reorder edge_len so that it is in the order of edge_df
        edge_len = c(node_edge_len, tip_edge_len)[match(edge_df$to, c(edge_df_node$to, tree_data$tip_id))]
        assertthat::assert_that(all(edge_len > 0))

        # computation begins here
        p_pos = log((1 - exp(- cbind(mut_rate) %*% edge_len)) * tree_data$allele_prob_mat[, edge_df$to])
        p_neg = -cbind(mut_rate) %*% edge_len
        colnames(p_neg) = colnames(p_pos) = edge_df$to
        # asserting no chance of mutation at root node
        p_neg[, "N-1"] = 0
        #1 matrix to store merged probabilities
        p_node = matrix(NA, length(mut_rate), nrow(edge_df))
        colnames(p_node) = edge_df$to
        #2 asssign tip probabilities
        p_node[, tree_data$tip_id] = reduce(map(1:length(tree_data$tip_id), function(i) {
                neg_vec = p_neg[, tree_data$tip_id[i]]
                neg_vec[sc_mat[i, ] != "0"] = -Inf
                matrixStats::rowLogSumExps(cbind(p_pos[, tree_data$tip_id[i]],
                                                 neg_vec))
        }), cbind)
        # merging node probabilites
        for(m_list in tree_data$tr_proc$m_seq) {
                for (m in names(m_list)) {
                        out_node1 = m_list[[m]][1]
                        out_node2 = m_list[[m]][2]
                        in_node = m
                        p_node[, in_node] = matrixStats::rowLogSumExps(cbind(p_node[, out_node1] + p_node[, out_node2] + p_neg[, in_node],
                                                                             p_pos[, in_node]))
                }
        }
        return(sum(p_node[, "N-1"]))
}
# constraint in the form of u %*% edge_len_node - c >= 0
build_constraint_matrices <- function(tree_data) {
        edge_df_node = tree_data$tr_proc$edge_df_node
        tip_node_seq = tree_data$tr_proc$tip_node_seq
        u_mat = matrix(0, length(tree_data$tip_id), nrow(edge_df_node))
        for (i in 1:length(tip_node_seq)) {
                u_mat[i, match(tip_node_seq[[i]], edge_df_node$to)] = -1
        }
        c_mat = - cbind(rep(tree_data$total_time, length(tree_data$tip_id)))
        tree_data$edge_len_constr = list(u = u_mat,
                                         c = c_mat)
        tree_data
}
# produce a new tree object using provided node_edge_len parameters
update_tree_edge_len <- function(tree_data, node_edge_len) {
        edge_df = tree_data$tr_proc$edge_df
        edge_df_node = tree_data$tr_proc$edge_df_node
        tip_node_seq = tree_data$tr_proc$tip_node_seq
        tip_edge_len = map_dbl(tip_node_seq, function(x) {
                tree_data$total_time - sum(node_edge_len[edge_df_node$to %in% x])
        })
        edge_len = c(node_edge_len, tip_edge_len)[match(edge_df$to, c(edge_df_node$to, tree_data$tip_id))]
        tr = rlang::duplicate(tree_data$tr)
        tr$edge.length = edge_len[2:length(edge_len)]
        tr$root.edge = edge_len[1]
        tr
}

# wrap up everything to evaluate a tree
compute_tree_loglike <- function(tr, sc_mat, total_time, root_edge_len) {
        tree_data = make_tree_data(tr, sc_mat, rownames(sc_mat), total_time, root_edge_len)
        tree_data = compute_allele_prob(tree_data, mut_p$recur_vec_list)
        compute_loglike(tree_data, mut_p$mut_rate)
}
# tree_data constructor
make_tree_data <- function(tr, sc_mat, tip_id, total_time, root_edge_len, preprocess = T) {
        out = list(tr = tr,
                   sc_mat = sc_mat,
                   tip_id = rownames(sc_mat),
                   total_time = total_time,
                   root_edge_len = root_edge_len)
        if(preprocess) {
                out = preproccess_tr(out)
        }
        out
}





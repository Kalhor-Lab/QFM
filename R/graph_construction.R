sample_node_ranking <- function(gr) {
        list_out = list_dd_and_tips(name_nodes(gr))
        gr_dd = list_out$dd
        gr_tip_lists = list_out$tips
        # genrate the draw for each node first
        node_shuffle = map(gr_dd, function(x) {
                s1 = pmax(1, length(gr_tip_lists[[x[1]]])) - 1
                s2 = pmax(1, length(gr_tip_lists[[x[2]]])) - 1
                sample(c(rep(x[1], s1), rep(x[2], s2)))
        })
        root_node = names(gr_dd)[1]
        ranking = rep(root_node, gr$Nnode)
        assign_dd <- function(dd) {
                if (dd %in% gr$tip.label) {
                        return()
                } else {
                        dd_indices = which(ranking == dd)
                        if (length(dd_indices) <= 1) {
                                return()
                        } else {
                                ranking[dd_indices[-1]] <<- node_shuffle[[dd]]
                                assign_dd(gr_dd[[dd]][1])
                                assign_dd(gr_dd[[dd]][2])
                        }
                }
        }
        assign_dd(root_node)
        ranking
}
# Steps for fate map generation:
# 1. generate topology
# 2. sample node time ranking
# 3. sample node time in interval
# this function attempts to fit in time interval based on the random node ordering sampled
# if failed it returns the minimum intervals needed
make_graph_from_phy_mod2 <- function(phy,
                                     node_time_start = 2.5,
                                     node_time_end = 10.9,
                                     target_time = 11.5,
                                     node_rank = NULL) {
        num_tip = length(phy$tip.label)
        phy$node.label = paste0("N-", 1:(phy$Nnode))
        if (is.null(node_rank)) {
                phy_node_rank = sample_node_ranking(phy)
        } else {
                phy_node_rank = node_rank
        }
        node_doubletime = map(seq(0.55, 0.4, length = length(phy_node_rank)), function(x) {
                runif(1, min = x - 0.05, max = x + 0.05)
        })
        names(node_doubletime) = phy_node_rank
        node_doubletime[["N-1"]] = 0.6
        tip_doubletime = map(phy$tip.label, function(x) runif(1, 0.35, 0.45))
        names(tip_doubletime) = phy$tip.label
        type_doubletime = c(node_doubletime, tip_doubletime)

        # type_doubletime = map(c(phy$node.label, phy$tip.label), function(x) {
        #         if (x == "N-1") {
        #                 return(0.6)
        #                 # return(double_time + root_padding)
        #         } else {
        #                 if (x %in% phy$tip.label) {
        #                         return(runif(1, 0.35, 0.4))
        #                 }
        #                 if (match(x, phy_node_rank) <= floor(length(phy$node.label) * 0.35)) {
        #                         # earliest 20% of nodes
        #                         return(runif(1, 0.45, 0.5))
        #                 } else {
        #                         return()
        #                 }
        #                 # return(double_time)
        #         }
        #         # if(x < 0) {
        #         #         return(0.35)
        #         # } else {
        #         #         if (type_time_list[[x]] < node_time_start) {
        #         #                 return(0.6)
        #         #         } else {
        #         #                 return(0.35)
        #         #         }
        #         # }
        # })
        # names(type_doubletime) = c(phy$node.label, phy$tip.label)
        # not all intervals need to be double_time apart
        # now get the intervals that are actually edges
        phy_edges = make_edge_df(phy)
        phy_node_rank_map = 1:length(phy_node_rank)
        names(phy_node_rank_map) = phy_node_rank
        # try listing all constraints
        # linear program to solve for the minimum constrained intervals
        # parameters: x_1, ... x_(Nnode-1)
        Amat = matrix(0, nrow(phy_edges), (length(phy_node_rank) - 1))
        bvec = numeric(nrow(phy_edges))
        for (i in 1:nrow(phy_edges)) {
                if (!grepl("t", phy_edges$to[i])) {
                        start_index = phy_node_rank_map[phy_edges$from[i]]
                        end_index = phy_node_rank_map[phy_edges$to[i]]
                        intervals_index = start_index:(end_index-1)
                        Amat[i, intervals_index] = 1
                        bvec[i] = type_doubletime[[phy_edges$from[i]]]
                }
        }
        ind = rowSums(Amat) > 0
        Amat = Amat[ind, ]
        bvec = bvec[ind]
        # bvec = rep(double_time, nrow(Amat))
        cvec = rep(1, ncol(Amat))
        lp_out = lp(objective.in = cvec, const.mat = Amat, const.rhs = bvec, const.dir = ">=")
        min_interval = lp_out$solution

        if ((node_time_end - node_time_start) <= sum(min_interval)) {
                return(sum(min_interval))
        }
        node_intervals = extraDistr::rdirichlet(n = 1, rep(5, length(phy_node_rank) + 1))[1, ] *
                (node_time_end - node_time_start - sum(min_interval))
        node_intervals[2:Nnode(phy)] = node_intervals[2:Nnode(phy)] + min_interval
        node_time = node_time_start + cumsum(node_intervals)
        node_time = node_time[1:(length(node_time)-1)]
        names(node_time) = phy_node_rank
        tip_time = rep(Inf, num_tip)
        names(tip_time) = phy$tip.label
        type_time = c(node_time, tip_time)

        tr_r = phy
        root_edge_len = node_time["N-1"]
        # list tips for each node
        # get edge_df with node names
        edge_mapper = character(max(tr_r$edge))
        edge_mapper[(1:Nnode(tr_r))+Ntip(tr_r)] = tr_r$node.label
        edge_mapper[1:Ntip(tr_r)] = tr_r$tip.label

        edge_df = tibble(from = edge_mapper[tr_r$edge[, 1]],
                         to = edge_mapper[tr_r$edge[, 2]])
        edge_df$length = tr_r$edge.length
        edge_df = bind_rows(tibble(from = "root", to = "N-1", length = root_edge_len),
                            edge_df)
        m_seq = make_merge_sequences(dplyr::rename(edge_df, out_node = to, in_node = from),
                                     tip_id = tr_r$tip.label)
        merge_list = do.call(c, m_seq)
        node_mapper = as.character(c(1:length(merge_list), c(-1:-(Ntip(phy)))))
        names(node_mapper) = c(names(merge_list), phy$tip.label)

        merge_list = map(merge_list, function(x) node_mapper[x])
        names(merge_list) = node_mapper[names(merge_list)]

        type_time_list = lapply(type_time, function(x) x)
        names(type_time_list) = node_mapper[names(type_time_list)]

        type_diff_mode_probs = map(type_time[phy$node.label], function(x) {
                a = 5
                p = rbeta(n = 1, a, a)
                c(p, 1-p, 0)
        })
        names(type_diff_mode_probs) = node_mapper[names(type_diff_mode_probs)]

        root_id = as.character(num_tip-1)
        node_id = as.character(1:(num_tip-1))
        tip_id = as.character(-1:-(num_tip))

        # rename double time
        names(type_doubletime) = node_mapper[names(type_doubletime)]

        big_graph = list(
                root_id = root_id,
                node_id = node_id,
                tip_id = tip_id,
                all_id = c(node_id, tip_id),
                merge = merge_list,
                diff_time = type_time_list,
                diff_mode_probs = type_diff_mode_probs,
                lifetime = type_doubletime,
                founder_size = 1,
                target_time = target_time
        )
        big_graph$prob_loss = map(big_graph$all_id, function(x) {
                if (x == big_graph$root_id) {
                        return(0)
                } else {
                        return(runif(1, 0.02, 0.08))
                }
        })
        names(big_graph$prob_loss) = big_graph$all_id
        big_graph$prob_pm = map(big_graph$all_id, function(x) return(0.00))
        names(big_graph$prob_pm) = big_graph$all_id

        class(big_graph) = "type_graph"
        # min_diff = 0.01
        # for (x in forward_merge_sequence(big_graph)) {
        #         print(x)
        #         if((big_graph$diff_time[[as.character(big_graph$merge[x, 1])]] - big_graph$diff_time[[x]]) < big_graph$lifetime[[x]] + min_diff) {
        #                 big_graph$diff_time[[as.character(big_graph$merge[x, 1])]] = big_graph$diff_time[[x]] + big_graph$lifetime[[x]] + min_diff
        #         }
        #         if((big_graph$diff_time[[as.character(big_graph$merge[x, 2])]] - big_graph$diff_time[[x]]) < big_graph$lifetime[[x]] + min_diff) {
        #                 big_graph$diff_time[[as.character(big_graph$merge[x, 2])]] = big_graph$diff_time[[x]] + big_graph$lifetime[[x]] + min_diff
        #         }
        # }
        big_graph = generate_edge_df_mod2(big_graph)
        big_graph
}
adjust_pec_v1 <- function(phy, move = TBR) {
        temp = make_graph_from_phy_mod2(phy)
        if (class(temp) == "type_graph") {
                return(temp)
        } else {
                message("attemping adjust topology.")
                assertthat::assert_that(class(temp) == "numeric")
                cur_min_total_time = temp
        }
        message(cur_min_total_time)
        fail_count = 0
        while(TRUE) {
                phy_new = move(phy) # perturb phy
                new_temp = make_graph_from_phy_mod2(phy_new)
                if (class(new_temp) == "type_graph") {
                        return(new_temp)
                } else {
                        # accept perturb if total min_time has decreased
                        assertthat::assert_that(class(new_temp) == "numeric")
                        new_min_total_time = new_temp
                        if (new_min_total_time < cur_min_total_time) {
                                message(new_min_total_time)
                                cur_min_total_time = new_min_total_time
                                phy = phy_new
                                fail_count = 0
                        } else {
                                fail_count = fail_count + 1
                                if (fail_count >= 500) {
                                        message("failed")
                                        return(NULL)
                                }
                        }
                }
        }
}

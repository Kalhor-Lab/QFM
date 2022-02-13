source("analysis/sample_MARC1_hgRNA.R")

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
make_graph_from_phy <- function(phy,
                                node_time_start = 1.8,
                                node_time_end = 9.8,
                                double_time = 0.35,
                                target_time = 15.,
                                root_padding = 0.1) {
        num_tip = length(phy$tip.label)
        phy$node.label = paste0("N-", 1:(phy$Nnode))
        phy_node_rank = sample_node_ranking(phy)
        # not all intervals need to be double_time apart
        # now get the intervals that are actually edges
        phy_edges = make_edge_df(phy)

        phy_node_rank_map = 1:length(phy_node_rank)
        names(phy_node_rank_map) = phy_node_rank
        # try listing all constraints
        # linear program to solve for the minimum constrained intervals
        # parameters: x_1, ... x_(Nnode-1)
        Amat = matrix(0, nrow(phy_edges), (length(phy_node_rank) - 1))
        for (i in 1:nrow(phy_edges)) {
                if (!grepl("t", phy_edges$to[i])) {
                        start_index = phy_node_rank_map[phy_edges$from[i]]
                        end_index = phy_node_rank_map[phy_edges$to[i]]
                        intervals_index = start_index:(end_index-1)
                        Amat[i, intervals_index] = 1
                }
        }
        Amat = Amat[rowSums(Amat) > 0, ]
        bvec = rep(double_time, nrow(Amat))
        cvec = rep(1, ncol(Amat))
        lp_out = lp(objective.in = cvec, const.mat = Amat, const.rhs = bvec, const.dir = ">=")
        min_interval = lp_out$solution
        min_interval[1] = min_interval[1] + root_padding

        assertthat::assert_that((node_time_end - node_time_start) > sum(min_interval))
        node_intervals = extraDistr::rdirichlet(n = 1, rep(5, length(phy_node_rank) + 1))[1, ] *
                (node_time_end - node_time_start - sum(min_interval))
        node_intervals[2:Nnode(phy)] = node_intervals[2:Nnode(phy)] + min_interval
        node_time = node_time_start + cumsum(node_intervals)
        node_time = node_time[1:(length(node_time)-1)]
        names(node_time) = phy_node_rank
        tip_time = rep(Inf, num_tip)
        names(tip_time) = phy$tip.label
        type_time = c(node_time, tip_time)

        # make merge_mat
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
        m_seq = make_merge_sequences(rename(edge_df, out_node = to, in_node = from),
                                     tip_id = tr_r$tip.label)
        merge_mat = do.call(rbind, map(m_seq, function(x) do.call(rbind, x)))
        node_mapper = as.character(c(1:nrow(merge_mat), c(-1:-(Ntip(phy)))))
        names(node_mapper) = c(rownames(merge_mat), phy$tip.label)
        merge_mat = matrix(as.numeric(node_mapper[c(merge_mat)]), ncol = 2)
        rownames(merge_mat) = 1:nrow(merge_mat)

        type_time_list = lapply(type_time, function(x) x)
        names(type_time_list) = node_mapper[names(type_time_list)]

        type_diff_mode_probs = map(type_time[phy$node.label], function(x) {
                a = 5
                p = rbeta(n = 1, a, a)
                c(p, 1-p, 0, 0, 0)
        })
        names(type_diff_mode_probs) = node_mapper[names(type_diff_mode_probs)]
        type_doubletime = lapply(as.character(c(1:(Ntip(phy)-1), -1:-Ntip(phy))), function(x) {
                if (x == as.character(Ntip(phy) - 1)) {
                        return(double_time + root_padding)
                } else {
                        return(double_time)
                }
                # if(x < 0) {
                #         return(0.35)
                # } else {
                #         if (type_time_list[[x]] < node_time_start) {
                #                 return(0.6)
                #         } else {
                #                 return(0.35)
                #         }
                # }
        })
        names(type_doubletime) = as.character(c(1:(num_tip-1), -1:-(num_tip)))

        root_id = as.character(num_tip-1)
        node_id = as.character(1:(num_tip-1))
        tip_id = as.character(-1:-(num_tip))
        big_graph = list(
                root_id = root_id,
                node_id = node_id,
                tip_id = tip_id,
                all_id = c(node_id, tip_id),
                merge = merge_mat,
                diff_time = type_time_list,
                diff_mode_probs = type_diff_mode_probs,
                lifetime = type_doubletime,
                founder_size = 1,
                target_time = 15.0
        )
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
        big_graph = generate_diff_outcome(big_graph)
        big_graph = generate_edge_df(big_graph)
        big_graph
}
calc_bsum <- function(phy) {
        sum(apply(balance(phy), 1, function(x) abs(x[1] - x[2])))
}
# *****
output_dir = "./intermediate_data/panel/"
# dir.create(output_dir)
# *****

library(ape)
library(lpSolve)
num_total = 10000
set.seed(73)
graph_list_all = map(c(16, 32, 64), function(num_tip) {
        random_trees = map(1:num_total, function(i) {
                rtree(n = num_tip, equiprob = F)
        })
        bsum_vec = map_dbl(random_trees, calc_bsum)
        bsum_bal = calc_bsum(stree(num_tip, "balanced"))
        bsum_pec = calc_bsum(stree(num_tip, "left"))
        bsum_cat = cut(bsum_vec, breaks = seq(from = bsum_bal,
                                              to = bsum_pec,
                                              by = 5))
        random_indices = reduce(map(levels(bsum_cat), function(x) {
                indices = which(bsum_cat == x)
                sample(indices, size = pmin(length(indices), 1), replace = F)
        }), c)
        # length(random_indices)
        # hist(c(bsum_bal, bsum_vec[random_indices], bsum_pec), breaks = 200)
        phy_list = c(list(stree(num_tip, "balanced")),
                     random_trees[random_indices])
        # if (num_tip <= 10) {
        #         phy_list = c(phy_list, list(stree(num_tip, "left")))
        # }
        graph_list = map(1:length(phy_list), function(i) {
                print(i)
                make_graph_from_phy(phy_list[[i]])
        })
        graph_list
})
graph_list_all = reduce(graph_list_all, c)
saveRDS(graph_list_all, file = paste0(output_dir, "all_graphs.rds"))

# simulating the experiments
library(furrr)
plan(multisession, workers = 12)
graph_list_all = readRDS(paste0(output_dir, "all_graphs.rds"))
exp_params = expand.grid(big_graph_id = 1:length(graph_list_all),
                         sample_size = c(100),
                         num_element = c(100),
                         sampling = c("fixed", "proportional"),
                         i_sim = 1:2) %>% as_tibble()
exp_params = mutate(exp_params, mut_p = map(num_element, function(i) {
        sample_mutp_mid_or_fast(i, c(0.6, 15.0))
}))
saveRDS(exp_params, file = paste0(output_dir, "exp_data.rds"))

# library(furrr)
# plan(multisession, workers = 12)
# rep_dir = paste0(output_dir, "exp_data_rep2/")
# dir.create(rep_dir)
# future_walk(1:nrow(exp_params), function(i) {
#         big_graph_id = exp_params$big_graph_id[i]
#         sample_size = exp_params$sample_size[i]
#         mut_p = exp_params$mut_p[[i]]
#         sampling = exp_params$sampling[i]
#
#         if (sampling == "fixed") {
#                 ss = sample_size
#         }
#         if (sampling == "proportional") {
#                 fm = graph_list_all[[big_graph_id]]
#                 gens0 = make_gens(fm)
#                 tip_size = map_dbl(gens0[fm$tip_id], "end_count")
#                 tip_sample = extraDistr::rmvhyper(nn = 1, k = length(fm$tip_id) * sample_size, n = tip_size)[1, ]
#                 tip_sample = pmax(tip_sample, 5)
#                 names(tip_sample) = fm$tip_id
#                 ss = tip_sample
#         }
#         out_file = paste0(rep_dir, str_pad(i, width = 4, pad = "0"), ".rds")
#         if (!file.exists(out_file)) {
#                 try({
#                         out = simulate_sc_data(graph_list_all[[big_graph_id]], mut_p, sample_size = ss)
#                         saveRDS(out, file = out_file)
#                 })
#         }
# }, .progress = T, .options = furrr_options(seed = T))
#
# out_files = map_chr(1:nrow(exp_params), function(i) {
#         out_file = paste0(rep_dir, str_pad(i, width = 4, pad = "0"), ".rds")
#         out_file
# })
# exp_params[!file.exists(out_files), ]
# exp_params = exp_params %>% mutate(data = future_pmap(., function(big_graph_id, sample_size, mut_p, sampling, ...) {
#         if (sampling == "fixed") {
#                 ss = sample_size
#         }
#         if (sampling == "proportional") {
#                 fm = graph_list_all[[big_graph_id]]
#                 gens0 = make_gens(fm)
#                 tip_size = map_dbl(gens0[fm$tip_id], "end_count")
#                 tip_sample = extraDistr::rmvhyper(nn = 1, k = length(fm$tip_id) * sample_size, n = tip_size)[1, ]
#                 tip_sample = pmax(tip_sample, 5)
#                 names(tip_sample) = fm$tip_id
#                 ss = tip_sample
#         }
#         out = simulate_sc_data(graph_list_all[[big_graph_id]], mut_p, sample_size = ss)
#         saveRDS(out, file = paste0(output_dir, "exp_data_rep1_", , ".rds"))
# }, .progress = T, .options = furrr_options(seed = T)))


exp_params = readRDS(paste0(output_dir, "exp_data.rds"))
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))

exp_params = exp_params %>% mutate(sc = map(data, "sc"))
exp_params = exp_params %>% mutate(tr = map(data, "tr"))

total_time = 15.0 - 0.6
for (j in 1:nrow(exp_params)) {
        message(j)
        print(paste0("Ntip: ", length(all_graphs[[exp_params$big_graph_id[j]]]$tip_id)))
        tr3 = name_nodes(phylotime(exp_params$sc[[j]], total_time, mut_p = exp_params1$mut_p[[j]]))
        saveRDS(tr3,
                file = paste0(output_dir, "temp_tr3/", stringr::str_pad(j, width = 3, side = "left", "0")))
}
exp_params = mutate(exp_params, tr3 = map2(sc, mut_p, function(x, m) {
        name_nodes(phylotime(x, total_time, mut_p = m))
}))
saveRDS(exp_params, file = paste0(output_dir, "exp_data_proc.rds"))




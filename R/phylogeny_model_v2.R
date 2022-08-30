#' current version of cell phylogeny model

# A function that plots the number of cells for each state at any given time point
get_cell_type_counts <- function(type_graph, tp) {
        gens0 = make_gens_mod2(type_graph)
        purrr::reduce(map(gens0, function(x) {
                if (tp < x$start_time | tp >= (x$end_time + x$double_time)) {
                        out = 0
                } else {
                        out = sum(x$cell_mode_counts[(floor((tp - x$start_time) / x$double_time) + 1), c(1, 3)])
                }
                names(out) = x$cell_type
                out
        }), c)
}
# end function
# a function that generates the offspring matrix
offspring_mat = function(n) {
        rbind(diag(2, nrow = n),
              do.call(rbind, map(1:(n-1), function(i) {
                      do.call(rbind, map((i+1):n, function(j) {
                              out = numeric(n)
                              out[i] = 1
                              out[j] = 1
                              out
                      }))
              })))
}
distribute_mode_counts_mod2 <- function(cell_type, mode_counts, type_graph) {
        m = length(type_graph$merge[[cell_type]])
        assertthat::assert_that(length(mode_counts) == (m + choose(m, 2)))
        daughter_mode_counts = mode_counts * offspring_mat(m)
        colnames(daughter_mode_counts) =
                type_graph$merge[[cell_type]]
        daughter_mode_counts
}
# next_cell_type_gens_mod2 <- function (x, type_graph) {
#         assertthat::assert_that(x$active)
#         # TODO: naming should be changed for the future
#         # the mode here is in the sense of commitment not, singlet vs doublet
#         prob_loss = type_graph$prob_loss[[x$cell_type]]
#         prob_pm = type_graph$prob_pm[[x$cell_type]]
#         daughter_mode_counts = distribute_mode_counts_mod2(x$cell_type,
#                                                            x$end_mode_counts,
#                                                            type_graph)
#         out_gens = map(1:ncol(daughter_mode_counts), function(j) {
#                 cell_count_last = sum(daughter_mode_counts[, j]) # this is the total for the first generation
#                 assertthat::assert_that(cell_count_last > 0)
#                 n_loss = floor(cell_count_last * prob_loss)
#                 n_pm = floor(cell_count_last * prob_pm)
#                 n_double = cell_count_last - n_loss - n_pm
#                 assertthat::assert_that((n_double + n_pm) > 0)
#                 make_cell_type_gens_mod2(colnames(daughter_mode_counts)[j],
#                                          start_time = x$end_time + x$double_time,
#                                          start_count = NULL,
#                                          start_mode_counts = c(n_double, n_loss, n_pm),
#                                          commit_mode_counts = daughter_mode_counts[, j],
#                                          type_graph = type_graph)
#         })
#         return(out_gens)
# }
next_cell_type_gens_mod2 <- function (x, type_graph) {
        assertthat::assert_that(x$active)
        # TODO: naming should be changed for the future
        # the mode here is in the sense of commitment not, singlet vs doublet
        daughter_mode_counts = distribute_mode_counts_mod2(x$cell_type,
                                                           x$end_mode_counts,
                                                           type_graph)
        out_gens = map(1:ncol(daughter_mode_counts), function(j) {
                cell_type_j = colnames(daughter_mode_counts)[j]
                cell_count_last = sum(daughter_mode_counts[, j]) # this is the total for the first generation
                assertthat::assert_that(cell_count_last > 0)
                prob_loss = type_graph$prob_loss[[cell_type_j]]
                prob_pm = type_graph$prob_pm[[cell_type_j]]
                n_loss = floor(cell_count_last * prob_loss)
                n_pm = floor(cell_count_last * prob_pm)
                n_double = cell_count_last - n_loss - n_pm
                assertthat::assert_that((n_double + n_pm) > 0)
                make_cell_type_gens_mod2(cell_type_j,
                                         start_time = x$end_time + x$double_time,
                                         start_count = NULL,
                                         start_mode_counts = c(n_double, n_loss, n_pm),
                                         commit_mode_counts = daughter_mode_counts[, j],
                                         type_graph = type_graph)
        })
        return(out_gens)
}
make_cell_type_gens_mod2 <- function(cell_type,
                                     start_time,
                                     start_count,
                                     start_mode_counts,
                                     commit_mode_counts,
                                     type_graph) {
        assertthat::assert_that(start_time < type_graph$target_time)
        prob_loss = type_graph$prob_loss[[cell_type]]
        prob_pm = type_graph$prob_pm[[cell_type]]
        # message(cell_type)
        # FIXED rates, now life time is simply a discrete time interval
        gen_state = num_gen_single_type(cell_type, start_time, type_graph)
        # if (cell_type == "12") {
        #         print(gen_state)
        # }
        life_duration = type_graph$lifetime[[cell_type]]
        # not counting the division as the cells commit, as spliting is needed
        eff_gen = ifelse(gen_state$diff,
                         gen_state$num_gen - 1,
                         gen_state$num_gen)

        # setting up mode count matrix
        # three columns for three modes: Dividing, Loss, Post-mitotic, chance of three modes add up to one
        # fixed chacne per cell cycle
        # ncol is equal to eff_gen+1
        cell_mode_counts = matrix(NA, nrow = (eff_gen+1), ncol = 3)
        cell_mode_counts[1, ] = start_mode_counts
        if (eff_gen > 0) {
                for (i in 1:eff_gen) {
                        cell_count_last = cell_mode_counts[i, 1] * 2 + cell_mode_counts[i, 3]
                        # drawing number of L and P
                        # TODO: add cases of unsampled lineage
                        # TODO: make the process stochastic
                        # n_loss = rbinom(n = 1, size = cell_count_last, prob = prob_loss)
                        n_loss = floor(cell_count_last * prob_loss)
                        n_pm = floor(cell_count_last * prob_pm)
                        n_double = cell_count_last - n_loss - n_pm
                        assertthat::assert_that((n_double + n_pm) > 0)
                        cell_mode_counts[i+1, ] = c(n_double, n_loss, n_pm)
                }
                #0 these counts do not completely decide the shape of the population tree, extra step is needed
        }
        end_count = cell_mode_counts[eff_gen+1, 1] + cell_mode_counts[eff_gen+1, 3] # not including cells that are dead at the last generation
        # here the modes are mode of commitment, those commiting to single types always doubles
        if (gen_state$diff) {
                end_mode_counts = get_mode_counts(cell_type, end_count, type_graph)
        } else {
                end_mode_counts = NULL
        }
        out = list(cell_type = cell_type,
                   start_time = start_time,
                   start_count = start_count,
                   start_mode_counts = start_mode_counts,
                   commit_mode_counts = commit_mode_counts,
                   num_gen = eff_gen, # NOTE: this num_gen is actually eff_gen, need to rename the previous variable to avoid confusion
                   double_time = life_duration,
                   cell_counts = NULL,
                   cell_mode_counts = cell_mode_counts,
                   end_time = start_time + life_duration * eff_gen,
                   end_count = end_count,
                   end_mode_counts = end_mode_counts,
                   active = gen_state$diff
        )
        class(out) = "cell_type_gen"
        out
}
make_gens_mod2 = function(type_graph) {
        init_gens = make_cell_type_gens_mod2(type_graph$root_id, start_time = 0,
                                             start_count = type_graph$founder_size,
                                             start_mode_counts = c(type_graph$founder_size, 0, 0),
                                             commit_mode_counts = NA,
                                             type_graph)
        collect_gens = list()
        collect_gens = append(collect_gens, list(init_gens))
        active_gens = list(init_gens)
        while (length(active_gens) > 0) {
                next_gens = unlist(lapply(active_gens, function(x)
                        next_cell_type_gens_mod2(x,
                                                 type_graph = type_graph)),
                        recursive = F)
                collect_gens = append(collect_gens, next_gens)
                active_ind = sapply(next_gens, "[[", "active")
                active_gens = next_gens[active_ind]
        }
        names(collect_gens) = sapply(collect_gens, "[[", "cell_type")
        collect_gens
}
backward_merge_sequence_mod2 <- function (type_graph)
{
        node_seq = character()
        cur_nodes = type_graph$tip_id
        while (length(node_seq) < length(type_graph$node_id)) {
                for (node_id in type_graph$node_id) {
                        node_dau = type_graph$merge[[node_id]]
                        if (all(node_dau %in% cur_nodes)) {
                                cur_nodes = cur_nodes[!cur_nodes %in% node_dau]
                                cur_nodes = c(cur_nodes, node_id)
                                node_seq = append(node_seq, node_id)
                        }
                }
        }
        return(node_seq)
}
generate_sample_size_mod2 <- function (x, sample_size) {
        assertthat::assert_that("matrix" %in% class(x$cell_mode_counts))
        cell_mode_sample_size = matrix(NA, nrow = x$num_gen, ncol = 3)

        # redesigned V1
        # meaning for columns in the cell_mode_sample_size:
        # column1: number of doublets
        # column2: number of singlets (failed to merge singlets + sure singlets)

        # death cells contribute to doublet and singlet and is make up the total of
        # x$cell_mode_counts[x$num_gen, 1] * 2 +x$cell_mode_counts[x$num_gen, 3]
        sample_size_vec = sample_size
        sample_size_cur = sample_size
        if (x$num_gen > 0) {
                for (i in x$num_gen:1) {
                        temp = extraDistr::rmvhyper(nn = 1,
                                                    k = sample_size_cur,
                                                    c(x$cell_mode_counts[i, 1] * 2,
                                                      x$cell_mode_counts[i, 3]))[1, ]
                        potential_double = temp[1]
                        sure_single = temp[2]
                        num_merges = sample_doublet_branch(x$cell_mode_counts[i, 1], potential_double)
                        potential_single = potential_double - 2 * num_merges
                        cell_mode_sample_size[i, ] = c(num_merges, potential_single, sure_single)

                        sample_size_cur = sum(cell_mode_sample_size[i, ])
                        sample_size_vec = c(sample_size_vec, sample_size_cur)
                }
                x$sample_size_mode = cell_mode_sample_size
        }
        x$sample_size = sample_size_vec # this is appenped in reverse order
        # start mode sample size is split into commitment modes
        if (all(!is.na(x$commit_mode_counts))) { # only root
                commit_mode_sample_size = extraDistr::rmvhyper(nn = 1, k = x$sample_size[x$num_gen + 1],
                                                               x$commit_mode_counts)[1, ]
                x$commit_mode_sample_size = commit_mode_sample_size
        } else {
                x$commit_mode_sample_size = NA
        }
        x
}
merge_sample_size_mod2 <- function (x_ls, y)
{
        m = length(x_ls) # multiplicity
        commit_mode_ss_mat = do.call(cbind, map(x_ls, "commit_mode_sample_size"))
        # symmetric
        s_vec = diag(commit_mode_ss_mat[1:m, ])
        z_vec = map_dbl(1:m, function(i) {
                sample_doublet_branch(y$end_mode_counts[i], s_vec[i])
        })
        m_vec = s_vec - z_vec
        # assymetric
        sa_ls = map(1:choose(m, 2), function(j) {
                x_ind = offspring_mat(m)[m + j, ] > 0 #pair of downstream states for this mode
                s_a = commit_mode_ss_mat[m + j, x_ind]
                s_a
        })
        za_vec = map_dbl(1:choose(m, 2), function(j) {
                sample_doublet_twotype(n = y$end_mode_counts[m + j], sa_ls[[j]][1], sa_ls[[j]][2])
        })
        ma_vec = map_dbl(sa_ls, sum) - za_vec

        y$end_mode_sample_size = c(m_vec, ma_vec) # just for record keeping
        y = generate_sample_size_mod2(y, sample_size = sum(c(m_vec, ma_vec)))
        y
}
calculate_type_counts_mod2 <- function(type_graph) {
        gens0 = make_gens_mod2(type_graph = type_graph)
        return(map_dbl(gens0, "end_count"))
}
# sample size functions
sample_size_gens_mod2 <- function (collect_gens, type_graph, sample_size)
{
        if (length(sample_size) == 1) {
                tip_sample_size = rep(sample_size, length(type_graph$tip_id))
                names(tip_sample_size) = type_graph$tip_id
        }
        else {
                assertthat::assert_that(length(sample_size) == length(type_graph$tip_id))
                assertthat::assert_that(all(type_graph$tip_id %in% names(sample_size)))
                tip_sample_size = sample_size[type_graph$tip_id]
        }
        tip_sample_size = pmin(tip_sample_size, calculate_type_counts_mod2(type_graph)[type_graph$tip_id])
        for (tip_id in type_graph$tip_id) {
                collect_gens[[tip_id]] = generate_sample_size_mod2(collect_gens[[tip_id]],
                                                                   tip_sample_size[tip_id])
        }
        for (node_id in backward_merge_sequence_mod2(type_graph)) {
                node_dau = type_graph$merge[[node_id]]
                collect_gens[[node_id]] = merge_sample_size_mod2(collect_gens[node_dau],
                                                                 collect_gens[[node_id]])
        }
        collect_gens
}
# Under construction !!
forward_merge_sequence_mod2 <- function (type_graph)
{
        node_seq = character()
        cur_nodes = type_graph$root_id
        while (length(node_seq) < length(type_graph$node_id)) {
                for (node_id in cur_nodes) {
                        if (!node_id %in% type_graph$tip_id) {
                                node_dau = type_graph$merge[[node_id]]
                                cur_nodes = cur_nodes[cur_nodes != node_id]
                                cur_nodes = c(cur_nodes, node_dau)
                                node_seq = append(node_seq, node_id)
                        }
                }
        }
        return(node_seq)
}
simulate_sc_muts_mod2 <- function (collect_gens, type_graph, mut_param)
{
        root_gens = collect_gens[[type_graph$root_id]]
        collect_gens[[type_graph$root_id]]$start_barcodes = init_cell_barcode(mut_param,
                                                                              num_cell = root_gens$sample_size[root_gens$num_gen+1])
        rownames(collect_gens[[type_graph$root_id]]$start_barcodes) = paste0("type_",
                                                                             type_graph$root_id, "_gen_", 0, "_", 1:nrow(collect_gens[[type_graph$root_id]]$start_barcodes))
        collect_gens[[type_graph$root_id]] = sample_sc_mut_history_gens_mod2(collect_gens[[type_graph$root_id]],
                                                                             mut_param = mut_param,
                                                                             target_time = type_graph$target_time)
        for (node_id in forward_merge_sequence_mod2(type_graph)) {
                node_dau = type_graph$merge[[node_id]]
                temp = split_barcodes_mod2(collect_gens[[node_id]],
                                           collect_gens[node_dau],
                                           mut_param = mut_param,
                                           target_time = type_graph$target_time)
                ###### TODO: fix output
                for (i in 1:length(node_dau)) {
                        collect_gens[[node_dau[i]]] = temp$x_ls[[i]]
                }
                collect_gens[[node_id]] = temp$y
        }
        collect_gens
}
sample_sc_mut_history_gens_mod2 <- function (x, mut_param, target_time = NA)
{
        message(x$cell_type)
        assertthat::assert_that(nrow(x$start_barcodes) == x$sample_size[x$num_gen + 1])
        b0 = x$start_barcodes
        b_history = list()
        e_history = list()
        n_history = list()
        if (nrow(b0) > 0) {
                if (x$num_gen >= 1) {
                        for (i in 1:x$num_gen) {
                                # print("mut gen")
                                # tt = proc.time()
                                b1 = sample_all_barcodes(b0, mut_param = mut_param,
                                                         start_time = x$start_time + x$double_time * (i - 1),
                                                         end_time = x$start_time + x$double_time * i)
                                # print(proc.time() - tt)
                                # print("remaining")
                                # tt = proc.time()
                                n0 = x$sample_size[x$num_gen + 2 - i] * 2 - x$sample_size[x$num_gen + 1 - i]
                                n1 = x$sample_size[x$num_gen + 2 - i] - n0 # n1 this doubling count
                                b1_dist = distribute_barcodes(b1, c(n0, n1))
                                b_history = append(b_history, list(b1_dist))
                                b0 = double_barcodes(b1_dist[[1]], b1_dist[[2]])
                                rownames(b0) = paste0("type_", x$cell_type,
                                                      "_gen_", i, "_", 1:nrow(b0))
                                edges = data.frame(in_node = c(rownames(b1_dist[[1]]),
                                                               rep(rownames(b1_dist[[2]]), each = 2)),
                                                   out_node = rownames(b0),
                                                   stringsAsFactors = F)
                                e_history = append(e_history, list(edges))
                                nodes = data.frame(id = c(rownames(b1_dist[[1]]),
                                                          rownames(b1_dist[[2]])),
                                                   degree = c(rep(2, nrow(b1_dist[[1]])),
                                                              rep(3, nrow(b1_dist[[2]]))),
                                                   type = x$cell_type,
                                                   gen = i - 1,
                                                   time = x$start_time + i * x$double_time,
                                                   double_time = x$double_time,
                                                   leaf = FALSE,
                                                   stringsAsFactors = F)
                                n_history = append(n_history, list(nodes))
                                # print(proc.time() - tt)
                        }
                }
        }
        b_end = b0
        m_end_time = ifelse(!x$active & !is.na(target_time), target_time,
                            x$end_time + x$double_time)
        b_end_mut = sample_all_barcodes(b_end, mut_param = mut_param,
                                        start_time = x$end_time,
                                        end_time = m_end_time)
        x$barcode_history = b_history
        x$edges_history = e_history
        x$nodes_history = n_history
        x$end_barcodes = b_end_mut
        if (x$active) {
                b_mode_dist = distribute_barcodes(b_end_mut, x$end_mode_sample_size)
                x$end_barcodes_mode = b_mode_dist
        }
        x
}
get_gen_from_id <- function (id_name) {
        out = map_chr(strsplit(id_name, "_"), function(x) {
                ifelse(length(x) == 5, x[4], "Unknown")
        })
        names(out) = id_name
        out
}
# m and m_vec are different things, should change naming
split_barcodes_mod2 <- function (y, x_ls, mut_param, target_time) {
        # message(y$cell_type)
        m = length(x_ls)

        # remap node names, the purpose is so that the commiting generation gets a label of the downstream states
        node_name_mapper = rownames(y$end_barcodes)
        names(node_name_mapper) = node_name_mapper
        # only remapping symmetric dividing cell names here
        for (i in 1:m) {
                old_node_names = rownames(y$end_barcodes_mode[[i]])
                new_node_names = paste0("type_", x_ls[[i]]$cell_type, "_gen_-1_",
                                        1:length(old_node_names))
                node_name_mapper[old_node_names] = new_node_names
        }
        # TODO: DELETE HERE
        # node_name_map = map2(y$end_barcodes_mode[1:m], x_ls, function(b, x) {
        #         old_node_names = rownames(b)
        #         new_node_names = paste0("type_", x$cell_type, "_gen_-1_",
        #                                 1:nrow(b))
        #         names(new_node_names) = old_node_names
        #         new_node_names
        # })
        # node_name_mapper = do.call(c, node_name_map)
        if (length(node_name_mapper) > 0) {
                if (y$num_gen == 0) {
                        y$transition_edges$out_node = node_name_mapper[y$transition_edges$out_node]
                        assertthat::assert_that(!anyNA(y$transition_edges$out_node))
                } else {
                        y$edges_history[[y$num_gen]]$out_node = node_name_mapper[y$edges_history[[y$num_gen]]$out_node]
                        assertthat::assert_that(!anyNA(y$edges_history[[y$num_gen]]$out_node))
                }
                for (i in 1:m) {
                        rownames(y$end_barcodes_mode[[i]]) = node_name_mapper[rownames(y$end_barcodes_mode[[i]])]
                }
        }
        # end remap node names

        # distributing and double barcodes
        commit_mode_ss_mat = do.call(cbind, map(x_ls, "commit_mode_sample_size"))
        s_vec = diag(commit_mode_ss_mat[1:m, ])
        m_vec = y$end_mode_sample_size[1:m]
        z_vec = s_vec - m_vec # number of doubles for symmetric

        n0_vec = 2 * m_vec - s_vec
        n1_vec = s_vec - m_vec # same as z_vec
        b_dist_ls = map(1:m, function(i) {
                distribute_barcodes(y$end_barcodes_mode[[i]],
                                    c(n0_vec[i], n1_vec[i]))

        })
        b_ls = map(1:m, function(i) {
                b_d = b_dist_ls[[i]]
                if ((nrow(b_d[[1]]) + nrow(b_d[[2]] > 0)) > 0) {
                        b_out = double_barcodes(b_d[[1]], b_d[[2]])
                        rownames(b_out) = paste0("type_", x_ls[[i]]$cell_type, "_gen_",
                                                 0, "_", 1:nrow(b_out))
                        b_out
                } else {
                        b_out = double_barcodes(b_d[[1]], b_d[[2]])
                        b_out
                }
                b_out
        })
        nodes_ls = map(1:m, function(i) {
                b_d = b_dist_ls[[i]]
                if ((nrow(b_d[[1]]) + nrow(b_d[[2]] > 0)) > 0) {
                        nodes = data.frame(id = c(rownames(b_d[[1]]), rownames(b_d[[2]])),
                                           degree = c(rep(2, nrow(b_d[[1]])), rep(3, nrow(b_d[[2]]))),
                                           type = x_ls[[i]]$cell_type, gen = -1, time = y$end_time + y$double_time,
                                           double_time = y$double_time, leaf = FALSE, stringsAsFactors = F)
                } else {
                        nodes = NULL
                }
                nodes
        })
        edges_ls = map(1:m, function(i) {
                b_d = b_dist_ls[[i]]
                if ((nrow(b_d[[1]]) + nrow(b_d[[2]] > 0)) > 0) {
                        edges = data.frame(in_node = c(rownames(b_d[[1]]),
                                                       rep(rownames(b_d[[2]]), each = 2)),
                                           out_node = rownames(b_ls[[i]]),
                                           stringsAsFactors = F)
                } else {
                        edges = NULL
                }
                edges
        })

        # for the asymmetric, keeping the old cell label, **effective commitment time one generation later**
        # for each downstream state, need to copy cells from two modes
        # a sampled asymmetric parent, needs to be distributed into three categories: A, B, and A + B
        # pair of downstream states for this mode
        xind_mode = map(1:choose(m, 2), function(j) {
                which(offspring_mat(m)[m + j, ] > 0)
        })
        n_ls = map(1:choose(m, 2), function(j) {
                s_a = commit_mode_ss_mat[(m + j), xind_mode[[j]]]
                m = y$end_mode_sample_size[(m + j)]
                z_a = sum(s_a) - m
                dist_vec = c(s_a - z_a, z_a) # sizes to be distributed into, in the order A, B and A + B
        })
        ba_dist_ls = map(1:choose(m, 2), function(j) {
                distribute_barcodes(y$end_barcodes_mode[[m + j]],
                                    n_ls[[j]])
        })
        ba_list = vector("list", m)
        nodes_asym_ls = vector("list", m)
        edges_asym_ls = vector("list", m)
        for (i in 1:m) {
                # looping over all modes, selecting corresponding cells
                b_mode_sel = map2(xind_mode, ba_dist_ls, function(indices, b_d) {
                        sel = which(indices == i)
                        if (length(sel) == 0) {
                                return(NULL)
                        } else {
                                assertthat::assert_that(sel == 1 | sel == 2)
                                return(list(b_d[[sel]],
                                            b_d[[3]]))
                        }
                })
                b_degree = do.call(c, map(b_mode_sel, function(x) {
                        if(!is.null(x)) {
                                c(rep(2, nrow(x[[1]])),
                                  rep(3, nrow(x[[2]])))
                        } else {
                                return(NULL)
                        }
                }))
                b_out = do.call(rbind, map(b_mode_sel, function(x) {
                        if(!is.null(x)) {
                                rbind(x[[1]], x[[2]])
                        } else {
                                return(NULL)
                        }
                }))
                if (nrow(b_out) == 0) {
                        ba_list[[i]] = NULL
                        nodes_asym_ls[[i]] = NULL
                        edges_asym_ls[[i]] = NULL
                        next
                }
                in_cell_names = rownames(b_out)
                out_cell_names = paste0("type_", x_ls[[i]]$cell_type, "_gen_",
                                        0, "_asym_", 1:nrow(b_out))
                rownames(b_out) = out_cell_names
                nodes = data.frame(id = in_cell_names,
                                   degree = b_degree,
                                   type = y$cell_type,
                                   gen = as.numeric(get_gen_from_id(in_cell_names)),
                                   time = y$end_time + y$double_time,
                                   double_time = y$double_time, leaf = FALSE, stringsAsFactors = F)
                edges = data.frame(in_node = in_cell_names,
                                   out_node = out_cell_names,
                                   stringsAsFactors = F)
                ba_list[[i]] = b_out
                nodes_asym_ls[[i]] = nodes
                edges_asym_ls[[i]] = edges
        }
        # combine barcodes, edges and nodes
        if (length(ba_list) > 0) {
                b_ls_all = map2(b_ls, ba_list, function(x ,y) {
                        rbind(x, y)
                })
                # nodes of degree 3 will be duplicated in both types
                nodes_ls_all = map2(nodes_ls, nodes_asym_ls, function(x, y) {
                        bind_rows(x, y)
                })
                edges_ls_all = map2(edges_ls, edges_asym_ls, function(x, y) {
                        bind_rows(x, y)
                })
        } else {
                b_ls_all = b_ls
                nodes_ls_all = nodes_ls
                edges_ls_all = edges_ls
        }
        for (i in 1:m) {
                assertthat::assert_that(nrow(b_ls_all[[i]]) == x_ls[[i]]$sample_size[x_ls[[i]]$num_gen + 1])
                x_ls[[i]]$start_barcodes = b_ls_all[[i]]
                x_ls[[i]]$transition_edges = edges_ls_all[[i]]
                x_ls[[i]]$transition_nodes = nodes_ls_all[[i]]
                x_ls[[i]] = sample_sc_mut_history_gens_mod2(x_ls[[i]],
                                                            mut_param = mut_param,
                                                            target_time = target_time)
        }
        # browser()
        return(list(y = y,
                    x_ls = x_ls))
}
construct_true_lineage_mod2 <- function (gens1, tip_id)
{
        nodes_df = collect_nodes(gens1, tip_id)
        nodes_df = unique(nodes_df)
        edges_df = collect_edges(gens1)
        edges_df = add_edge_length(edges_df, nodes_df)
        edges_nest = edges_df %>% mutate(in_node_dup = in_node) %>% nest(out = -in_node_dup)
        edges_list = edges_nest$out
        names(edges_list) = edges_nest$in_node_dup
        # edges_list = split(as.data.table(edges_df), edges_df$in_node)
        x_cell = edges_df[1, "in_node"]
        while (nrow(edges_list[[x_cell]]) == 1) {
                x_cell = edges_list[[x_cell]]$out_node
        }
        sp0 = list(edges_list[[x_cell]][1, ], edges_list[[x_cell]][2,
        ])
        n_tip = 0
        node_count = 0
        edges_simple_collect = list()
        while (TRUE) {
                p0 = lapply(sp0, get_next, edges_list = edges_list)
                edges_simple_collect = append(edges_simple_collect, list(do.call(rbind,
                                                                                 mapply(function(x, y) c(in_node = as.character(x$in_node),
                                                                                                         out_node = y$node, length = y$length), sp0, p0,
                                                                                        SIMPLIFY = F))))
                p0_node = p0[!sapply(p0, "[[", "tip")]
                tip_count = sum(sapply(p0, "[[", "tip"))
                node_count = length(p0_node)
                n_tip = n_tip + tip_count
                # message(paste0("merging_node: ", node_count))
                if (length(p0_node) > 0) {
                        sp1 = do.call(c, lapply(edges_list[as.character(sapply(p0_node,
                                                                               "[[", "node"))], function(x) list(x[1,
                                                                               ], x[2, ])))
                        sp0 = sp1
                }
                else {
                        assertthat::assert_that(sum(nodes_df$leaf) == n_tip)
                        (break)()
                }
        }
        edges_df_simple = data.frame(do.call(rbind, edges_simple_collect))
        edges_df_simple$length = as.numeric(as.character(edges_df_simple$length))
        edges_df_simple$in_node = as.character(edges_df_simple$in_node)
        edges_df_simple$out_node = as.character(edges_df_simple$out_node)
        tr_nw = edges_df_to_newick(edges_df_simple, nodes_df)
        tr_nw = paste0(tr_nw, nodes_df[1, "double_time"], ";")
        tr = ape::read.tree(text = tr_nw)
        tr
}
collect_nodes <- function(gens, tip_id) {
        # gens = gens1
        nodes_leaf = lapply(gens[tip_id], function(x) {
                if (nrow(x$end_barcodes) == 0) {
                        return(NULL)
                }
                nodes = data.frame(id = rownames(x$end_barcodes),
                                   degree = 1,
                                   gen = x$num_gen+1,
                                   type = x$cell_type,
                                   time = x$end_time,
                                   double_time = x$double_time,
                                   leaf = TRUE, stringsAsFactors = F)
                nodes
        })
        nodes_internal = lapply(gens, function(x) {
                rbind(x$transition_nodes,
                      do.call(rbind, x$nodes_history))
        })
        rbind(do.call(rbind, nodes_internal),
              do.call(rbind, nodes_leaf))
}
collect_edges <- function(gens) {
        out = do.call(rbind, lapply(gens, function(x) {
                rbind(x$transition_edges,
                      do.call(rbind, x$edges_history))
        }))
        out
}

simulate_sc_data_mod2 <- function (type_graph, mut_p, sample_size, delay = F) {
        gens0 = make_gens_mod2(type_graph)
        gens1 = sample_size_gens_mod2(gens0, type_graph, sample_size = sample_size)
        true_sampled_sizes = get_true_sample_size_mod2(gens1, type_graph, delay = delay)

        mut_param = do.call(make_mut_param_by_rate_rvec, mut_p)
        gens2  = simulate_sc_muts_mod2(gens1, type_graph, mut_param)

        tr = construct_true_lineage_mod2(gens2, type_graph$tip_id)
        x = get_barcodes(gens2, type_graph$tip_id)
        x = remove_allele_ver(x)

        out = list(tr = tr,
                   sc = x,
                   true_sampled_sizes = true_sampled_sizes)
        out
}
as.phylo_mod2.type_graph <- function (type_graph) {
        root_id = type_graph$root_id
        diff_edge = type_graph$edges
        diff_edge$in_node = as.character(diff_edge$in_node)
        diff_edge$out_node = as.character(diff_edge$out_node)

        # record tips as merging
        tips_collect = map(type_graph$tip_id, function(x) x)
        names(tips_collect) = type_graph$tip_id

        leaf_newick = lapply(type_graph$tip_id, function(x) {
                paste0(x, ":", diff_edge$length[diff_edge$out_node == x])
        })
        eligible_nodes = type_graph$tip_id
        eligible_newick = leaf_newick
        names(eligible_newick) = eligible_nodes
        while (length(eligible_nodes) > 1) {
                # print(eligible_nodes)
                p0_df = diff_edge[diff_edge$out_node %in% eligible_nodes, ]
                p0_df$out_newick = unlist(eligible_newick[p0_df$out_node])
                p0_df = nest(p0_df, out_data = -in_node)
                merging_ind = map2_lgl(p0_df$in_node, p0_df$out_data, function(x, dat) {
                        all(type_graph$merge[[x]] %in% dat$out_node)
                })
                p0_merging = p0_df[merging_ind, ]
                # make newick # shi gong xian chang

                # p0_count = table(p0_df$in_node)
                # p0_merging = as.numeric(names(p0_count)[p0_count == 2])
                # p0_df_merge = p0_df[p0_df$in_node %in% p0_merging, ]
                # d0_newick = eligible_newick[match(p0_df_merge$out_node,
                #                                   eligible_nodes)]
                tips_new = map(p0_merging$out_data, function(x) purrr::reduce(tips_collect[x$out_node], c))
                names(tips_new) = p0_merging$in_node
                tips_collect = c(tips_collect, tips_new)

                p0_merging$in_newick = map_chr(p0_merging$out_data, function(tb) {
                        paste0("(",
                               paste0(tb$out_newick, collapse = ","),
                               ")")
                })
                p0_merging$length = map_dbl(p0_merging$in_node, function(x) {
                        if (x == type_graph$root_id) {
                                return(type_graph$root_time)
                        }
                        diff_edge$length[diff_edge$out_node == x]
                })
                p0_merging$in_newick = paste0(p0_merging$in_node, p0_merging$in_newick)
                p0_merging$in_newick = paste0(p0_merging$in_newick, ":", p0_merging$length)

                merging_nodes = unlist(purrr::map(p0_merging$out_data, "out_node"))
                eligible_nodes1 = c(eligible_nodes[!eligible_nodes %in% merging_nodes],
                                    p0_merging$in_node)
                eligible_newick1 = c(eligible_newick[!eligible_nodes %in% merging_nodes],
                                     p0_merging$in_newick)
                names(eligible_newick1) = eligible_nodes1

                eligible_nodes = eligible_nodes1
                eligible_newick = eligible_newick1

                # p0_newick0 = lapply(split(d0_newick, p0_df_merge$in_node),
                #                     function(x) {
                #                             paste0("(", x[[1]], ",", x[[2]],
                #                                    ")")
                #                     })
                # p0_length = sapply(as.numeric(names(p0_newick0)), function(x) diff_edge[diff_edge$out_node ==
                #                                                                                 x, "length"])
                # p0_newick = mapply(function(z, x, y) paste0(x, z, ":",
                #                                             y), names(p0_newick0), p0_newick0, p0_length)
                # eligible_nodes1 = c(eligible_nodes[!eligible_nodes %in%
                #                                            p0_df_merge$out_node], as.numeric(names(p0_newick)))
                # eligible_newick1 = c(eligible_newick[!eligible_nodes %in%
                #                                              p0_df_merge$out_node], p0_newick)
                # eligible_nodes = eligible_nodes1
                # eligible_newick = eligible_newick1
        }
        tr = ape::read.tree(text = paste0(eligible_newick[[1]],
                                          # type_graph$diff_time[[root_id]],
                                          ";"))
        tr = name_nodes(tr)
        tr_tips = list_dd_and_tips_mod2(tr)$tips
        match_indices = match(map_chr(tr_tips, function(x) paste0(sort(x), collapse = ",")),
                              map_chr(tips_collect, function(x) paste0(sort(x), collapse = ",")))
        tr_node_mapper = names(tips_collect)[match_indices]
        names(tr_node_mapper) = names(tr_tips)
        tr$node.label = tr_node_mapper[tr$node.label]
        # tr$root.edge = type_graph$root_time
        tr
}
make_merge_sequences_mod2 <- function (edge_mat, dd_list, tip_id) {
        merge_sequence = list()
        eligible_nodes = tip_id
        while (length(eligible_nodes) > 1) {
                p0_df = edge_mat[edge_mat$out_node %in% eligible_nodes, ]
                p0_df = nest(p0_df, out_data = -in_node)
                p0_merging_ind = map2_lgl(p0_df$in_node, p0_df$out_data, function(x, tb) {
                        all(dd_list[[x]] %in% tb$out_node)
                })
                p0_df_merge = unnest(p0_df[p0_merging_ind, ], cols = out_data)
                cur_merge_list = split(p0_df_merge$out_node, p0_df_merge$in_node)
                merge_sequence = append(merge_sequence, list(cur_merge_list))
                eligible_nodes1 = c(eligible_nodes[!eligible_nodes %in%
                                                           p0_df_merge$out_node], unique(p0_df_merge$in_node))
                eligible_nodes = eligible_nodes1
        }
        merge_sequence
}
list_dd_and_tips_mod2 <- function (tr_r) {
        edge_mapper = character(max(tr_r$edge))
        edge_mapper[(1:Nnode(tr_r)) + Ntip(tr_r)] = tr_r$node.label
        edge_mapper[1:Ntip(tr_r)] = tr_r$tip.label
        edge_df = tibble(from = edge_mapper[tr_r$edge[, 1]], to = edge_mapper[tr_r$edge[,
                                                                                        2]])
        edge_df$length = tr_r$edge.length
        tr_dd_new = split(edge_df$to, edge_df$from)
        tr_dd_new = tr_dd_new[tr_r$node.label]
        m_seq = make_merge_sequences_mod2(dplyr::rename(edge_df, out_node = to,
                                                        in_node = from),
                                          dd_list = tr_dd_new,
                                          tip_id = tr_r$tip.label)
        tr_tip_list_new = map(tr_r$tip.label, function(x) x)
        names(tr_tip_list_new) = tr_r$tip.label
        for (x in m_seq) {
                y = map(x, function(dd) {
                        purrr::reduce(tr_tip_list_new[dd], c)
                })
                tr_tip_list_new = append(tr_tip_list_new, y)
        }
        tr_tip_list_new = tr_tip_list_new[tr_r$node.label]
        list(tips = tr_tip_list_new, dd = tr_dd_new)
}
# Deprecated, the diff time is off by one generation
# get_true_diff_time_mod2 <- function (type_graph) {
#         gens0 = make_gens_mod2(type_graph = type_graph)
#         diff_time = map_dbl(names(gens0),
#                             function(x) {
#                                     # message(x)
#                                     x_gens = gens0[[x]]
#                                     if (x_gens$active) {
#                                             if (x_gens$num_gen >= 1) {
#                                                     return(x_gens$start_time + (x_gens$num_gen -
#                                                                                         1) * x_gens$double_time)
#                                             }
#                                             else {
#                                                     assertthat::assert_that(x_gens$num_gen == 0)
#                                                     parent_node = get_parent_node(x, type_graph)
#                                                     x_parent_gens = gens0[[parent_node]]
#                                                     out = x_parent_gens$start_time + x_parent_gens$num_gen *
#                                                             x_parent_gens$double_time
#                                                     return(out)
#                                             }
#                                     }
#                                     else {
#                                             out = x_gens$start_time + x_gens$num_gen * x_gens$double_time
#                                             return(out)
#                                     }
#                             })
#         names(diff_time) = names(gens0)
#         diff_time[type_graph$tip_id] = type_graph$target_time
#         diff_time
# }
# delayed time if commitment is purely asymmetric
get_true_diff_time_mod3 <- function (type_graph, delay = F, target_tip_time = F) {
        gens0 = make_gens_mod2(type_graph = type_graph)
        if (delay) {
                diff_time = map_dbl(gens0, function(x) x$end_time + x$double_time)
                names(diff_time) = names(gens0)
        } else {
                diff_time = map_dbl(gens0, "end_time")
                names(diff_time) = names(gens0)
        }
        if (target_tip_time) {
                diff_time[type_graph$tip_id] = type_graph$target_time
        }
        diff_time
}
get_true_start_time_mod3 <- function (type_graph) {
        gens0 = make_gens_mod2(type_graph = type_graph)
        diff_time = map_dbl(gens0, "start_time")
        names(diff_time) = names(gens0)
        diff_time
}

#' this function returns the incoming size as well as the outgoing sizes
get_true_size_mod2 <- function (type_graph, delay = F) {
        gens0 = make_gens_mod2(type_graph = type_graph)
        # map_dbl(gens0, "num_gen")
        diff_count = map(names(gens0),
                         function(x) {
                                 # message(x)
                                 x_gens = gens0[[x]]
                                 if (x_gens$active) {
                                         if (delay) {
                                                 in_size = x_gens$end_count
                                         } else {
                                                 if (x_gens$num_gen >= 1) {
                                                         in_size = sum(x_gens$cell_mode_counts[x_gens$num_gen, c(1, 3)])
                                                 } else {
                                                         assertthat::assert_that(x_gens$num_gen == 0)
                                                         parent_node = get_parent_node(x, type_graph)
                                                         x_parent_gens = gens0[[parent_node]]
                                                         if (length(x_parent_gens$end_mode_counts) == 3) {
                                                                 in_size = sum(x_parent_gens$end_mode_counts[c(get_parent_split_index_mod2(x, parent_node, type_graph))]) +
                                                                         x_parent_gens$end_mode_counts[3] / 2
                                                         }
                                                         if (length(x_parent_gens$end_mode_counts) == 6) {
                                                                        parent_split_index = get_parent_split_index_mod2(x, parent_node, type_graph)
                                                                        n_symm = x_parent_gens$end_mode_counts[parent_split_index]
                                                                        if (parent_split_index == 1) {
                                                                                n_asymm = sum(x_parent_gens$end_mode_counts[c(4, 5)] /2)
                                                                        }
                                                                        if (parent_split_index == 2) {
                                                                                n_asymm = sum(x_parent_gens$end_mode_counts[c(4, 6)] /2)
                                                                        }
                                                                        if (parent_split_index == 3) {
                                                                                n_asymm = sum(x_parent_gens$end_mode_counts[c(5, 6)] /2)
                                                                        }
                                                                        in_size = n_symm + n_asymm
                                                         }
                                                 }
                                         }
                                         # out_size not affected by delay
                                         if (length(x_gens$end_mode_counts) == 3) {
                                                 out_size = x_gens$end_mode_counts[1:2] + x_gens$end_mode_counts[3]/2
                                         }
                                         if (length(x_gens$end_mode_counts) == 6) {
                                                 assertthat::assert_that(all(x_gens$end_mode_counts[4:6] == 0))
                                                 out_size = x_gens$end_mode_counts[1:3]
                                         }
                                 } else {
                                         in_size = sum(x_gens$cell_mode_counts[x_gens$num_gen + 1, c(1, 3)])
                                         out_size = rep(NA, 2)
                                 }
                                 out = c(in_size, out_size)
                                 names(out) = c(x, type_graph$merge[[x]])
                                 out
                         })
        names(diff_count) = names(gens0)
        diff_count
}
get_true_sample_size_mod2 <- function(gens1, type_graph, delay = F) {
        assertthat::assert_that(!is.null(gens1[[1]]$sample_size_mode))
        diff_count = map(names(gens1),
                         function(x) {
                                 # message(x)
                                 x_gens = gens1[[x]]
                                 if (x_gens$active) {
                                         if (delay) {
                                                 in_size = x_gens$sample_size[1]
                                                 dd_types = type_graph$merge[[x_gens$cell_type]]
                                                 if (length(dd_types) == 2) {
                                                         x_gens_d1 = gens1[[dd_types[1]]]
                                                         x_gens_d2 = gens1[[dd_types[2]]]
                                                         out_size = c(x_gens_d1$sample_size[x_gens_d1$num_gen+1],
                                                                      x_gens_d2$sample_size[x_gens_d2$num_gen+1])
                                                 }
                                                 if (length(dd_types) == 3) {
                                                         x_gens_d1 = gens1[[dd_types[1]]]
                                                         x_gens_d2 = gens1[[dd_types[2]]]
                                                         x_gens_d3 = gens1[[dd_types[3]]]
                                                         out_size = c(x_gens_d1$sample_size[x_gens_d1$num_gen+1],
                                                                      x_gens_d2$sample_size[x_gens_d2$num_gen+1],
                                                                      x_gens_d3$sample_size[x_gens_d3$num_gen+1])
                                                 }
                                         } else {
                                                 if (x_gens$num_gen >= 1) {
                                                         in_size = x_gens$sample_size[2]
                                                 } else {
                                                         assertthat::assert_that(x_gens$num_gen == 0)
                                                         parent_node = get_parent_node(x, type_graph)
                                                         x_parent_gens = gens1[[parent_node]]
                                                         parent_split_index = get_parent_split_index_mod2(x, parent_node, type_graph)
                                                         if (length(x_parent_gens$end_mode_sample_size) == 3) {
                                                                 in_size = x_parent_gens$end_mode_sample_size[parent_split_index]+
                                                                         x_parent_gens$end_mode_sample_size[3]/2
                                                         }
                                                         if (length(x_parent_gens$end_mode_sample_size) == 6) {
                                                                 n_symm = x_parent_gens$end_mode_sample_size[parent_split_index]
                                                                 if (parent_split_index == 1) {
                                                                         n_asymm = sum(x_parent_gens$end_mode_sample_size[c(4, 5)] /2)
                                                                 }
                                                                 if (parent_split_index == 2) {
                                                                         n_asymm = sum(x_parent_gens$end_mode_sample_size[c(4, 6)] /2)
                                                                 }
                                                                 if (parent_split_index == 3) {
                                                                         n_asymm = sum(x_parent_gens$end_mode_sample_size[c(5, 6)] /2)
                                                                 }
                                                                 in_size = n_symm + n_asymm
                                                         }
                                                 }
                                                 if (length(x_gens$end_mode_sample_size) == 3) {
                                                         out_size = x_gens$end_mode_sample_size[1:2] + x_gens$end_mode_counts[3]/2
                                                 }
                                                 if (length(x_gens$end_mode_sample_size) == 6) {
                                                         assertthat::assert_that(all(x_gens$end_mode_sample_size[4:6] == 0))
                                                         out_size = x_gens$end_mode_sample_size[1:3]
                                                 }
                                         }
                                 }
                                 else {
                                         in_size = sum(x_gens$sample_size[1])
                                         out_size = rep(NA, 2)
                                 }
                                 out = c(in_size, out_size)
                                 names(out) = c(x, type_graph$merge[[x]])
                                 out
                         })
        names(diff_count) = names(gens1)
        diff_count
}

get_parent_split_index_mod2 <- function(node, parent_node, type_graph) {
        which(as.character(type_graph$merge[[parent_node]]) == node)
}

generate_edge_df_mod2 <- function (type_graph) {
        diff_edge = bind_rows(purrr::map(names(type_graph$merge), function(x) {
                tibble(in_node = x,
                       out_node = type_graph$merge[[x]])
        }))

        # diff_edge = do.call(rbind, lapply(1:nrow(type_graph$merge),
        #                                   function(i) {
        #                                           data.frame(in_node = i, out_node = type_graph$merge[i,
        #                                           ])
        #                                   }))
        diff_time_edge = list()
        gens0 = make_gens_mod2(type_graph)
        type_graph$edges = diff_edge

        gens0_time = get_true_diff_time_mod3(type_graph)

        diff_edge$length = gens0_time[diff_edge$out_node] -
                gens0_time[diff_edge$in_node]
        assertthat::assert_that(all(diff_edge$length > 0))
        type_graph$edges = diff_edge
        type_graph$root_time = gens0_time[[type_graph$root_id]]
        type_graph
}
get_mode_counts <- function(cell_type, count, type_graph) {
        mode_prob = type_graph$diff_mode_probs[[cell_type]]
        if (length(mode_prob) == 3) {
                assertthat::assert_that(count >= 8)
                if (mode_prob[3] == 0) {
                        mode_min_val = c(4, 4, 0)
                } else {
                        mode_min_val = c(0, 0, 8)
                }
        } else if (length(mode_prob) == 6) {
                assertthat::assert_that(count >= 12)
                if (all(mode_prob[4:6] == 0)) {
                        mode_min_val = c(4, 4, 4, 0, 0, 0)
                } else {
                        mode_min_val = c(0, 0, 0, 4, 4, 4)
                }
        } else {
                mode_min_val = rep(0, length(mode_prob))
        }
        mode_counts = round_frac(count, mode_prob, mode_min_val)
        mode_counts
}
round_frac <- function(total, prob, min_val = rep(0, length(prob))) {
        assertthat::assert_that(length(min_val) == length(prob))
        assertthat::assert_that(sum(min_val) <= total)
        n_ind = which.max(prob)[1]
        rounded = round(total * prob[-n_ind])
        rounded = pmax(rounded, min_val[-n_ind])
        out = numeric(length(prob))
        out[-n_ind] = rounded
        out[n_ind] = total - sum(rounded)
        assertthat::assert_that(out[n_ind] >= min_val[n_ind])
        names(out) = names(prob)
        out
}
plot_type_graph_mod2 <- function (type_graph, node_col_mapper = white_col_mapper, node_label_mapper = NULL,
                                  effective_time = T, show_node_text = T)
{
        gens0 = make_gens_mod2(type_graph)
        type_ig = ape::as.igraph.phylo(as.phylo_mod2.type_graph(type_graph))
        V(type_ig)$name = str_replace_all(V(type_ig)$name, ":",
                                          "")
        V(type_ig)$Type = V(type_ig)$name
        V(type_ig)$Time = get_true_diff_time_mod3(type_graph)[V(type_ig)$name]
        V(type_ig)$counts = map_dbl(get_true_size_mod2(type_graph)[V(type_ig)$name], 1)
        double_time = unlist(type_graph$lifetime)
        V(type_ig)$double_time = double_time[match(V(type_ig)$name,
                                                   names(double_time))]
        V(type_ig)$label = paste0(ifelse(V(type_ig)$counts > 1000,
                                         format(V(type_ig)$counts, digits = 2), V(type_ig)$counts))
        if (!is.null(node_label_mapper)) {
                V(type_ig)$node_label = node_label_mapper(V(type_ig)$name)
        }
        else {
                V(type_ig)$node_label = V(type_ig)$name
        }
        node_col = node_col_mapper(sort(V(type_ig)$name))
        type_ig = type_ig %>% add_vertices(1, name = "Founder",
                                           node_label = "F", Time = 0) %>% add_edges(c("Founder",
                                                                                       as.character(type_graph$root_id)))
        V(type_ig)[["Founder"]]$label = paste0("count: ",
                                               type_graph$founder_size)
        node_col = c(node_col, "#FFFFFF")
        names(node_col)[length(node_col)] = "Founder"
        for (x_n in names(type_graph$trans_sel)) {
                x = type_graph$trans_sel[[x_n]]
                x_gens = gens0[[x_n]]
                node_out = x_n
                if (x_n == type_graph$root_id) {
                        node_in = "Founder"
                }
                else {
                        node_in = type_graph$edges$in_node[type_graph$edges$out_node ==
                                                                   x_n]
                }
                for (j in 1:length(x)) {
                        y = x[[j]]
                        trans_name = paste0(x_n, " unsampled ", j)
                        type_ig = type_ig - edge(paste0(node_in, "|",
                                                        node_out))
                        trans_time = ifelse(effective_time, x_gens$start_time +
                                                    x_gens$double_time * (which(x_gens$cell_mode_counts[,
                                                                                                        2] > 0)[j] - 3), x_gens$start_time + x_gens$double_time *
                                                    (which(x_gens$cell_mode_counts[, 2] > 0)[j] -
                                                             2))
                        type_ig = type_ig %>% add_vertices(1, name = trans_name,
                                                           Time = trans_time) %>% add_edges(c(node_in, trans_name)) %>%
                                add_edges(c(trans_name, node_out))
                        V(type_ig)[[trans_name]]$label = NA
                        node_col = c(node_col, "#FFFFFF")
                        names(node_col)[length(node_col)] = trans_name
                        node_in = trans_name
                }
        }
        V(type_ig)$node_label[is.na(V(type_ig)$node_label)] = V(type_ig)$name[is.na(V(type_ig)$node_label)]
        g_out = ggraph(type_ig, layout = "dendrogram", height = Time) +
                geom_edge_bend(aes(width = factor(1), fontface = "plain"),
                               arrow = arrow(type = "closed", length = unit(3,
                                                                            "pt"), angle = 45), start_cap = ggraph::circle(5,
                                                                                                                           "pt"), end_cap = ggraph::circle(5, "pt")) +
                geom_node_label(aes(label = node_label, fill = name,
                                    size = factor(1)), label.padding = unit(6, "pt"))
        if (show_node_text) {
                g_out = g_out + geom_node_text(aes(label = label), nudge_x = -0.2,
                                               nudge_y = 0.3)
        }
        g_out = g_out + scale_fill_manual(values = node_col) + ylim(c(type_graph$target_time,
                                                                      0)) + ylab("Time") + scale_edge_width_manual(values = 0.5,
                                                                                                                   guide = "none") + scale_size_manual(values = 4,
                                                                                                                                                       guide = "none") + theme(axis.line.y = element_line(),
                                                                                                                                                                               axis.ticks.y = element_line(), axis.text.y = element_text(size = 12),
                                                                                                                                                                               axis.title = element_text(size = 12), panel.background = element_rect(fill = NA,
                                                                                                                                                                                                                                                     color = NA, size = 10), legend.position = "none",
                                                                                                                                                                               axis.title.x = element_blank(), text = element_text(size = 12))
        g_out
}
plot_type_graph_clean_mod2 <- function(type_graph, node_col_mapper = white_col_mapper,
                                       node_label_mapper = NULL, effective_time = T, add_founder = F,
                                       show_node_text = F, node_size = 3) {
        gens0 = make_gens_mod2(type_graph)
        type_ig = ape::as.igraph.phylo(as.phylo_mod2.type_graph(type_graph))
        V(type_ig)$name = str_replace_all(V(type_ig)$name, ":", "")
        V(type_ig)$Type = V(type_ig)$name
        V(type_ig)$Time = get_true_diff_time_mod3(type_graph)[V(type_ig)$name]
        V(type_ig)$Time[V(type_ig)$name %in% type_graph$tip_id] = type_graph$target_time
        # node_frac = do.call(rbind,
        #                     lapply(type_graph$node_id, function(id)
        #                             data.frame(
        #                                     id = as.character(type_graph$merge[id,]),
        #                                     fraction = type_graph$diff_mode_probs[[id]][1:2]
        #                             )))
        # V(type_ig)$fraction = node_frac[match(V(type_ig)$name, node_frac$id), "fraction"]
        V(type_ig)$counts = map_dbl(get_true_size_mod2(type_graph)[V(type_ig)$name], 1)
        double_time = unlist(type_graph$lifetime)
        V(type_ig)$double_time = double_time[match(V(type_ig)$name,
                                                   names(double_time))]
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
        V(type_ig)$tip = ifelse(!V(type_ig)$name %in% type_graph$tip_id, "Node", "Tip")

        if (!is.null(node_label_mapper)) {
                V(type_ig)$node_label = node_label_mapper(V(type_ig)$name)
        } else {
                V(type_ig)$node_label =V(type_ig)$name
        }
        node_col = node_col_mapper(sort(V(type_ig)$name))

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

        if (add_founder) {
                type_ig = type_ig %>% add_vertices(1,
                                                   name = "Founder",
                                                   Type = "Founder",
                                                   node_label = "F",
                                                   Time = 0) %>%
                        add_edges(c("Founder", as.character(type_graph$root_id)))
                V(type_ig)[["Founder"]]$label = paste0("count: ", type_graph$founder_size)
                # node_col = c(node_col, node_col[type_graph$root_id])
                node_col = c(node_col, "#FFFFFF")
                names(node_col)[length(node_col)] = "Founder"
        }

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
                geom_node_point(aes(color = Type, shape = tip), size = node_size)
        if (show_node_text) {
                g_out = g_out +
                        geom_node_text(aes(label = node_label), nudge_x = - 0.2, nudge_y = 0.3)
        }
        g_out = g_out +
                scale_color_manual(values = node_col) +
                ylim(c(type_graph$target_time, 0)) + ylab('Time') +
                scale_edge_width_manual(values = 0.5, guide = "none") +
                scale_size_manual(values = 3, guide = "none") +
                scale_shape_manual(values = c(15, 17)) +
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

plot_type_graph_clean_mod2 <- function(type_graph, node_col_mapper = white_col_mapper,
                                       node_label_mapper = NULL, effective_time = T, add_founder = F,
                                       show_node_text = F, node_size = 3,
                                       show_node_counts = F, ylim = c(type_graph$target_time, 0)) {
        gens0 = make_gens_mod2(type_graph)
        type_ig = ape::as.igraph.phylo(as.phylo_mod2.type_graph(type_graph))
        V(type_ig)$name = str_replace_all(V(type_ig)$name, ":", "")
        V(type_ig)$Type = V(type_ig)$name
        V(type_ig)$Time = get_true_diff_time_mod3(type_graph)[V(type_ig)$name]
        V(type_ig)$Time[V(type_ig)$name %in% type_graph$tip_id] = type_graph$target_time
        # node_frac = do.call(rbind,
        #                     lapply(type_graph$node_id, function(id)
        #                             data.frame(
        #                                     id = as.character(type_graph$merge[id,]),
        #                                     fraction = type_graph$diff_mode_probs[[id]][1:2]
        #                             )))
        # V(type_ig)$fraction = node_frac[match(V(type_ig)$name, node_frac$id), "fraction"]
        V(type_ig)$counts = map_dbl(get_true_size_mod2(type_graph)[V(type_ig)$name], 1)
        double_time = unlist(type_graph$lifetime)
        V(type_ig)$double_time = double_time[match(V(type_ig)$name,
                                                   names(double_time))]
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
        V(type_ig)$count_label = paste0(ifelse(V(type_ig)$counts > 1000,
                                               format(V(type_ig)$counts, digits = 2),
                                               V(type_ig)$counts))
        V(type_ig)$tip = ifelse(!V(type_ig)$name %in% type_graph$tip_id, "Node", "Tip")

        if (!is.null(node_label_mapper)) {
                V(type_ig)$node_label = node_label_mapper(V(type_ig)$name)
        } else {
                V(type_ig)$node_label =V(type_ig)$name
        }
        node_col = node_col_mapper(sort(V(type_ig)$name))

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

        if (add_founder) {
                type_ig = type_ig %>% add_vertices(1,
                                                   name = "Founder",
                                                   Type = "Founder",
                                                   node_label = "F",
                                                   Time = 0) %>%
                        add_edges(c("Founder", as.character(type_graph$root_id)))
                V(type_ig)[["Founder"]]$label = paste0("count: ", type_graph$founder_size)
                # node_col = c(node_col, node_col[type_graph$root_id])
                node_col = c(node_col, "#FFFFFF")
                names(node_col)[length(node_col)] = "Founder"
        }

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
                geom_node_point(aes(color = Type, shape = tip), size = node_size)
        if (show_node_text) {
                g_out = g_out +
                        geom_node_text(aes(label = node_label), nudge_x = - 0.2, nudge_y = 0.3)
        }
        if (show_node_counts) {
                g_out = g_out +
                        geom_node_text(aes(label = count_label), nudge_x = - 0.2, nudge_y = -0.3)
        }
        g_out = g_out +
                scale_color_manual(values = node_col) +
                ylim(ylim) + ylab('Time') +
                scale_edge_width_manual(values = 0.5, guide = "none") +
                scale_size_manual(values = 3, guide = "none") +
                scale_shape_manual(values = c(15, 17)) +
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

get_type_from_id <- function(id_name) {
        out = map_chr(strsplit(id_name, "_"), function(x) {
                ifelse(length(x) %in% c(5, 6), x[2], "Unknown")
        })
        names(out) = id_name
        out
}

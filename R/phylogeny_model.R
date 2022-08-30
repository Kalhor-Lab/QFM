# STEP1: Generate count based on type graph
#' Given a single cell type, generate the number of generations it will divide until it reaches target time or differentiate
num_gen_single_type <- function(cell_type, rtime, type_graph) {
        # get relevant info from graph
        diff_time = type_graph$diff_time[[cell_type]]
        life_duration = type_graph$lifetime[[cell_type]]
        diff_indicator = F
        if (type_graph$target_time <= diff_time) {
                # This situation could be case 2 or case 3, behavior decided by \
                # target time
                # NOTE: do not count the division that happens after target time
                diff_indicator = F
                num_gen = floor((type_graph$target_time - rtime)/life_duration)
        } else {
                # In this case the last division during differentiation may be
                # after target time, which means still no differentiation
                # expected number of generations regardless of target time
                num_gen = ceiling((diff_time - rtime)/life_duration)
                if (rtime + num_gen * life_duration > type_graph$target_time) {
                        num_gen = num_gen - 1
                        diff_indicator = F
                } else {
                        diff_indicator = T
                }
        }
        return(list(num_gen = num_gen,
                    diff = diff_indicator))
}
# deprecated
# get_mode_counts <- function(cell_type, count, type_graph) {
#         assertthat::assert_that(count >= 2)
#         mode_prob = type_graph$diff_mode_probs[[cell_type]]
#         if (length(mode_prob) == 3) {
#                 mode_min_val = c(1, 1, 0)
#         } else if (length(mode_prob) == 6) {
#                 mode_min_val = c(1, 1, 1, 0, 0, 0)
#         } else {
#                 mode_min_val = rep(0, length(mode_prob))
#         }
#         mode_counts = round_frac(count, mode_prob, mode_min_val)
#         mode_counts
# }
distribute_mode_counts <- function(cell_type, mode_counts, type_graph) {
        # no mode 4 and 5
        assertthat::assert_that(all(mode_counts[4:5] == 0))
        daughter_mode_counts = mode_counts * mode_offspring
        colnames(daughter_mode_counts) =
                c(cell_type, type_graph$merge[cell_type, ])
        daughter_mode_counts
}
make_cell_type_gens <- function(cell_type,
                                start_time,
                                start_count,
                                start_mode_counts,
                                type_graph) {
        gen_state = num_gen_single_type(cell_type, start_time, type_graph)
        life_duration = type_graph$lifetime[[cell_type]]
        if (all(!is.na(start_mode_counts))) {
                assertthat::assert_that(sum(start_mode_counts) == start_count)
        }
        # not counting the division as the cell differentiates, as spliting is needed
        eff_gen = ifelse(gen_state$diff,
                         gen_state$num_gen - 1,
                         gen_state$num_gen)
        # TODO: how to fix too many cells?
        cell_counts = start_count * 2^(0:eff_gen)
        # transition selections
        if (cell_type %in% names(type_graph$trans_sel)) {
                # transitions within the same type
                type_trans = type_graph$trans_sel[[cell_type]]
                # transition type and symmetric vs asymetric fraction
                type_trans_time = sapply(type_trans, "[[", "time")
                type_trans_frac = sapply(type_trans, "[[", "fraction") # columns are transitions
                cell_gen_time = start_time + life_duration * (1:eff_gen)
                # the generation at which transition happens
                type_trans_gen = sapply(type_trans_time, function(rt) max(which(rt > cell_gen_time)))+1
                assertthat::assert_that(all(type_trans_gen <= eff_gen))
                # number of cells transitioning via three different
                cell_mode_counts = matrix(NA, nrow = (eff_gen+1), ncol = 3)
                cell_mode_counts[1, 1] = start_count
                cell_mode_counts[1, 2:3] = 0
                for (i in 1:eff_gen) {
                        cell_count_last = cell_mode_counts[i, 1] * 2 + cell_mode_counts[i, 3]
                        if (i %in% type_trans_gen) {
                                m = type_trans_frac[, which(type_trans_gen==i)]
                                assertthat::assert_that(length(m) == 2)
                                assertthat::assert_that(sum(m) <= 1)
                                cell_mode_counts[i+1, ] = round_frac(cell_count_last, c(m[1], 1 - sum(m), m[2]))
                        } else {
                                cell_mode_counts[i+1, 1] = cell_count_last
                                cell_mode_counts[i+1, 2:3] = 0
                        }
                }
                end_count = cell_mode_counts[eff_gen+1, 1] + cell_mode_counts[eff_gen+1, 3]
        } else {
                cell_mode_counts = NULL
                end_count = cell_counts[eff_gen+1]
        }
        end_mode_counts = get_mode_counts(cell_type, end_count, type_graph)
        out = list(cell_type = cell_type,
                   start_time = start_time,
                   start_count = start_count,
                   start_mode_counts = start_mode_counts,
                   num_gen = eff_gen,
                   double_time = life_duration,
                   cell_counts = cell_counts,
                   cell_mode_counts = cell_mode_counts,
                   end_time = start_time + life_duration * eff_gen,
                   end_count = end_count,
                   end_mode_counts = end_mode_counts,
                   active = gen_state$diff
        )
        class(out) = "cell_type_gen"
        return(out)
}
next_cell_type_gens <- function(x, type_graph) {
        assertthat::assert_that(x$active)
        # counts are doubled here via mode counts
        daughter_mode_counts = distribute_mode_counts(x$cell_type,
                                                      x$end_mode_counts,
                                                      type_graph)
        cg1 = make_cell_type_gens(colnames(daughter_mode_counts)[2],
                                  start_time = x$end_time + x$double_time,
                                  start_count = sum(daughter_mode_counts[, 2]),
                                  start_mode_counts = daughter_mode_counts[, 2],
                                  type_graph = type_graph)
        cg2 = make_cell_type_gens(colnames(daughter_mode_counts)[3],
                                  start_time = x$end_time + x$double_time,
                                  start_count = sum(daughter_mode_counts[, 3]),
                                  start_mode_counts = daughter_mode_counts[, 3],
                                  type_graph = type_graph)
        return(list(cg1, cg2))
}
generate_sample_size <- function(x, sample_size) {
        # message(x$cell_type)
        if (class(x$cell_mode_counts) == "matrix") {
                cell_mode_sample_size = matrix(NA, nrow = (x$num_gen+1), ncol = 3)
                cell_mode_sample_size[, 2] = 0
                cell_mode_sample_size[x$num_gen+1, c(1, 3)] = extraDistr::rmvhyper(nn = 1, k = sample_size, x$cell_mode_counts[x$num_gen+1, c(1, 3)])[1, ]
                for (i in x$num_gen:1) {
                        # assymetric needs testing
                        s1 = cell_mode_sample_size[i+1, 1]
                        s1a = cell_mode_sample_size[i+1, 3]
                        # TODO: this needs to be double checked
                        z1 = sample_doublet_branch(ceiling(x$cell_mode_counts[i+1, 1]/2), s1)
                        m1 = s1 - z1
                        cell_mode_sample_size[i, 1] = m1
                        cell_mode_sample_size[i, 3] = s1a
                }
                x$sample_size_mode = cell_mode_sample_size
                # ignore asymmetric divisions for now
                x$sample_size = x$sample_size_mode[, 1]
        } else {
                out = numeric(x$num_gen+1)
                out[x$num_gen+1] = sample_size
                if (x$num_gen > 0) {
                        k0 =  sample_size
                        n0 = x$end_count
                        for (i in x$num_gen:1) {
                                n0 = n0/2
                                z = sample_doublet_branch(n0, k0)
                                k1 = k0 - z
                                k0 = k1
                                out[i] = k0
                        }
                }
                x$sample_size = out
        }
        if (all(!is.na(x$start_mode_counts))) {
                mode_x = x$start_mode_counts[1:2]
                assertthat::assert_that(mode_x[1] == 0 | mode_x[2] == 0)
                mode_x = mode_x[mode_x > 0]
                mode_y = x$start_mode_counts[3]
                # split into two start modes
                mode_sample_size = extraDistr::rmvhyper(nn = 1, k = x$sample_size[1], c(mode_x, mode_y))[1, ]
                x$start_mode_sample_size = mode_sample_size
        } else {
                x$start_mode_sample_size = NA
        }
        x
}
merge_sample_size <- function(x1, x2, y) {
        # message(paste0(x1$cell_type, x2$cell_type, y$cell_type))
        assertthat::assert_that(y$end_mode_counts[1] * 2 == x1$start_mode_counts[1])
        assertthat::assert_that(y$end_mode_counts[3] == x1$start_mode_counts[3])
        assertthat::assert_that(y$end_mode_counts[2] * 2 == x2$start_mode_counts[2])
        assertthat::assert_that(y$end_mode_counts[3] == x2$start_mode_counts[3])
        s1 = x1$start_mode_sample_size[1]
        s2 = x2$start_mode_sample_size[1]
        s1a = x1$start_mode_sample_size[2]
        s2a = x2$start_mode_sample_size[2]
        z1 = sample_doublet_branch(y$end_mode_counts[1], s1)
        z2 = sample_doublet_branch(y$end_mode_counts[2], s2)
        z3 = sample_doublet_twotype(n = y$end_mode_counts[3], s1a, s2a)
        m1 = s1 - z1
        m2 = s2 - z2
        m3 = s1a + s2a - z3
        y$end_mode_sample_size = c(m1, m2, m3)
        y = generate_sample_size(y, sample_size = sum(c(m1, m2, m3)))
        y
}
double_mut_counts <- function(y) {
        # first element of y is singleton mut_counts
        lapply(y, function(x) x[[1]] + x[[2]]*2)
}
get_new_muts <- function(c0, c1) {
        # fetch new mutations comparing c1 to c0
        assertthat::assert_that(all(sapply(c0, sum) == sapply(c1, sum)))
        lapply(1:length(c0), function(j) {
                mut0 = names(c0[[j]])
                mut1 = names(c1[[j]])
                mut1[!mut1 %in% mut0]
        })
}
generate_gen_cell_id <- function(type, gen, n) {
        paste0(type, "_", gen, 1:n)
}
sample_mut_history_gens <- function(x, mut_param, target_time=NA) {
        # target_time for stopping the doubling at some time
        # generate history of mutations from "start_mut_counts"
        assertthat::assert_that(all(sapply(x$start_mut_counts, sum) == x$sample_size[1]))
        c0 = x$start_mut_counts
        c_history = list()
        m_history = list()
        if (x$num_gen >= 1) {
                for (i in 1:x$num_gen) {
                        c1 = sample_all_mut_counts(c0,
                                                   mut_param = mut_param,
                                                   start_time = x$start_time + x$double_time * (i - 1),
                                                   end_time = x$start_time + x$double_time * i)
                        # distribute into singleton vs otherwise
                        # singleton
                        n0 = x$sample_size[i]*2 - x$sample_size[i+1]
                        # non-singleton
                        n1 = x$sample_size[i] - n0
                        # distributed and mutated counts
                        if (sum(c1[[1]])!=(n0+n1)) {
                                print(x)
                        }
                        c1_dist = distribute_all_mut_counts(c1, n0, n1)
                        c_history = append(c_history, list(c1_dist))
                        # m_history = append(m_history, list(get_new_muts(c0, c1)))

                        c0 = double_mut_counts(c1_dist)
                }
        }
        c_end = c0
        # if ((target_time - x$end_time) < 0) {
        #         print(target_time - x$end_time)
        # }
        m_end_time = ifelse(!x$active & !is.na(target_time),
                            target_time,
                            x$end_time + x$double_time)
        c_end_mut = sample_all_mut_counts(c_end,
                                          mut_param = mut_param,
                                          start_time = x$end_time,
                                          end_time = m_end_time)
        # m_history = append(m_history, list(get_new_muts(c_end, c_end_mut)))
        x$mut_counts_history = c_history
        # x$mut_history = m_history
        x$end_mut_counts = c_end_mut
        if (x$active) {
                # distribute between modes
                c_mode_dist1 = distribute_all_mut_counts(c_end_mut,
                                                         x$end_mode_sample_size[1],
                                                         sum(x$end_mode_sample_size[2:3]))
                c_mode_dist2 = distribute_all_mut_counts(lapply(c_mode_dist1, "[[", 2),
                                                         x$end_mode_sample_size[2],
                                                         x$end_mode_sample_size[3])
                c_mode = mapply(function(x, y) c(x[1], y), c_mode_dist1, c_mode_dist2, SIMPLIFY = F)
                x$end_mut_counts_mode = c_mode
        }
        x
}
split_mut_counts <- function(y, x1, x2, mut_param, target_time) {
        # child 1
        # singletons
        c1_n0_m1 = 2 * y$end_mode_sample_size[1] - x1$start_mode_sample_size[1]
        c1_n1_m1 = y$end_mode_sample_size[1] - c1_n0_m1
        c1_m1_dist = distribute_all_mut_counts(lapply(y$end_mut_counts_mode, "[[", 1),
                                               c1_n0_m1,
                                               c1_n1_m1)
        c1_m1 = double_mut_counts(c1_m1_dist)
        # TODO: need to add singletons for assymetric mode
        c1 = double_mut_counts(mapply(function(x, z) c(list(x), z[3]), c1_m1, y$end_mut_counts_mode, SIMPLIFY = F))

        # child 2
        c2_n0_m1 = 2 * y$end_mode_sample_size[2] - x2$start_mode_sample_size[1]
        c2_n1_m1 = y$end_mode_sample_size[2] - c2_n0_m1
        c2_m1_dist = distribute_all_mut_counts(lapply(y$end_mut_counts_mode, "[[", 2),
                                               c2_n0_m1,
                                               c2_n1_m1)
        c2_m1 = double_mut_counts(c2_m1_dist)
        # TODO: need to add singletons for assymetric mode
        c2 = double_mut_counts(mapply(function(x, z) c(list(x), z[3]), c2_m1, y$end_mut_counts_mode, SIMPLIFY = F))
        # browser()

        # print(c1)
        # print(x1)
        assertthat::assert_that(all(sapply(c1, sum) == x1$sample_size[1]))
        assertthat::assert_that(all(sapply(c2, sum) == x2$sample_size[1]))
        x1$start_mut_counts = c1
        x2$start_mut_counts = c2
        x1 = sample_mut_history_gens(x1, mut_param = mut_param, target_time = target_time)
        x2 = sample_mut_history_gens(x2, mut_param = mut_param, target_time = target_time)
        return(list(x1, x2))
}
extract_allele_matrix <- function(y, slot = "end_mut_counts") {
        if (!is.null(slot)) {
                cells_states_mut = lapply(y, "[[", slot)
        } else {
                cells_states_mut = y
        }
        num_elements = sapply(cells_states_mut, length)
        assertthat::assert_that(is_zero_range(num_elements))
        m_counts_list = lapply(1:num_elements[1], function(j) {
                m_alleles = unique(do.call(c, lapply(cells_states_mut, function(x) names(x[[j]]))))
                # m_alleles = m_alleles[order(as.numeric(m_alleles))]
                zz = do.call(rbind, lapply(cells_states_mut, function(x) {
                        m_counts = x[[j]]
                        out = m_counts[match(m_alleles, names(m_counts))]
                        names(out) = m_alleles
                        out
                }))
                rownames(zz) = names(y)
                zz[is.na(zz)] = 0
                zz
        })
        m_counts_list
}
filter_allele_matrix <- function(mut_mat, frac_cut = 0.05) {
        frac = mut_mat/rowSums(mut_mat)
        frac[is.nan(frac)] = 0
        fil = apply(frac > frac_cut, 2, any)
        mut_mat[, fil, drop = F]
}
strip_recur_ver <- function(x) {
        sapply(strsplit(x, "\\("), "[[", 1)
}
collapse_reucr_ver <- function(mut_mat) {
        allele_name = strip_recur_ver(colnames(mut_mat))
        out = t(mat_colsum_fac(t(mut_mat), allele_name))
        rownames(out) = rownames(mut_mat)
        out
}
mat_colsum_fac = function (mat, fac)
{
        if (class(fac) != "factor")
                fac <- factor(fac)
        rown <- length(levels(fac))
        coln <- dim(mat)[2]
        out <- matrix(NA, rown, coln)
        ind <- as.numeric(fac)
        for (i in 1:rown) {
                out[i, ] <- colSums(mat[ind == i, , drop = F], na.rm = T)
        }
        rownames(out) <- levels(fac)
        return(out)
}
get_mutated_frac <- function(mut_mat) {
        apply(mut_mat, 1, function(x) {
                sum(x[names(x) != "0"])/sum(x)
        })
}
make_gens <- function(type_graph) {
        init_gens = make_cell_type_gens(type_graph$root_id,
                                        start_time = 0,
                                        start_count = type_graph$founder_size,
                                        start_mode_counts = NA,
                                        type_graph)
        collect_gens = list()
        collect_gens = append(collect_gens, list(init_gens))
        active_gens = list(init_gens)
        while(length(active_gens) > 0) {
                next_gens = unlist(lapply(active_gens, function(x) next_cell_type_gens(x, type_graph = type_graph)),
                                   recursive = F)
                collect_gens = append(collect_gens, next_gens)
                active_ind = sapply(next_gens, "[[", "active")
                active_gens = next_gens[active_ind]
        }
        names(collect_gens) = sapply(collect_gens, "[[", "cell_type")
        collect_gens
}
sample_size_gens <- function(collect_gens, type_graph, sample_size = 5000) {
        if (length(sample_size) == 1) {
                tip_sample_size = rep(sample_size, length(type_graph$tip_id))
                names(tip_sample_size) = type_graph$tip_id
        } else {
                assertthat::assert_that(length(sample_size) == length(type_graph$tip_id))
                assertthat::assert_that(all(type_graph$tip_id %in% names(sample_size)))
                tip_sample_size = sample_size[type_graph$tip_id]
        }
        tip_sample_size = pmin(tip_sample_size, calculate_type_counts_new(type_graph)[type_graph$tip_id])
        for (tip_id in type_graph$tip_id) {
                collect_gens[[tip_id]] = generate_sample_size(collect_gens[[tip_id]],
                                                              tip_sample_size[tip_id])
        }

        for (node_id in backward_merge_sequence(type_graph)) {
                node_dau = as.character(type_graph$merge[node_id, ])
                collect_gens[[node_id]] =
                        merge_sample_size(
                                collect_gens[[node_dau[1]]],
                                collect_gens[[node_dau[2]]],
                                collect_gens[[node_id]]
                        )
        }
        collect_gens
}
simulate_muts <- function(collect_gens, type_graph, mut_param) {
        collect_gens[[type_graph$root_id]]$start_mut_counts = init_mut_counts_list(mut_param, num_cells = type_graph$founder_size)
        collect_gens[[type_graph$root_id]] = sample_mut_history_gens(collect_gens[[type_graph$root_id]],
                                                                     mut_param = mut_param,
                                                                     target_time = type_graph$target_time)
        for (node_id in forward_merge_sequence(type_graph)) {
                node_dau = as.character(type_graph$merge[node_id, ])
                temp = split_mut_counts(collect_gens[[node_id]],
                                        collect_gens[[node_dau[1]]],
                                        collect_gens[[node_dau[2]]],
                                        mut_param = mut_param,
                                        target_time = type_graph$target_time)
                collect_gens[[node_dau[1]]] = temp[[1]]
                collect_gens[[node_dau[2]]] = temp[[2]]
        }
        collect_gens
}
simulate_sc_data <- function(type_graph, mut_p, sample_size, somatic = F) {
        gens0 = make_gens(type_graph = type_graph)
        gens1 = sample_size_gens(gens0, type_graph, sample_size = sample_size)
        if (!somatic) {
                mut_param = do.call(make_mut_param_by_rate_rvec, mut_p)
        }
        # message("simluating sc..")
        gens2 = simulate_sc_muts(gens1,
                                 type_graph = type_graph,
                                 mut_param = mut_param,
                                 somatic = somatic)
        # message("constructing true tree..")
        tr = construct_true_lineage(gens2, type_graph$tip_id)
        if (somatic) {
                x = NULL
        } else {
                x = get_barcodes(gens2, type_graph$tip_id)
                x = remove_allele_ver(x)
        }
        # this is the pre-split field size
        sampled_field_size = map_dbl(c(type_graph$node_id, type_graph$tip_id), function(x) {
                s = gens1[[x]]$sample_size
                if (gens1[[x]]$active) {
                        if(length(s) >= 2) {
                                return(s[length(s) - 1])
                        } else {
                                assertthat::assert_that(length(s) == 1)
                                parent_node = get_parent_node(x, type_graph)
                                s_out = gens1[[parent_node]]$end_mode_sample_size[get_parent_split_index(x, parent_node, type_graph)]
                                return(s_out)
                        }
                } else {
                        return(s[length(s)])
                }
        })
        names(sampled_field_size) = c(type_graph$node_id, type_graph$tip_id)
        sampled_field_split = map(type_graph$node_id, function(node) {
                out = gens1[[node]]$end_mode_sample_size[1:2]
                names(out) = as.character(type_graph$merge[node, ])
                out
        })
        # below are the split field size
        # sampled_field_size = map_dbl(c(type_graph$node_id, type_graph$tip_id), function(x) {
        #         s = gens1[[x]]$sample_size
        #         return(s[length(s)])
        # })
        out = list(tr = tr,
                   sc = x,
                   sample_size = sample_size,
                   sampled_field_size = sampled_field_size,
                   sampled_field_split = sampled_field_split)
        return(out)
}

# plots for bulk testing
# plot_type_graph(type_graph0, type_col_mapper)
# plot_type_graph(type_graph1, type_col_mapper)

# testing
# mut_p = list(mut_rate = rep(1/4, 10),
#              recur_prob = 1.,
#              recur_vec = extraDistr::rdirichlet(1, alpha = rep(0.15, 200)))
# type_graph = type_graph0
# mut_param = do.call(make_mut_param_by_rate, mut_p)
#
# gens0 = make_gens(type_graph = type_graph)
# gens1 = sample_size_gens(gens0, type_graph, sample_size = 5000)
# gens1 = simulate_muts(gens1, type_graph, mut_param = mut_param)
#
# cell_history = lapply(gens0, function(x) {
#         t0 = x$start_time + (0:x$num_gen)*x$double_time
#         t1 = t0 + x$double_time
#         cell_count = x$start_count * 2^(0:x$num_gen)
#         list(cell_type = x$cell_type,
#              t0 = t0,
#              t1 = t1,
#              count = cell_count)
# })
# birth_events = do.call(rbind, lapply(cell_history, function(x) cbind(x$t0, x$count)))
# birth_events = data.frame(birth_events[order(birth_events[, 1]), ])
# death_events = do.call(rbind, lapply(cell_history, function(x) {
#         out = cbind(x$t1, x$count)
#         if (x$cell_type %in% type_graph$tip_id) {
#                 out = out[1:(nrow(out)-1), ]
#         }
#         out
# }))
# death_events = data.frame(death_events[order(death_events[, 1]), ])
# death_events$X2 = -death_events$X2
# collapse_events = rbind(birth_events, death_events)
# collapse_events$X1 = round(collapse_events$X1, digits = 3)
# collapse_events = collapse_events %>% group_by(X1) %>% summarise(increment = sum(X2))
# collapse_events$count = cumsum(collapse_events$increment)
# collapse_events = collapse_events[1:(nrow(collapse_events)-1), ]
# plot(collapse_events$X1, log2(collapse_events$count))

# produce_phy <- function(type_graph, mut_p) {
#         mut_param = do.call(make_mut_param_by_rate, mut_p)
#         gens0 = make_gens(type_graph = type_graph)
#         gens1 = sample_size_gens(gens0, type_graph, sample_size = 5000)
#
#         gens1 = simulate_muts(gens1, type_graph, mut_param = mut_param)
#         tip_mat = extract_allele_matrix(gens1[type_graph$tip_id])
#         # tip_mat_fil = lapply(tip_mat, filter_allele_matrix, frac_cut = 0.01)
#         # tip_mat_scaled = lapply(tip_mat_fil, function(x) clr(normalize_row(x, m = 5000)))
#         # collapsed version
#         tip_mat_collapsed = lapply(tip_mat, collapse_reucr_ver)
#         tip_mat_collapsed_fil = lapply(tip_mat_collapsed, filter_allele_matrix, frac_cut = 0.01)
#         tip_mat_collapsed_scaled = lapply(tip_mat_collapsed_fil, function(x) clr(normalize_row(x, m = 5000)))
#         tip_mat_collasped_scaled_all = do.call(cbind, tip_mat_collapsed_scaled)
#
#         phy = upgma(dist(tip_mat_collasped_scaled_all))
#         return(list(mat = tip_mat_collasped_scaled_all, phy = phy))
#
#         # phy$tip.label
#         # phy$edge
#         # as.igraph(phy)
#         #
#         # plot(phy, use.edge.length = T)
#         # edgelabels(format(phy$edge.length, digit = 3), 1:10, cex = 1.5)
# }
# check_phy_and_get_edge_len <- function(phy) {
#         phy$edge
# }
# correct_ind = sapply(dat0_phy, function(x) all.equal(as.phylo(type_graph0), x$phy, use.edge.length = F))
# plot(dat0_phy[[4]]$phy)
# Heatmap(dat0_phy[[4]]$mat, name = "CLR", show_column_names = F, show_column_dend = F, left_annotation = type_annt)
#
# dat0_edge = lapply(dat0_phy[correct_ind], function(x) x$phy$edge)
# reference_edge = dat0_edge[[1]]
# edge_ident = which(sapply(dat0_edge, function(y) all((reference_edge - y)==0)))
# edge_len0 = do.call(rbind, lapply(dat0_phy[correct_ind][edge_ident], function(x) x$phy$edge.length))
#
# correct_ind1 = sapply(dat1_phy, function(x) all.equal(as.phylo(type_graph1), x$phy, use.edge.length = F))
# dat1_edge = lapply(dat1_phy[correct_ind1], function(x) x$phy$edge)
# edge_ident1 = which(sapply(dat1_edge, function(y) all((reference_edge - y)==0)))
# edge_len1 = do.call(rbind, lapply(dat1_phy[correct_ind1][edge_ident1], function(x) x$phy$edge.length))
#
# phy_node_map = c(-1:-6, c(5, 3, 4, 2, 1))
# edge_name = apply(matrix(phy_node_map[reference_edge], ncol = 2), 1, function(x) paste0(x[1], "_", x[2]))
# edge_df = phy$edge
#
# dat0_phy = lapply(1:100, function(i) {
#         produce_phy(type_graph0, mut_p)
# })
# dat1_phy = lapply(1:100, function(i) produce_phy(type_graph1, mut_p))
#
#
# # dat_test0 = do.call(rbind, lapply(1:50, function(i) produce_edge_len(type_graph0, mut_p)))
# # dat_test1 = do.call(rbind, lapply(1:50, function(i) produce_edge_len(type_graph1, mut_p)))
#
#
# col_fun = colorRamp2(seq(0, 12, length = 3), c("blue", "#EEEEEE", "red"))
# colnames(edge_len0) = edge_name
# colnames(edge_len1) = edge_name
# Heatmap(edge_len0, col = col_fun, cluster_columns = F,
#         name = "Distance", show_row_dend = F)
# Heatmap(edge_len1, col = col_fun, cluster_columns = F, name = "Distance", show_row_dend = F)
#
# par(mfrow = c(1, 3))
# for (i in c(5, 6, 10)) {
#         boxplot(cbind(G0 = edge_len0[, i], G1 = edge_len1[, i]), main = edge_name[i])
# }
#
# boxplot(cbind(G0 = edge_len0[, 5], G1 = edge_len1[, 5]))


# make_mut_history <- function(gens, type_graph, mut_param) {
#         mut_history_map = lapply(1:length(mut_param), function(element_id){
#                 out = do.call(rbind, lapply(c(forward_merge_sequence(type_graph), type_graph$tip_id), function(node_id) {
#                         mut_gen_seq = lapply(gens1[[node_id]]$mut_history, "[[", element_id)
#                         do.call(rbind, lapply(1:length(mut_gen_seq), function(g) {
#                                 if (length(mut_gen_seq[[g]]) > 0) {
#                                         data.frame(Type = node_id,
#                                                    Mut = as.character(mut_gen_seq[[g]]),
#                                                    Gen = factor(g))
#                                 }
#                         }))
#                 }))
#                 out = rbind(data.frame(Type = "0", Mut = "0", Gen = 0), out)
#                 out$element_id = element_id
#                 out
#         })
#         allele_his = lapply(1:length(mut_param), function(j) {
#                 his_map = mut_history_map[[j]]
#                 mat_names = colnames(tip_allele_mat_fil[[j]])
#                 # assertthat::assert_that(length(mat_names) == nrow(his_map))
#                 assertthat::assert_that(all(mat_names %in% his_map$Mut))
#                 his_map[match(mat_names, his_map$Mut), ]
#         })
#         allele_his_df = do.call(rbind, allele_his)
# }
# allele_his_df = make_mut_history(gens1, type_graph0, mp0)
#
# library(RColorBrewer)
# node_names = c(0, type_graph$node_id, type_graph$tip_id)
# col_vec = type_col_mapper(node_names)
# col_vec["0"] = "#878787"
# names(col_vec) = node_names
# gen_col = brewer.pal(8, "YlGnBu")
# names(gen_col) = 0:(length(gen_col)-1)
# allele_his_map = lapply(1:length(tip_mat_scaled), function(j) {
#         his_df = subset(allele_his_df, allele_his_df$element_id == j)
#         his_df[match(colnames(tip_mat_scaled[[j]]), his_df$Mut), ]
# })
#
#
# allele_annt = HeatmapAnnotation(df=do.call(rbind, allele_his_map)[, c("Type", "Gen")],
#                                 col = list(Type = col_vec, Gen = gen_col))
# type_annt = rowAnnotation(df = data.frame(Type = rownames(tip_mat_scaled[[1]])),
#                           col = list(Type = type_col_mapper(rownames(tip_mat_scaled[[1]]))))
# Heatmap(do.call(cbind, tip_mat_scaled), top_annotation = allele_annt,
#         left_annotation = type_annt,
#         column_order = order(allele_his_df$Type, allele_his_df$Gen),
#         show_column_names = F, name = "clr")
# Heatmap(do.call(cbind, tip_mat_scaled), top_annotation = allele_annt,
#         left_annotation = type_annt, name = "clr", show_column_names = F)
# Heatmap(do.call(cbind, tip_mat_collasped_scaled),
#         left_annotation = type_annt, show_column_names = F,
#         name = "clr")
# emb_spacer_mat = do.call(cbind, lapply(emb_spacer$E12p5_t3689_emb1$spacer_count, function(x) clr(normalize_row(filter_allele_matrix(x, frac = 0.01), m = 5000))))
# Heatmap(emb_spacer_mat, show_column_names = F, show_column_dend = F, name = "clr")




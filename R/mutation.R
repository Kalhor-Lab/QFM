#' Initialize the global mut_element counter for non-recurring alleles
initialize_mut_element_counter <- function(mut_element, element_dim = 1) {
        assertthat::assert_that(mut_element$max_depth == 1)
        return(rep(1, element_dim))
}
#' Constructor for mutation parameter for a single element
#' Only supports max_depth of 1
#' @param recur chance of recurring
#' @export
make_mut_element <- function(element_id, mut_rate, active_time, recur_prob = 0,
                             recur_vec = NULL, max_depth = 1) {
        assertthat::assert_that(recur_prob >= 0 & recur_prob <= 1)
        out = list(element_id = element_id,
                   rate = mut_rate,
                   active_time = active_time,
                   recur_prob = recur_prob,
                   recur_vec = recur_vec,
                   max_depth = max_depth)
        class(out) = "mut_element"
        if (recur_prob < 1) {
                if(!exists('global_mut_counter', envir = .GlobalEnv)) {
                        global_mut_counter <<- list()
                }
                global_mut_counter[[element_id]] <<- initialize_mut_element_counter(out)
        }
        if (recur_prob > 0) {
                assertthat::assert_that(!is.null(recur_vec))
                # assertthat::assert_that(sum(recur_vec) == 1)
                if(!exists('global_mut_recur_counter', envir = .GlobalEnv)) {
                        global_mut_recur_counter <<- list()
                }
                global_mut_recur_counter[[element_id]] <<- initialize_mut_element_counter(out, length(recur_vec))
                names(global_mut_recur_counter[[element_id]]) <<- paste0("R-", 1:length(recur_vec))
        }
        out
}
#' Make a mut_param by a sequence of mutation rates
#' @export
make_mut_param_by_rate <- function(mut_rate,
                                   active_time,
                                   recur_prob = 0,
                                   recur_vec = NULL) {
        # reset global
        global_mut_counter <<- list()
        global_mut_recur_counter <<- list()
        out = lapply(1:length(mut_rate), function(j) {
                make_mut_element(element_id = j,
                                 mut_rate = mut_rate[j],
                                 active_time = active_time,
                                 recur_prob = recur_prob,
                                 recur_vec = recur_vec)
        })
        class(out) = "mut_param"
        out
}
#' Make a mut_param by a sequence of mutation rates and multinomial allele probabilities
#' @export
make_mut_param_by_rate_rvec <- function(mut_rate,
                                        active_time,
                                        recur_prob_vec,
                                        recur_vec_list) {
        # reset global
        global_mut_counter <<- list()
        global_mut_recur_counter <<- list()
        out = lapply(1:length(mut_rate), function(j) {
                assertthat::assert_that(length(recur_vec_list) == length(mut_rate))
                make_mut_element(element_id = j,
                                 mut_rate = mut_rate[j],
                                 active_time = active_time,
                                 recur_prob = recur_prob_vec[j],
                                 recur_vec = recur_vec_list[[j]])
        })
        class(out) = "mut_param"
        out
}
#' Initialize unmutated matrix for single cell model
#' @export
initialize_mut_state <- function(mut_param, num_cells) {
        return(matrix("0", num_cells, length(mut_param)))
}
#' Initialize mutant allele count list for bulk counting
#' @export
init_mut_counts_list <- function(mut_param, num_cells) {
        assertthat::assert_that(all(mut_param$max_depth == 1))
        out = lapply(1:length(mut_param), function(j) {
                c("0"=num_cells)
        })
        class(out) = "mut_counts_list"
        return(out)
}
#' Gives the counts for mut_counts_list
#' @export
mut_counts_total <- function(mut_counts_list) {
        m_total = sapply(mut_counts_list, sum)
        assertthat::assert_that(is_zero_range(m_total))
        return(m_total[1])
}

# generate_mut <- function(m, mut_element) {
#         #TODO: needs to fix
#         m_depth = str_count(m, ",")
#         assertthat::assert_that(m_depth == 1)
#         if (m_depth < mut_param$max_depth[j]) {
#                 # old version
#                 # m_letter = global_mut_counter[[j]][[m]]
#                 # global_mut_counter[[j]][[m]] <<- global_mut_counter[[j]][[m]] + 1
#                 # new_m = paste0(m, ",", m_letter)
#                 # global_mut_counter[[j]][[new_m]] <<- 0
#                 # return(new_m)
#                 # end old version
#                 recur_ind = sample(c(T, F), size = 1, replace = T,
#                                    prob = mut_param[[j]]$recur_prob)
#                 if (recur_ind) {
#                         return(generate_recur_mut_ids(1, j, mut_param))
#                 }
#         } else {
#                 return(m)
#         }
# }

# ------------------------------------------------------------------------------
# Functions for generating new alleles given that a number of mutation has occurred
#' Generate a batch of non recurring depth-1 mutante allele
#' @export
generate_nonrecur_mut_ids <- function(n, mut_element) {
        if (n <= 0) return(NULL)
        j = mut_element$element_id
        new_ids = (global_mut_counter[[j]]):(global_mut_counter[[j]]+n-1)
        global_mut_counter[[j]] <<- global_mut_counter[[j]] + n
        new_ids
}
generate_recur_mut_ver_ids <- function(n, element_id, recur_id) {
        if (n <= 0) return(NULL)
        j = element_id
        cur = global_mut_recur_counter[[j]][recur_id]
        new_ids = (cur):(cur+n-1)
        global_mut_recur_counter[[j]][recur_id] <<- global_mut_recur_counter[[j]][recur_id] + n
        paste0("(", new_ids, ")")
}
#' Generate a batch of recurring depth-1 mutant allele
#' @export
generate_recur_mut_ids <- function(n, mut_element, version = TRUE) {
        if (n <= 0) return(NULL)
        recur_ids = paste0("R-", 1:length(mut_element$recur_vec))
        out = sample(recur_ids, size = n, replace = T, prob = mut_element$recur_vec)
        if (version) {
                for (x in recur_ids) {
                        idx = which(out == x)
                        out[idx] = paste0(x, generate_recur_mut_ver_ids(length(idx),
                                                                        element_id = mut_element$element_id,
                                                                        recur_id = x))
                }
        }
        out
}
#' Generate a batch of depth-1 mutante allele
#' @export
generate_mut_ids = function(n, mut_element, recur_ver = TRUE) {
        # sample if recurring allele is generated
        recur_prob = mut_element$recur_prob
        recur_ind = sample(c(T, F), size = n,
                           prob = c(recur_prob, 1 - recur_prob),
                           replace = T)
        out = character(n)
        out[recur_ind] = generate_recur_mut_ids(sum(recur_ind), mut_element, version = recur_ver)
        out[!recur_ind] = generate_nonrecur_mut_ids(sum(!recur_ind), mut_element)
        out
}
sample_mut_num <- function(n, mut_element, duration) {
        mut_rate = mut_element$rate
        # mut_prob = 1 - dpois(0, lambda = duration * mut_rate)
        mut_prob = 1 - exp(- (duration * mut_rate))
        mut_num = rbinom(n = 1, size = n, prob = mut_prob)
        mut_num
}
# This functrion is deprecated, somatic mutation conditional on at least one event,
# direct forward simulation is used instead.
# sample_somatic_mut_num <- function(n, mut_element, duration, unmut_duration, t_total) {
#         mut_rate = mut_element$rate
#         mut_prob = (exp(-mut_rate * unmut_duration) - exp(-mut_rate * (unmut_duration + mut_duration))) /
#                 (exp(-mut_rate * unmut_duration) - exp(-mut_rate * (t_total)))
#         mut_num = rbinom(n = 1, size = n, prob = mut_prob)
#         mut_num
# }

#' Draw one sample of the mutation counts after some time
#' @param life_duration length of the generation
#' @param double_counts whether the counts of the alleles should be doubled as a reulst of division
#' @export
sample_mut_counts <- function(mut_counts,
                              mut_element,
                              duration) {
        # do not collapse recurring here
        assertthat::assert_that(mut_element$max_depth == 1)
        # only the unmutated can mutate, sample number of mutation events
        mut_num = sample_mut_num(mut_counts[1],
                                 mut_element,
                                 duration)
        mut_counts[1] = mut_counts[1] - mut_num

        new_muts = generate_mut_ids(mut_num, mut_element)
        assertthat::assert_that(all(new_muts != "0"))
        assertthat::assert_that(anyDuplicated(new_muts) == 0)
        assertthat::assert_that(!any(new_muts %in% names(mut_counts)))
        new_muts_count = rep(1, length(new_muts))
        names(new_muts_count) = new_muts
        mut_counts = c(mut_counts, new_muts_count)
        mut_counts
}
sample_all_barcodes <- function(barcodes, mut_param, start_time, end_time) {
        assertthat::assert_that(ncol(barcodes) == length(mut_param))
        for (j in 1:length(mut_param)) {
                u_indices = which(barcodes[, j] == "0")
                if (length(u_indices) > 0) {
                        mut_num = sample_mut_num(length(u_indices),
                                                 mut_param[[j]],
                                                 duration = get_mut_duration(mut_param[[j]]$active_time,
                                                                             start_time, end_time))
                        m_indices = sample(u_indices, size = mut_num, replace = F)
                        barcodes[m_indices, j] = generate_mut_ids(length(m_indices),
                                                                  mut_param[[j]])
                }
        }
        barcodes
}
sample_all_barcodes_v1 <- function(barcodes, mut_param, start_time, end_time, parallel = T) {
        assertthat::assert_that(ncol(barcodes) == length(mut_param))
        out = do.call(cbind, future_map(1:length(mut_param), function(j) {
                mut_el = mut_param[[j]]
                barcode_vec = barcodes[, j]
                u_indices = which(barcode_vec == "0")
                if (length(u_indices) > 0) {
                        mut_num = sample_mut_num(length(u_indices),
                                                 mut_el,
                                                 duration = get_mut_duration(mut_el$active_time,
                                                                             start_time, end_time))
                        m_indices = sample(u_indices, size = mut_num, replace = F)
                        barcode_vec[m_indices] = generate_mut_ids(length(m_indices),
                                                                  mut_el)
                }
                barcode_vec
        }, .options = furrr_options(seed = T)))
        out
}
#' get overlaps of intervals using trick on IRanges
#' TODO: dependes on IRanges
get_mut_duration <-  function(active_time, start_time, end_time) {
        mul = 1e5
        ir = IRanges::IRanges(start = sapply(active_time, "[[", 1) * mul,
                              end = sapply(active_time, "[[", 2) * mul)
        sum(IRanges::width(IRanges::intersect(IRanges::IRanges(start = start_time * mul,
                                                               end = end_time * mul),
                                              ir))-1)/mul
}
#' deprecated function, used for somatic mutation
#' get the total amount of time passed that the site could have mutated
# get_unmut_duration <- function(active_time, start_time) {
#         mul = 1e5
#         ir = IRanges::IRanges(start = sapply(active_time, "[[", 1) * mul,
#                               end = sapply(active_time, "[[", 2) * mul)
#         sum(IRanges::width(IRanges::intersect(IRanges::IRanges(start = 0,
#                                                                end = start_time * mul),
#                                               ir))-1)/mul
# }
get_t_total <- function(active_time) {
        sum(map_dbl(active_time, function(x) {
                out = x[2] - x[1]
                assertthat::assert_that(out > 0)
                out
        }))
}
# old version of sampling somatic barcode
# sample_somatic_barcodes_old <- function(barcodes, mut_param, start_time, end_time) {
#         assertthat::assert_that(ncol(barcodes) == length(mut_param))
#         # for each element
#         for (j in 1:length(mut_param)) {
#                 # cells capable of mutating
#                 u_indices = which(barcodes[, j] == "0")
#                 if (length(u_indices) > 0) {
#                         # number of cells mutated
#                         mut_num = sample_somatic_mut_num(length(u_indices),
#                                                          mut_param[[j]],
#                                                          unmut_duration = get_unmut_duration(mut_param[[j]]$active_time,
#                                                                                              start_time),
#                                                          duration = get_mut_duration(mut_param[[j]]$active_time,
#                                                                                      start_time, end_time),
#                                                          t_total = get_t_total(mut_param[[j]]$active_time))
#                         mut_num = sample_mut_num(length(u_indices),
#                                                  mut_param[[j]],
#                                                  duration = get_mut_duration(mut_param[[j]]$active_time,
#                                                                              start_time, end_time))
#                         m_indices = sample(u_indices, size = mut_num, replace = F)
#                         barcodes[m_indices, j] = generate_mut_ids(length(m_indices),
#                                                                   mut_param[[j]])
#                 }
#         }
#         barcodes
# }


#' Draw one sample of mutation counts for all elements
#' @export
sample_all_mut_counts <- function(mut_counts_list, mut_param, start_time, end_time,
                                  double_counts = F) {
        assertthat::assert_that(length(mut_counts_list) == length(mut_param))
        out = lapply(1:length(mut_counts_list), function(j) {
                m_counts = sample_mut_counts(mut_counts_list[[j]],
                                             mut_param[[j]],
                                             get_mut_duration(mut_param[[j]]$active_time,
                                                              start_time, end_time))
                if (double_counts) {
                        m_counts = m_counts * 2
                }
                return(m_counts)
        })
        out
}
#' Draw one sample of the mutation counts after a number of generatioins
#' @param num_gen the number of generation
#' @export
sample_all_mut_counts_gen <- function(mut_counts_list, mut_param,
                                  life_duration, num_gen,
                                  double_counts = T) {
        assertthat::assert_that(num_gen >= 1)
        for (i in 1:num_gen){
                mut_counts_list = sample_all_mut_counts(mut_counts_list,
                                                        mut_param,
                                                        life_duration,
                                                        double_counts = double_counts)
        }
        return(mut_counts_list)
}
#' Draw a subsample of mutant allele counts
#' @export
subsample_mut_counts <- function(mut_counts, sample_size) {
        if (length(mut_counts) == 1){
                assertthat::assert_that(sample_size <= mut_counts)
                sampled_counts = sample_size
                names(sampled_counts) = names(mut_counts)
                return(sampled_counts)
        } else {
                sampled_counts = extraDistr::rmvhyper(1, n = mut_counts, k = sample_size)[1, ]
                names(sampled_counts) = names(mut_counts)
                return(sampled_counts)
        }
}
#' Draw subsample of allele counts for all elements
#' @export
subsample_all_mut_counts <- function(mut_counts_list, sample_size) {
        return(lapply(mut_counts_list, function(x) subsample_mut_counts(x, sample_size)))
}
#' Distribute mutant allele counts between two groups
#' @importFrom extraDistr rmvhyper
#' @export
distribute_mut_counts <- function(mut_counts, n1, n2) {
        assertthat::assert_that(sum(mut_counts) == (n1 + n2),
                                msg = paste0("n1: ", n1, "\nn2: ", n2, "\nmut_counts: ", sum(mut_counts)))
        if (length(mut_counts) == 1){
                mut_counts1 = n1
        } else {
                mut_counts1 = extraDistr::rmvhyper(1, n = mut_counts, k = n1)[1, ]
        }
        mut_counts2 = mut_counts - mut_counts1
        names(mut_counts1) = names(mut_counts)
        names(mut_counts2) = names(mut_counts)
        return(list(mut_counts1,
                    mut_counts2))
}
#' Distribute allele counts for all elements
#' @export
distribute_all_mut_counts <- function(mut_counts_list, n1, n2) {
        return(lapply(mut_counts_list, function(x) distribute_mut_counts(x, n1, n2)))
}
#' Generate mutant allele for single cell mutation states
#' @export
mutate_barcodes <- function(mut_state, life_duration, mut_param) {
        assertthat::assert_that(length(mut_param) == ncol(mut_state))
        out_mut_state = do.call(cbind, lapply(1:length(mut_param), function(j) {
                in_state = mut_state[, j]
                unmut_idx = which(in_state == "0")
                mut_time = rexp(n = length(unmut_idx),
                                rate = mut_param[[j]]$rate)
                idx = which((life_duration[unmut_idx] - mut_time) > 0)
                mut_idx = unmut_idx[idx]
                out_state = in_state
                out_state[mut_idx] = generate_mut_ids(length(mut_idx),
                                                      mut_param[[j]])
                return(out_state)
                # OLD multi-depth version BELOW
                # mut_time = rexp(n = length(life_duration),
                #                 rate = mut_param[[j]]$rate)
                # time_remain = life_duration - mut_time
                # out_state = in_state
                # print(time_remain)
                # while (any(time_remain > 0)) {
                #         mut_ind = time_remain > 0
                #         for (i in which(mut_ind)) {
                #                 out_state[i] = generate_mut(out_state[i], j, mut_param)
                #         }
                #         time_remain[mut_ind] = time_remain[mut_ind] - rexp(n = sum(mut_ind), rate = mut_param[[j]]$rate)
                # }

        }))
        out_mut_state
}
#' Mutate the barcodes for cells object
#' @export
mutate_cells <- function(cells, mut_param) {
        cells$mut_state = mutate_barcodes(cells$mut_state,
                                          cells$life_duration,
                                          mut_param)
        cells
}

#' Constructor for bulk cell state counts of a single type
#' @export
make_cells_state <- function(time, count, active, mut_counts_list=NULL) {
        out = list(time = time,
                   count = count,
                   active = active,
                   mut_counts_list = mut_counts_list)
        class(out) = "cells_state"
        check_valid_cells_state(out)
        out
}
#' Check validity for cells state
#' @export
check_valid_cells_state <- function(cells_state) {
        assertthat::assert_that(length(cells_state$time) == 1)
        assertthat::assert_that(length(cells_state$count) == 1)
        assertthat::assert_that(length(cells_state$active) == 1)
        if (!is.null(cells_state$mut_counts_list)) {
                assertthat::assert_that(all(sapply(cells_state$mut_counts_list, sum) == cells_state$count))
        }
        return(TRUE)
}

#' Extend the counts for a single cell type until differentiation or target time
#' Cases:
#' 1. cells can self-renew and differentiate for one generation before target time
#' 2. cells will self-renew until target time
#' 3. cells cannot divide but can mutate, make it inactive
#' For mutation:
#' Draw one sample for count of mutant alleles of the cell type at the end time
#' @param differentiate if the cell can differentiate before target time, whether the result for before or after differentiation should be returned.
#' @return a list of cell states, each element is a cell type
#' @export
count_single_type <- function(cell_type, cells_state, type_graph,
                              mut_param=NULL,
                              differentiate = T) {
        mut_counts_list = cells_state$mut_counts_list
        if (!is.null(mut_counts_list)) {
                assertthat::assert_that(!is.null(mut_param))
        }
        diff_time = type_graph$diff_time[[cell_type]]
        diff_indicator = F
        life_duration = type_graph$lifetime[[cell_type]]
        if (type_graph$target_time <= diff_time) {
                # This situation could be case 2 or case 3, behavior decided by \
                # target time
                # NOTE: do not count the division that happens after target time
                diff_indicator = F
                num_gen = floor((type_graph$target_time - cells_state$time)/life_duration)
                if (num_gen == 0){
                        # case 3
                        time_remain = type_graph$target_time - cells_state$time
                        if (!is.null(mut_counts_list)) {
                                m_list = sample_all_mut_counts_gen(mut_counts_list,
                                                                   mut_param = mut_param,
                                                                   life_duration = time_remain,
                                                                   num_gen = 1,
                                                                   double_counts = F)
                        } else {
                                m_list = NULL
                        }
                        cells_state$time = type_graph$target_time
                        cells_state$mut_counts_list = m_list
                        cells_state$active = F
                        out_list = list()
                        out_list[[cell_type]] = cells_state
                        return(out_list)
                }
        } else {
                # In this case the last division before differentiation may be
                # after target time, which means still no differentiation
                # expected number of generations regardless of target time
                num_gen = ceiling((diff_time - cells_state$time)/life_duration)
                if (cells_state$time + num_gen * life_duration > type_graph$target_time) {
                        num_gen = num_gen - 1
                        if (num_gen == 0){
                                # case 3
                                time_remain = type_graph$target_time - cells_state$time
                                if (!is.null(mut_counts_list)) {
                                        m_list = sample_all_mut_counts_gen(mut_counts_list,
                                                                           mut_param = mut_param,
                                                                           life_duration = time_remain,
                                                                           num_gen = 1,
                                                                           double_counts = F)
                                } else {
                                        m_list = NULL
                                }
                                cells_state$time = type_graph$target_time
                                cells_state$mut_counts_list = m_list
                                cells_state$active = F
                                out_list = list()
                                out_list[[cell_type]] = cells_state
                                return(out_list)
                        }
                        diff_indicator = F
                } else {
                        diff_indicator = T
                }
        }
        if (diff_indicator & differentiate) {
                mode_counts = round_frac(cells_state$count * 2^(num_gen-1),
                                         type_graph$diff_mode_probs[[cell_type]])
                daughter_counts = (mode_counts %*% mode_offspring)[1, ]
                names(daughter_counts) =
                        c(cell_type, type_graph$merge[cell_type, ])
                ## NOTE: no mode 1 for now
                assertthat::assert_that(daughter_counts[1] == 0)
                daughter_counts = daughter_counts[-1]
                assertthat::assert_that(length(daughter_counts) == 2)
                if (!is.null(mut_counts_list)) {
                        m_list = sample_all_mut_counts_gen(mut_counts_list,
                                                           mut_param = mut_param,
                                                           life_duration = life_duration,
                                                           num_gen = num_gen)
                        assertthat::assert_that(mut_counts_total(m_list) == sum(daughter_counts))
                        m_list_split = distribute_all_mut_counts(m_list,
                                                                 n1 = daughter_counts[1],
                                                                 n2 = daughter_counts[2])
                }
                out = lapply(1:2, function(i) {
                        if (!is.null(mut_counts_list)) {
                                m_list1 = lapply(m_list_split, "[[", i)
                        } else {
                                m_list1 = NULL
                        }
                        make_cells_state(time = cells_state$time + num_gen * life_duration,
                                         count =  daughter_counts[i],
                                         active = T,
                                         mut_counts_list = m_list1)
                })
                names(out) = names(daughter_counts)
                return(out[!sapply(out, "[[", "count") == 0])
        } else {
                if ((!differentiate) & diff_indicator) {
                        num_gen = num_gen - 1
                }
                if (!is.null(mut_counts_list)) {
                        m_list = sample_all_mut_counts_gen(mut_counts_list,
                                                           mut_param = mut_param,
                                                           life_duration = life_duration,
                                                           num_gen = num_gen)
                } else {
                        m_list = NULL
                }
                out = list()
                out[[cell_type]] = make_cells_state(time = cells_state$time + num_gen * life_duration,
                                                    count = cells_state$count * 2 ^ num_gen,
                                                    active = T,
                                                    mut_counts_list = m_list)
                return(out)
        }
}
#' Advance a list of cell states
#' @export
advance_cells_state <- function(cells_state_list, type_graph, mut_param = NULL) {
        do.call(c, lapply(names(cells_state_list), function(cell_type) {
                cells_state = cells_state_list[[cell_type]]
                if (!cells_state$active) {
                        # pass on inactive cells state
                        out_list = list()
                        out_list[[cell_type]] = cells_state
                        return(out_list)
                } else {
                        out_list = count_single_type(cell_type = cell_type,
                                                     cells_state = cells_state,
                                                     type_graph = type_graph,
                                                     mut_param = mut_param)

                        return(out_list)
                }
        }))
}
#' Get the results at target time
#' @export
get_target_states <- function(init_state_list, type_graph, mut_param=NULL) {
        # This function needs further testing to ensure consistency with non-sampling version
        # NOTE: there may be bugs when target time is close to latest differentiation
        cur_states_list = init_state_list
        while (any(sapply(cur_states_list, "[[", "active"))) {
                new_state_list = advance_cells_state(cur_states_list,
                                                     type_graph = type_graph,
                                                     mut_param = mut_param)
                cur_states_list = new_state_list
        }
        return(cur_states_list)
}
#' Convert single cell results to a list of cell states
#' @export
cells_to_state_list <- function(cells, convert_mut = F) {
        types_queue = cells$type_state
        state_list = list()
        for (cell_type in types_queue) {
                type_time = cells$birth_time[cells$type_state == cell_type]
                assertthat::assert_that(all(type_time == type_time[1]))
                if (convert_mut) {
                        # TODO: test consistency between bulk and sc
                        m_list = lapply(1:ncol(cells$mut_state),
                                        function(j) {
                                                table(cells$mut_state[, j])
                                        })
                } else {
                        m_list = NULL
                }
                state_list[[cell_type]] =
                        make_cells_state(time = type_time[1],
                                         count = sum(cells$type_state == cell_type),
                                         active = T,
                                         mut_counts_list = m_list)
        }
        return(state_list)
}
#' Get counts from state list
#' @export
get_counts_from_states_list <- function(states_list, types=NULL) {
        counts = sapply(states_list, "[[", "count")
        names(counts) = names(states_list)
        if (is.null(types)) {
                return(counts)
        } else {
                assertthat::assert_that(all(names(counts) %in% types))
                type_counts = numeric(length(types))
                names(type_counts) = types
                for (ctype in names(counts)) {
                        type_counts[ctype] = counts[ctype]
                }
                return(type_counts)
        }
}
#' Make initial states list
#' @export
initial_states_list_with_mut <- function(type_graph, mut_param, count = 1) {
        states_list = list()
        states_list[[type_graph$root_id]] =  make_cells_state(time = 0,
                                                              count = count,
                                                              active = T,
                                                              mut_counts_list = init_mut_counts_list(mut_param = mut_param,
                                                                                                     num_cells = count))
        states_list
}












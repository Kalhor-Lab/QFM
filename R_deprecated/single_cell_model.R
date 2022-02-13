# Single cell sampling model
# Under the case of sample size equal to target counts,
# should be consistent with bulk counting functions

#' Simulate a new generation of cells given last generation
#' A generation contains active and inactive cells, which have reached target time
#' 1. For the active cells:
#' 1a. decide daughter cell types and double time for dividing cells
#' 1b. sampling: 1b1. decide daughter counts at target time
#'               1b2. decide daughter sampled cell counts based on parent sampled counts
#'               1b3. remove unsampled daughter cells
#' 2. mutate daughter barcodes
#' 3. split sampled daughter cells into active vs. inactive
#' 4. combine current inactive with previous inactives
#' @return a new generation with inactive cells and daughter cells for the active cells
simulate_generation <- function(prev_gen, type_graph, mut_param) {
        # prev_gen = sim_history[[7]]
        daughter_cells = generate_daughter_cells(prev_gen$active,
                                                 type_graph = type_graph)
        assertthat::assert_that(all(prev_gen$active$id %in% daughter_cells$parent))
        daughter_cells = mutate_cells(daughter_cells, mut_param = mut_param)
        cells_sp = split_cells(daughter_cells, daughter_cells$active)
        names(cells_sp) = c('active', 'inactive')
        if (num_cells(prev_gen$inactive) > 0) {
                cells_sp$inactive = combine_cells(list(prev_gen$inactive,
                                                       cells_sp[['inactive']]))
        }
        return(cells_sp)
}
#' Distribute sample size by type
distribute_sample_size <- function(parent_target_counts,
                                   daughter_target_counts,
                                   parent_sample_size) {
        parent_sample_size = rbind(c("-1"=1000,
                                     "-2"=1000,
                                     "-3"=1000,
                                     "-4"=1000,
                                     "-5"=1000,
                                     "-6"=1000))
        # parent_target_counts = cbind(c("-1"=10000,
        #                                "-2"=10000,
        #                                "-3"=10000,
        #                                "-4"=10000,
        #                                "-5"=20000,
        #                                "-6"=20000))
        daughter_target_counts = rbind(c("-1"=0,
                                         "-2"=5000,
                                         "-3"=7000,
                                         "-4"=7000,
                                         "-5"=5000,
                                         "-6"=0),
                                       c("-1"=10000,
                                         "-2"=5000,
                                         "-3"=0,
                                         "-4"=3000,
                                         "-5"=15000,
                                         "-6"=15000))
        # assertthat::assert_that()
        parent_target_counts = do.call(rbind, lapply(1:1, function(i) colSums(daughter_target_counts[(2*i-1):(2*i), ])))
        parent_target_expected = sum(parent_target_counts)/length(parent_target_counts)

}



#' Generate daughter cells object of the next generation
#' a. Decide types and double time
#' b. Get target counts and decide Sample size for daughters
generate_daughter_cells <- function(cells, type_graph) {
        # cells = prev_gen$active
        # process each type separately
        types_queue = sort(unique(cells$type_state))
        type_lifetime = unlist(type_graph$lifetime)
        cells_by_type = lapply(types_queue, function(cell_type) {
                # cell_type = types_queue[1]
                # a. decide daughter cell types and life time
                cells_subset = subset_cells(cells, cells$type_state == cell_type)
                death_time = cells_subset$birth_time + cells_subset$life_duration
                # death time should all be the same for the same cell type
                assertthat::assert_that(is_zero_range(death_time))
                # numeric problem, temporary fix, or avoid using fractions
                diff_ind = cells_subset$birth_time >= type_graph$diff_time[[cell_type]]
                # further narrow to differentiated cells and decide differentiated mode
                if (any(diff_ind > 0)) {
                        diff_modes = numeric(sum(cells$type_state == cell_type))
                        diff_modes_counts = round_frac(sum(diff_ind),
                                                       type_graph$diff_mode_prob[[cell_type]])
                        # Important: randomly shuffle modes among parent cells
                        diff_modes[diff_ind] = rep(2:6, diff_modes_counts)[sample(sum(diff_modes_counts))]
                        # old stochastic mode assignment
                        # diff_modes[diff_ind] = sample(2:6, size = sum(diff_ind),
                        #                                      prob = type_graph$diff_mode_prob[[cell_type]],
                        #                                      replace = T)
                        diff_modes[!diff_ind] = 1
                        daughter_types = c(t(type_graph$diff_outcome[[cell_type]][diff_modes, ]))
                } else {
                        daughter_types = rep(cell_type, 2 * sum(cells$type_state == cell_type))
                }
                # daughter life durations
                target_time_remain = type_graph$target_time - death_time
                daughter_life_durations =
                        pmin(target_time_remain,
                             type_lifetime[match(daughter_types, names(type_lifetime))])

                # b. decide daughter activity and target counts
                daughter_types_queue = sort(unique(daughter_types))
                daughter_type_activity = sapply(daughter_types_queue, function(dp) {
                        # this is for testing activity only
                        state_list = count_single_type(cell_type = dp,
                                                       cells_state = make_cells_state(time = death_time[1],
                                                                                     count = sum(daughter_types == dp),
                                                                                     active = T),
                                                       type_graph = type_graph)
                        if (length(state_list) == 1) {
                                if (!state_list[[1]]$active) {
                                        return(FALSE)
                                } else {
                                        return(TRUE)
                                }
                        } else {
                                return(TRUE)
                        }
                })
                daughter_activity = daughter_type_activity[match(daughter_types, daughter_types_queue)]
                # TODO: correctly handle non-integers
                daughter_target_counts = numeric(length(daughter_types))
                for (dp in daughter_types_queue) {
                        type_count = sum(daughter_types == dp)
                        state_list = list()
                        state_list[[dp]] = make_cells_state(time = death_time[1],
                                                            count = type_count,
                                                            active = T)
                        target_state = get_target_states(state_list,
                                                         type_graph = type_graph)
                        # averaged over possible cell fates
                        target_counts = round_frac(total = sum(sapply(target_state, "[[", "count")),
                                                prob = rep(1, type_count)/type_count)
                        daughter_target_counts[daughter_types == dp] = target_counts
                }
                # redistributes parent sample size when daughter has multiple types
                # has to do with how cell fate is decided in relation to symmetric vs assymetric division!
                parent_target_counts =
                        sapply(1:num_cells(cells_subset), function(i) sum(daughter_target_counts[(2*i-1):(2*i)]))
                parent_target_expected = sum(parent_target_counts)/length(parent_target_counts)
                #TODO: better handle non-integers, uses stochastic rounding for now
                cells_subset$sample_size = scale_counts(cells_subset$sample_size,
                                                        parent_target_counts/parent_target_expected)
                assertthat::assert_that(all(cells_subset$sample_size >= 1))
                daughter_sample_size_1 = sapply(1:num_cells(cells_subset), function(i) {
                        assertthat::assert_that(sum(daughter_target_counts[(2*i-1):(2*i)]) >=
                                            cells_subset$sample_size[i],
                                    msg = "sample size larger than total cells.")
                        rhyper(nn=1,
                               m=daughter_target_counts[2*i-1],
                               n=daughter_target_counts[2*i],
                               k = cells_subset$sample_size[i])
                })
                daughter_sample_size_2 = cells_subset$sample_size - daughter_sample_size_1
                daughter_sample_size = c(rbind(daughter_sample_size_1, daughter_sample_size_2))
                daughter_sample_ind = (daughter_sample_size > 0)
                return(make_cells(id = generate_cell_ids(sum(daughter_sample_ind)),
                                  parent = rep(cells_subset$id, each = 2)[daughter_sample_ind],
                                  birth_time = rep(death_time, each = 2)[daughter_sample_ind],
                                  life_duration = daughter_life_durations[daughter_sample_ind],
                                  type_state = daughter_types[daughter_sample_ind],
                                  sample_size = daughter_sample_size[daughter_sample_ind],
                                  active = daughter_activity[daughter_sample_ind],
                                  mut_state = cells_subset$mut_state[rep(1:nrow(cells_subset$mut_state),
                                                                         each = 2), ][daughter_sample_ind, ]))
        })
        daughters = combine_cells(cells_by_type)
        return(daughters)
}
#' Initialize single cell simulation
sim_init <- function(type_graph, mut_param, sample_size) {
        global_id_counter <<- 0
        return(list(active = make_cells(id = generate_cell_ids(1),
                                        parent=NA,
                                        birth_time = 0,
                                        life_duration=type_graph$lifetime[[type_graph$root_id]],
                                        type_state=type_graph$root_id,
                                        # sample_size = sample_size,
                                        sample_size=cbind(c("-1"=1000,
                                                            "-2"=1000,
                                                            "-3"=1000,
                                                            "-4"=1000,
                                                            "-5"=1000)),
                                        active = T,
                                        mut_state = mutate_barcodes(initialize_mut_state(mut_param, num_cells = 1),
                                                                    type_graph$lifetime[[type_graph$root_id]],
                                                                    mut_param = mut_param))))
}
#' Run simulation for once
#' @export
run_sim <- function(type_graph, mut_param, sample_size, max_gen = 64) {
        sim_history = list()
        sim_history[[1]] = sim_init(type_graph = type_graph,
                                    mut_param = mut_param,
                                    sample_size = sample_size)
        cur_gen = sim_history[[1]]
        # cur_gen = sim_history[[length(sim_history)]]
        for (i in 1:max_gen) {
                message(i)
                target_state = get_target_states(cells_to_state_list(cur_gen$active),
                                                 type_graph = type_graph)
                assert_that(sum(cur_gen$active$sample_size) <= sum(sapply(target_state, "[[", "count")),
                            msg = "sample size larger than total cells")
                next_gen = simulate_generation(cur_gen,
                                               type_graph = type_graph,
                                               mut_param = mut_param)
                sim_history[[i+1]] = next_gen
                if (num_cells(next_gen$active) == 0) {
                        message(paste0('generation: ', i))
                        break()
                }
                cur_gen = next_gen
        }
        assertthat::assert_that(get_sample_size(sim_history) == sample_size)
        assertthat::assert_that(get_target_time(sim_history) == type_graph$target_time)
        return(sim_history)
}

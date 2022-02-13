#' Constructs a cells object which contains
#' 1. Cell ids
#' 2. Parent ids
#' 3. Life duration
#' 4. Cell type states
#' 5. Mutant allele states, a list of mutating elements
#' 6. Number of sampled progenies
#' 7. Indicator of cell activity
#' @export
make_cells <- function(id, parent, birth_time, life_duration,
                       type_state, active, sample_size, mut_state=NULL) {
        out = list(id = id,
                   parent = parent,
                   birth_time = birth_time,
                   life_duration = life_duration,
                   type_state = type_state,
                   mut_state = mut_state,
                   sample_size = sample_size,
                   active = active)
        class(out) = "cells"
        check_valid_cells(out)
        out
}

#' Generate a new batch of cell ids
generate_cell_ids <- function(n) {
        new_ids = (global_id_counter):(global_id_counter+n-1)
        global_id_counter <<- global_id_counter + n
        new_ids
}

#' Check if cells object is NA
#' @export
check_na_cells <- function(cells) {
        if (length(cells) == 1) {
                assertthat::assert_that(is.na(cells))
                return(TRUE)
        } else {
                return(FALSE)
        }
}

#' Count number of cells in object
#' @export
num_cells <- function(cells) {
        check_valid_cells(cells)
        if (check_na_cells(cells)){
                return(0)
        } else {
                return(length(cells$id))
        }
}

#' Check validity of cells object
#' @export
check_valid_cells <- function(cells) {
        if (!check_na_cells(cells)) {
                assertthat::assert_that(length(cells$id) == length(cells$parent))
                assertthat::assert_that(length(cells$id) == length(cells$birth_time))
                assertthat::assert_that(length(cells$id) == length(cells$life_duration))
                assertthat::assert_that(length(cells$id) == length(cells$type_state))
                assertthat::assert_that(length(cells$id) == nrow(cells$sample_size))
                assertthat::assert_that(length(cells$id) == length(cells$active))
                if (!is.null(cells$mut_state)) {
                        assertthat::assert_that(length(cells$id) == nrow(cells$mut_state))
                }
        }
        return(TRUE)
}

#' Splits a cells object into two objects based on a indicator vector
#' @export
split_cells <- function(cells, ind) {
        check_valid_cells(cells)
        assertthat::assert_that(num_cells(cells) == length(ind))
        if (any(ind)) {
                batch0 = make_cells(id = cells$id[ind],
                                    parent = cells$parent[ind],
                                    birth_time = cells$birth_time[ind],
                                    life_duration = cells$life_duration[ind],
                                    type_state = cells$type_state[ind],
                                    mut_state = cells$mut_state[ind, ,drop=F],
                                    sample_size = cells$sample_size[ind],
                                    active = cells$active[ind])
        } else {
                batch0 = NA
        }
        if (all(ind)) {
                batch1 = NA
        } else {
                batch1 = make_cells(id = cells$id[!ind],
                                    parent = cells$parent[!ind],
                                    birth_time = cells$birth_time[!ind],
                                    life_duration = cells$life_duration[!ind],
                                    type_state = cells$type_state[!ind],
                                    mut_state = cells$mut_state[!ind, ,drop=F],
                                    sample_size=cells$sample_size[!ind],
                                    active = cells$active[!ind])
        }
        list(batch0, batch1)
}

#' subset cells based on indicator vector
#' @export
subset_cells <- function(cells, ind) {
        check_valid_cells(cells)
        assertthat::assert_that(num_cells(cells) == length(ind))
        make_cells(id = cells$id[ind],
                   parent = cells$parent[ind],
                   birth_time = cells$birth_time[ind],
                   life_duration = cells$life_duration[ind],
                   type_state = cells$type_state[ind],
                   mut_state = cells$mut_state[ind, ,drop=F],
                   sample_size=cells$sample_size[ind],
                   active = cells$active[ind])
}

#' Combine a list of cells
#' @export
combine_cells <- function(cells_list) {
        cells_list = lapply(cells_list, function(x) {
                if (check_na_cells(x)){
                        return(make_cells(NULL, NULL, NULL, NULL,
                                          NULL, NULL, NULL, NULL))
                } else {
                        return(x)
                }
        })
        make_cells(id = do.call(c, lapply(cells_list, "[[", "id")),
                   parent = do.call(c, lapply(cells_list, "[[", "parent")),
                   birth_time = do.call(c, lapply(cells_list, "[[", "birth_time")),
                   life_duration = do.call(c, lapply(cells_list, "[[", "life_duration")),
                   type_state = do.call(c, lapply(cells_list, "[[", "type_state")),
                   mut_state = do.call(rbind, lapply(cells_list, "[[", "mut_state")),
                   sample_size = do.call(c, lapply(cells_list, "[[", "sample_size")),
                   active = do.call(c, lapply(cells_list, "[[", "active")))
}


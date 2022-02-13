# this is a rebuild of the cosar count function that I lost
state_mapper <- c("N_N" = "N",
                  "A_A" = "A",
                  "B_B" = "B",
                  "A_B" = "R",
                  "B_A" = "R",
                  "A_N" = "A",
                  "N_A" = "A",
                  "B_N" = "B",
                  "N_B" = "B",
                  "R_R" = "R",
                  "R_A" = "R",
                  "A_R" = "R",
                  "R_B" = "R",
                  "B_R" = "R",
                  "R_N" = "R",
                  "N_R" = "R")

simulate_cosar_counts <- function(size, split_ratio, collect_size) {
        size_split = round_frac(size, split_ratio)
        size_split_sampled = c(colSums(rmultinom(n = 1, size = collect_size[1], prob = rep(1, size_split[1])/size_split[1]) > 0),
                               colSums(rmultinom(n = 1, size = collect_size[2], prob = rep(1, size_split[2])/size_split[2]) > 0))
        field_vec = c(c(rep("A", size_split_sampled[1]), rep("N", size_split[1] - size_split_sampled[1])),
                      c(rep("B", size_split_sampled[2]), rep("N", size_split[2] - size_split_sampled[2])))
        cosar_count_list = list()
        while (any(field_vec %in% c("A", "B"))) {
                field_vec = sample(field_vec)
                if (length(field_vec) %% 2 == 1) {
                        field_vec = c(field_vec, "N")
                }
                field_mat = matrix(field_vec, ncol = 2)
                field_coal = map_chr(1:nrow(field_mat), function(i) paste0(field_mat[i, ], collapse = "_"))
                count_disjoint = sum(field_coal %in% c("A_B", "B_A"))
                count_nondisjoint_A = sum(field_coal %in% c("A_R", "R_A"))
                count_nondisjoint_B = sum(field_coal %in% c("R_B", "B_R"))
                cosar_count_list = append(cosar_count_list,
                                          list(
                                                  c(count_disjoint, count_nondisjoint_A, count_nondisjoint_B))
                )
                field_vec = state_mapper[field_coal]
        }
        cosar_count_mat = do.call(rbind, cosar_count_list)
        cosar_count_mat
}
compute_mean_cosar_error <- function(num_sim = 100, ...) {
        map_dbl(1:num_sim, function(i) {
                count_mat = simulate_cosar_counts(...)
                disjoint_count_vec = count_mat[, 1]
                sum((0:-(length(disjoint_count_vec)-1)) *
                            (disjoint_count_vec/sum(disjoint_count_vec)))
        })
}








# adjustment curve
# this idea has failed ####
# allow variable division rate
# going to base this off the old code, by be relatively independent from the simulation itself
# start count is always 1
# always doubles but with custom interdivision time, no double_time parameter
make_variable_gens <- function(div_time, active, diff_mode_probs = NULL, start_time =0, start_count = 1, start_mode_counts = NA) {
        # div_time is a vector of interdivision times
        if (active) {
                assertthat::assert_that(!is.null(diff_mode_probs))
                n_gen = length(div_time) - 1
        } else {
                n_gen = length(div_time)
        }
        list(cell_type = NA,
             start_time = start_time,
             start_count = start_count,
             start_mode_counts = start_mode_counts,
             cell_counts = start_count * 2^(0:n_gen),
             num_gen = n_gen,
             div_time = div_time[1:n_gen],
             gen_time = start_time + c(0, cumsum(div_time[1:n_gen])),
             end_time = start_time + sum(div_time[1:n_gen]),
             end_count = start_count * 2^n_gen,
             end_mode_counts = get_mode_counts_simple(start_count * 2^n_gen, diff_mode_probs),
             active = active)
}
get_mode_counts_simple <- function(count, diff_mode_probs) {
        round_frac(count,
                   diff_mode_probs)
}
distribute_mode_counts_simple <- function(mode_counts) {
        # without naming the types
        # no mode 3 -5
        assertthat::assert_that(all(mode_counts[3:5] == 0))
        daughter_mode_counts = mode_counts * mode_offspring
        daughter_mode_counts
}
get_mode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
}

# input diff trajectory
# output COSAR distribution
# div_time_all = c(rep(0.6, 7), rep(0.35, 23))
# sample_size1 = 500
# sample_size2 = 500
# diff_mode_prob = c(0.5, 0.5, 0, 0, 0)
# sample_curve <- function(diff_gen_seq) {
#         out = map(diff_gen_seq, function(diff_gen) {
#                 div_time_undiff = div_time_all[1:diff_gen]
#                 div_time_diff1 = div_time_all[(diff_gen+1):length(div_time_all)]
#                 div_time_diff2 = div_time_all[(diff_gen+1):length(div_time_all)]
#                 gens_undiff = make_variable_gens(div_time_undiff, active = T, diff_mode_prob = diff_mode_prob)
#                 # differentiate
#                 diff_mode_counts = distribute_mode_counts_simple(gens_undiff$end_mode_counts)
#                 gens_diff1 = make_variable_gens(div_time_diff1,
#                                                 active = F,
#                                                 start_time = sum(div_time_undiff),
#                                                 start_count = sum(diff_mode_counts[, 2]),
#                                                 start_mode_counts = diff_mode_counts[, 2])
#                 gens_diff2 = make_variable_gens(div_time_diff2,
#                                                 active = F,
#                                                 start_time = sum(div_time_undiff),
#                                                 start_count = sum(diff_mode_counts[, 3]),
#                                                 start_mode_counts = diff_mode_counts[, 3])
#
#                 gens_diff1 = generate_sample_size(gens_diff1, sample_size = pmin(sample_size1, gens_diff1$end_count))
#                 gens_diff2 = generate_sample_size(gens_diff2, sample_size = pmin(sample_size2, gens_diff2$end_count))
#                 gens_undiff = merge_sample_size(gens_diff1, gens_diff2, gens_undiff)
#
#                 # simulate COSARs without barcodes
#                 # algorithm: keep track of current types in sample R-A,B
#                 apparent_types = c(rep("A", gens_undiff$end_mode_sample_size[1]),
#                                    rep("B", gens_undiff$end_mode_sample_size[2]))
#                 cosar_time = numeric()
#                 for (i in gens_undiff$num_gen:1) {
#                         n_double = gens_undiff$sample_size[i+1] - gens_undiff$sample_size[i]
#                         n_single = gens_undiff$sample_size[i] - n_double
#                         if (n_single > 0) {
#                                 single_group = paste0("S-", 1:n_single)
#                         } else {
#                                 single_group = character()
#                         }
#                         if (n_double > 0) {
#                                 double_group = rep(paste0("D-", 1:n_double), each = 2)
#                         } else {
#                                 double_group = character()
#                         }
#                         group_vec = c(single_group, double_group)
#                         apparent_types_group = split(sample(apparent_types), group_vec)
#                         apparent_types_new = map_chr(apparent_types_group, function(x) {
#                                 if (all(x == "A")) {
#                                         return("A")
#                                 }
#                                 if (all(x == "B")) {
#                                         return("B")
#                                 }
#                                 if ("R" %in% x) {
#                                         return("R")
#                                 }
#                                 if ("A" %in% x & "B" %in% x) {
#                                         return("C")
#                                 }
#                         })
#                         cosar_time = c(cosar_time, rep(gens_undiff$gen_time[i], sum(apparent_types_new == "C")))
#                         apparent_types_new[apparent_types_new == "C"] = "R"
#                         apparent_types = apparent_types_new
#                 }
#                 return(list(sample_size = gens_undiff$sample_size[gens_undiff$num_gen],
#                             mf_cosar_time = get_mode(cosar_time)))
#         })
#         return(out)
# }
#
# curve_tb = bind_rows(map(1:100, function(i) {
#         message(i)
#         diff_gen_seq = 5:30
#         curve_list = sample_curve(diff_gen_seq)
#         tibble(diff_time = cumsum(div_time_all)[diff_gen_seq-2],
#                mf_cosar_time = map_dbl(curve_list, "mf_cosar_time"),
#                log2_field_size = diff_gen_seq - 2,
#                log2_sample_size = log2(map_dbl(curve_list, "sample_size")))
# }))
# ggline(curve_tb, x = "diff_time", y = "mf_cosar_time", add = "mean_sd", numeric.x.axis = T) + geom_abline(col = "red")
# ggline(curve_tb, x = "log2_field_size", y = "log2_sample_size", add = "mean_sd", numeric.x.axis = T) + geom_abline(col = "red")
#
#








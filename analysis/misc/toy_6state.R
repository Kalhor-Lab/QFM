# # Graph Definition:
root_id = "5"
node_id = c("1", "2", "3", "4", "5")
tip_id = c("-1", "-2", "-3", "-4", "-5", "-6")

# # transition fraction
# trans_sel1 = list("5" = list(list(fraction = c(0.5, 0),
#                                   time = 3.5),
#                              list(fraction = c(0.5, 0),
#                                   time = 5.0),
#                              list(fraction = c(0.02, 0),
#                                   time = 6.5)))
# trans_sel0 = list("3" = list(list(fraction = c(0.2, 0),
#                                   time = 5.0)),
#                   "-4" = list(list(fraction = c(0.3, 0),
#                                    time = 9.5)))

# topology
merge0 = do.call(rbind, list("1" = c(-1, -2),
                             "2" = c(-3, -4),
                             "3" = c(-5, -6),
                             "4" = c(1, 2),
                             "5" = c(3, 4)))
merge1 = do.call(rbind, list("1" = c(-1, -2),
                             "2" = c(1, -3),
                             "3" = c(2, -4),
                             "4" = c(3, -5),
                             "5" = c(4, -6)))
# division type probs for the nodes c(11, 22, 12, 01, 02)
type_mode_probs0 = list("1" = c(0.9, 0.1, 0, 0, 0),
                        "2" = c(0.3125, 0.6875, 0, 0, 0),
                        "3" = c(0.5, 0.5, 0, 0, 0),
                        "4" = c(0.6, 0.4, 0, 0, 0),
                        "5" = c(0.725, 0.375, 0, 0, 0))
type_mode_probs1 = list("1" = c(0.8, 0.2, 0, 0, 0),
                        "2" = c(0.8, 0.2, 0, 0, 0),
                        "3" = c(0.8, 0.2, 0, 0, 0),
                        "4" = c(0.8, 0.2, 0, 0, 0),
                        "5" = c(0.8, 0.2, 0, 0, 0))
type_time0 = list("1" = 5.0,
                  "2" = 6.5,
                  "3" = 6.0,
                  "4" = 4.0,
                  "5" = 3.0)
tip_time = list("-1" = Inf,
                "-2" = Inf,
                "-3" = Inf,
                "-4" = Inf,
                "-5" = Inf,
                "-6" = Inf)
type_time0 = c(type_time0, tip_time)
type_time2 = list("1" = 13.0,
                  "2" = 12.0,
                  "3" = 10.5,
                  "4" = 9.0,
                  "5" = 8.0)
type_time2 = c(type_time2, tip_time)
# type_diff_time2 = list("1" = 6.5,
#                        "2" = 4.5,
#                        "3" = 3.5,
#                        "4" = 2.5,
#                        "-1" = Inf,
#                        "-2" = Inf,
#                        "-3" = Inf,
#                        "-4" = Inf,
#                        "-5" = Inf)
# type_diff_mode_probs1 = list("1" = c(0.33, 0.67, 0, 0, 0),
#                              "2" = c(0.4, 0.6, 0, 0, 0),
#                              "3" = c(0.3, 0.7, 0, 0, 0),
#                              "4" = c(0.5, 0.5, 0, 0, 0))
type_doubletime = list("1" = 0.6,
                       "2" = 0.6,
                       "3" = 0.6,
                       "4" = 0.6,
                       "5" = 0.6,
                       "-1" = 0.35,
                       "-2" = 0.35,
                       "-3" = 0.35,
                       "-4" = 0.35,
                       "-5" = 0.35,
                       "-6" = 0.35)
type_graph0 = make_type_graph(root_id = root_id,
                              node_id = node_id,
                              tip_id = tip_id,
                              merge_matrix = merge0,
                              differntiation_time = type_time0,
                              differntiation_mode_probs = type_mode_probs0,
                              double_time = type_doubletime,
                              founder_size = 1,
                              target_time = 10)
type_time1 = rlang::duplicate(type_time0)
type_time1[["2"]] = 6.0
type_time1 = c(type_time1, tip_time)
type_graph1 = make_type_graph(root_id = root_id,
                              node_id = node_id,
                              tip_id = tip_id,
                              merge_matrix = merge0,
                              differntiation_time = type_time1,
                              differntiation_mode_probs = type_mode_probs0,
                              double_time = type_doubletime,
                              founder_size = 1,
                              target_time = 10)


type_graph2 = make_type_graph(root_id = root_id,
                              node_id = node_id,
                              tip_id = tip_id,
                              merge_matrix = merge1,
                              differntiation_time = type_time2,
                              differntiation_mode_probs = type_mode_probs1,
                              double_time = type_doubletime,
                              founder_size = 1,
                              target_time = 15.0)
plot_type_graph_clean(type_graph1, node_col_mapper = type_col_mapper)

type_col_mapper <- function(types) {
        types = as.numeric(types)
        col_p0 = brewer.pal(11, 'RdYlBu')
        col_p1 = brewer.pal(9, 'RdYlGn')
        col_id = types
        col_id[types == 5] = col_p0[6]
        col_id[types == 4] = col_p0[5]
        col_id[types == -1] = col_p0[9]
        col_id[types == -2] = col_p0[10]
        col_id[types == -3] = col_p0[2]
        col_id[types == -4] = col_p0[1]
        col_id[types == -5] = col_p1[8]
        col_id[types == -6] = col_p1[9]
        col_id[types == 3] = col_p1[7]
        col_id[types == 2] = col_p0[3]
        col_id[types == 1] = col_p0[8]
        names(col_id) = types
        col_id
}
# plot_type_graph(type_graph2_a, type_col_mapper2)
plot_type_graph(type_graph0, type_col_mapper)

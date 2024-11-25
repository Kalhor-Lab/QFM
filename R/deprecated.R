# deprecated, replaced by new version that that involves multifurcation
# get_edge_diff <- function(x, data_obj, trans_time) {
#         diff_time = trans_time[x]
#         gr_tip_list = data_obj$gr_tip_list
#         gr_dd = data_obj$gr_dd
#         assertthat::assert_that(!is.null(data_obj$tr_edges_state_tb))
#         edge_tb = data_obj$tr_edges_state_tb
#
#         state_include = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[x]]))]
#         state_include_d1 = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[gr_dd[[x]][1]]]))]
#         state_include_d2 = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[gr_dd[[x]][2]]]))]
#
#         edge_tb_diff = mutate(edge_tb, ind = map_lgl(type_path, function(y) any(state_include %in% y)))
#         edge_tb_diff = edge_tb_diff %>% filter(from_time <= diff_time &
#                                                        to_time > diff_time &
#                                                        ind)
#         edge_tb_diff = mutate(edge_tb_diff, d1_ind = map_lgl(type_path, function(y) any(state_include_d1 %in% y))) # old version gr_dd[[x]][1] (old version has been overwritten)
#         edge_tb_diff = mutate(edge_tb_diff, d2_ind = map_lgl(type_path, function(y) any(state_include_d2 %in% y)))
#         assertthat::assert_that(all(!(edge_tb_diff$d1_ind & edge_tb_diff$d2_ind)))
#         edge_tb_diff
# }
# deprecated, replaced by new version that that involves multifurcation
# get_node_size <- function(data_obj, tr_node_assign, trans_time) {
#         # data_obj = update_edge_tb_state(data_obj, tr_node_assign)
#         node_size = map(data_obj$gr$node.label, function(x) {
#                 edge_diff = get_edge_diff(x, data_obj, trans_time)
#                 n0 = length(unique(edge_diff$from))
#                 n1 = sum(edge_diff$d1_ind)
#                 n2 = sum(edge_diff$d2_ind)
#                 out = c(n0, n1, n2)
#                 names(out) = c(x, data_obj$gr_dd[[x]])
#                 out
#         })
#         names(node_size) = data_obj$gr$node.label
#         node_size
# }
# deprecated, defined in new script
# correct_trans_time <- function(gr_trans_time, gr_tr_data) {
#         gr_node_dd = gr_tr_data$gr_dd[gr_tr_data$gr$node.label]
#         for (x in names(gr_node_dd)) {
#                 if(gr_trans_time[gr_node_dd[[x]][1]] < gr_trans_time[x]) {
#                         gr_trans_time[gr_node_dd[[x]][1]] = gr_trans_time[x] + rnorm(1, mean = 0.001, sd = 0.0001)
#                         # message(paste0('corrected: ', gr_node_dd[[x]][1]))
#                 }
#                 if(gr_trans_time[gr_node_dd[[x]][2]] < gr_trans_time[x]) {
#                         gr_trans_time[gr_node_dd[[x]][2]] = gr_trans_time[x] + rnorm(1, mean = 0.001, sd = 0.0001)
#                         # message(paste0('corrected: ', gr_node_dd[[x]][2]))
#                 }
#         }
#         return(gr_trans_time)
# }

# ice_fase <- function(tr,
#                      sc_celltypes,
#                      total_time,
#                      root_time = 0,
#                      theta = 0.0) {
#         # nrounds = 20,
#         # max_depth = 4,
#         # remove_uniform = F,
#         # time_stat_func = mean) {
#         # if (remove_uniform) {
#         #         cell_mat = remove_uniform_id(cell_mat)
#         # }
#         # mut_p = estimate_mut_p(cell_mat, total_time)
#         # mat_im = impute_characters(cell_mat, nrounds = nrounds, max_depth = max_depth)
#         #
#         # tr_upgma = phylotime(mat_im[sample(nrow(mat_im)), ],
#         #                      mut_p = mut_p,
#         #                      total_time)
#         tr = name_nodes(tr)
#         if (is.character(sc_celltypes)) {
#                 assertthat::assert_that(!is.null(names(sc_celltypes)),
#                                         msg = "cell types must be named.")
#                 assertthat::assert_that(all(tr$tip.label %in% names(sc_celltypes)))
#                 cell_type_vec = sc_celltypes
#         } else {
#                 assertthat::assert_that(class(sc_celltypes) == "function")
#                 cell_type_vec = sc_celltypes(tr$tip.label)
#                 names(cell_type_vec) = rownames(tr$tip.label)
#         }
#         gr = name_nodes(reconstruct_graph(tr,
#                                           sc_celltypes = cell_type_vec,
#                                           theta = theta,
#                                           total_time = total_time,
#                                           stat_func = mean))
#         data_obj = make_gr_tr_data(gr, tr, cell_type_vec)
#         tr_node_assign = assign_node_states(data_obj)
#
#         gr_trans_time = est_transition_time(data_obj, tr_node_assign, stat_func = mean)
#         # cases where no node is assigned, assign parent time
#         if (any(!gr$node.label %in% names(gr_trans_time))) {
#                 na_state = gr$node.label[!gr$node.label %in% names(gr_trans_time)]
#                 assertthat::assert_that(all(!na_state %in% tr_node_assign))
#                 gr_trans_time = gr_trans_time[gr$node.label]
#                 names(gr_trans_time) = gr$node.label
#                 gr_trans_time[na_state] = 0.
#         }
#         gr_tip_time = rep(total_time, length(unique(cell_type_vec)))
#         names(gr_tip_time) = unique(cell_type_vec)
#         gr_trans_time = c(gr_trans_time, gr_tip_time)
#         gr_trans_time = correct_trans_time(gr_trans_time, data_obj)
#         gr_trans_time = gr_trans_time + root_time
#
#         data_obj = update_edge_tb_state(data_obj, tr_node_assign)
#         gr_node_sizes = get_node_size(data_obj, tr_node_assign, gr_trans_time)
#
#         out = list(tr = tr,
#                    sc_celltypes = sc_celltypes,
#                    total_time = total_time,
#                    root_time = root_time,
#                    theta = theta,
#                    gr = gr,
#                    gr_trans_time = gr_trans_time,
#                    gr_node_sizes = gr_node_sizes,
#                    tr_node_assign = tr_node_assign,
#                    gr_tr_data = data_obj)
#         class(out) = "ice_fase_results"
#         out
# }
# plot_gr_data <- function(out_data, target_time, gr_col = NULL, gr_lab = NULL) {
#         if (is.null(gr_col)) {
#                 gr_col = gr_color(out_data$gr)
#         }
#         plot_gr(out_data$gr,
#                 out_data$gr_trans_time,
#                 gr_col,
#                 target_time,
#                 node_label = gr_lab)
# }

#' Impute missing characters from cell allele matrix
# impute_characters <- function(mat, nrounds = 100, max_depth = 4) {
#         op = options(na.action='na.pass')
#         on.exit(options(op))
#
#         # separate out part of the matrix that has no variability
#         col_div = apply(mat, 2, function(x) length(unique(x[!is.na(x)])))
#         im_indices = which(col_div > 1)
#         keep_indices = which(col_div == 1)
#
#         mat_keep = rlang::duplicate(mat[, keep_indices, drop =F])
#         if (ncol(mat_keep) > 0) {
#                 for (i in 1:ncol(mat_keep)) {
#                         char_vec = mat_keep[, i]
#                         char_vec_unique = unique(char_vec[!is.na(char_vec)])
#                         assertthat::assert_that(length(char_vec_unique) == 1)
#                         mat_keep[is.na(char_vec), i] = char_vec_unique
#                 }
#         }
#         assertthat::assert_that(all(!is.na(mat_keep)))
#
#         # impute identical elements with the same allele
#         mat_im = as.data.frame(rlang::duplicate(mat[, im_indices, drop =F]))
#         element_id = 1:ncol(mat_im)
#         na_frac = map_dbl(element_id, function(j) {
#                 mean(is.na(mat[, j]))
#         })
#         # imputation order
#         element_id = element_id[order(na_frac)]
#         for (i in element_id) {
#                 # message(i)
#                 design_mat = model.matrix(~., mat_im[, -i])[, -1]
#                 y_fac = factor(as.matrix(mat_im[i])[, 1])
#                 y_vec = as.numeric(y_fac)
#                 suppressWarnings(bst <- xgboost(data = design_mat[!is.na(y_vec), ],
#                                                 label = y_vec[!is.na(y_vec)] - 1,
#                                                 max_depth = max_depth,
#                                                 eta = 1,
#                                                 nthread = 12,
#                                                 nrounds = nrounds,
#                                                 objective = "multi:softmax",
#                                                 eval_metric = "mlogloss",
#                                                 num_class= max(y_vec[!is.na(y_vec)]),
#                                                 verbose = 0))
#                 y_pred = levels(y_fac)[predict(bst, design_mat[is.na(y_vec), ]) + 1]
#                 mat_im[is.na(y_vec), i] = y_pred
#         }
#         mat_im = as.matrix(mat_im)
#         rownames(mat_im) = rownames(mat)
#
#         mat_out = cbind(mat_im, mat_keep)
#         mat_out = mat_out[, colnames(mat)]
#         mat_out
# }


# deprecated
# set_color_palette <- function(res, palette = NULL) {
#         if (is.null(palette)) {
#                 res$col_pal = gr_color(res$gr)
#         } else {
#                 assertthat::assert_that(
#                         length(palette) == length(c(res$gr$node.label, res$gr$tip.label))
#                         )
#                 names(palette) = c(res$gr$node.label, res$gr$tip.label)
#                 res$col_pal = palette
#         }
#         res
# }


# plot_gr <- function(gr,
#                     gr_node_time = NULL,
#                     total_time = NULL,
#                     type_col = NULL,
#                     node_label = NULL) {
#         ig_gr = as.igraph(gr)
#         if (is.null(gr_node_time)) {
#                 gr_node_time = node.depth.edgelength(gr)
#                 names(gr_node_time) = c(gr$tip.label, gr$node.label)
#                 V(ig_gr)$Time = gr_node_time[V(ig_gr)$name]
#         }
#         V(ig_gr)$Time = gr_node_time[V(ig_gr)$name]
#         if (is.null(total_time)) {
#                 total_time = max(gr_node_time)
#         }
#
#         # V(ig_gr)$Type = V(ig_gr)$name
#         if (!is.null(node_label)) {
#                 V(ig_gr)$node_label = node_label[V(ig_gr)$name]
#         } else {
#                 V(ig_gr)$node_label = V(ig_gr)$name
#         }
#         # V(ig_gr)$node_label = rep("", length(V(ig_gr)$name))
#         if (is.null(type_col)) {
#                 type_col = gr_color_v1(gr)
#         }
#         g_out = ggraph(ig_gr, layout = 'dendrogram', height = Time) +
#                 geom_edge_bend(
#                         aes(width = factor(1),
#                             fontface = 'plain'),
#                         arrow = arrow(
#                                 type = "closed",
#                                 length = unit(2, "pt"),
#                                 angle = 45
#                         ),
#                         start_cap = circle(5, 'pt'),
#                         end_cap = circle(5, 'pt')
#                 ) +
#                 geom_node_label(aes(
#                         label = node_label,
#                         fill = name,
#                         size = factor(1)
#                 ),
#                 label.padding = unit(6, 'pt'))
#         g_out = g_out +
#                 scale_fill_manual(values = type_col) +
#                 ylim(c(total_time, 0)) + ylab('Time') +
#                 scale_edge_width_manual(values = 0.5, guide = "none") +
#                 scale_size_manual(values = 4, guide = "none") +
#                 theme(
#                         axis.line.y = element_line(),
#                         axis.ticks.y = element_line(),
#                         axis.text.y = element_text(size = 12),
#                         axis.title = element_text(size = 12),
#                         panel.background = element_rect(fill = NA, color = NA, size = 10),
#                         legend.position = 'none',
#                         axis.title.x = element_blank(),
#                         text = element_text(size = 12)
#                 )
#         g_out
# }

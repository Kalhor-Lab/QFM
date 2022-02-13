# NEW attempt at Fatogram
get_cell_type_comp <- function(data_obj, tt) {
        assertthat::assert_that(!is.null(data_obj$tr_edges_state_tb))
        assertthat::assert_that(!is.null(data_obj$gr_edges_trans))
        edge_tb = data_obj$tr_edges_state_tb
        # edges present at time
        edge_tb_tt = edge_tb %>% filter(from_time <= tt &
                                                to_time > tt)
        state_tt = filter(data_obj$gr_edges_trans, from_time  <= tt & to_time > tt)
        current_state = state_tt$to
        if (length(current_state) == 0) {
                current_state = "Node-1"
        }
        # extending type_path upstream
        edge_tt_states = map_chr(edge_tb_tt$type_path_ext, function(x) {
                out_state = current_state[current_state %in% x]
                assertthat::assert_that(length(out_state) <= 1)
                if (length(out_state) == 0) {
                        return(x[1])
                } else {
                        return(out_state)
                }
        })
        edge_tt_states_count = table(edge_tt_states)
        tibble(time = tt,
               cell_type = names(edge_tt_states_count),
               count = as.numeric(edge_tt_states_count))
}
get_raw_cell_type_comp <- function(data_obj, tt) {
        assertthat::assert_that(!is.null(data_obj$tr_edges_state_tb))
        edge_tb = data_obj$tr_edges_state_tb
        # edges present at time
        edge_tb_tt = edge_tb %>% filter(from_time <= tt &
                                                to_time > tt)
        edge_tt_states_count = table(edge_tb_tt$to_type)
        tibble(time = tt,
               cell_type = names(edge_tt_states_count),
               count = as.numeric(edge_tt_states_count))
}

gr_tr_data = update_gr_trans_time(gr_tr_data, gr_trans_time)
data_obj = gr_tr_data
cell_comp_tb = bind_rows(map(seq(1, 8, length = 200), function(tt) {
        get_cell_type_comp(gr_tr_data, tt)
}))
cell_comp_tb$cell_type = factor(cell_comp_tb$cell_type,
                                levels = lout$name[order(lout$x)])
cell_comp_tb = arrange(cell_comp_tb, time, cell_type)
(ggplot(data = cell_comp_tb,
       aes(x = time, y = count, fill = cell_type)) + geom_area(position = "stack") +
        scale_x_continuous() +
        scale_fill_manual(values = gr3_col) +
        ylab("Cell count") +
        xlab("Time") +
        coord_flip() +
                scale_x_reverse() +
        # geom_vline(xintercept = gr_trans_time)
        ggpubr::theme_pubr() +
        theme(legend.position = "none",
              text = element_text(size = 12))) %>%
        push_pdf(file_name = "fatogram",
                 width = 4.5,
                 height = 2.5,
                 ps = 12,
                 dir = "./plots/fate_map/")

# diff_time = trans_time[x]
# gr_tip_list = data_obj$gr_tip_list
# gr_dd = data_obj$gr_dd
# assertthat::assert_that(!is.null(data_obj$tr_edges_state_tb))
#
# state_include = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[x]]))]
# state_include_d1 = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[gr_dd[[x]][1]]]))]
# state_include_d2 = names(gr_tip_list)[map_lgl(gr_tip_list, function(y) all(y %in% gr_tip_list[[gr_dd[[x]][2]]]))]
#
# edge_tb_diff = mutate(edge_tb, ind = map_lgl(type_path, function(y) any(state_include %in% y)))
# edge_tb_diff = edge_tb_diff %>% filter(from_time <= diff_time &
#                                                to_time > diff_time &
#                                                ind)
# edge_tb_diff = mutate(edge_tb_diff, d1_ind = map_lgl(type_path, function(y) any(state_include_d1 %in% y))) # old version gr_dd[[x]][1] (old version has been overwritten)
# edge_tb_diff = mutate(edge_tb_diff, d2_ind = map_lgl(type_path, function(y) any(state_include_d2 %in% y)))
# assertthat::assert_that(all(!(edge_tb_diff$d1_ind & edge_tb_diff$d2_ind)))
# edge_tb_diff

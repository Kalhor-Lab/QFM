process_split_ratio <- function(eval_tb) {
        # node split is the reference for ordering left and right split
        out_tb = mutate(eval_tb,
                        node_split_order = map_dbl(node_split, function(x) {
                                if (is.null(x)) {
                                        return(NA)
                                }
                                return(x[1])
                        }))
        # mapping sampled ratio to the same order
        out_tb = mutate(out_tb, node_split_sampled_order = map2_dbl(node_split, node_split_sampled, function(x, y) {
                if (is.null(x)) {
                        return(NA)
                }
                return(y[names(x)[1]]/sum(y))
        }))
        out_tb = mutate(out_tb, gr_node_split_order = map2_dbl(node_split, gr_node_size_mapped, function(x, y) {
                if (is.null(y)) {
                        return(NA)
                }
                if (is.null(x)) {
                        return(NA)
                }
                y = y[2:3]
                return((1+y[names(x)[1]])/(2+sum(y)))
        }))
        out_tb
}
process_node_size <- function(eval_tb) {
        out_tb = mutate(eval_tb,
                        gr_node_size_in = map_dbl(gr_node_size, function(x) {
                                if (is.null(x)) {
                                        return(NA)
                                }
                                return(pmax(1, x[1]))
                        }))
        out_tb
}

process_res <- function(eval_tb_list, return_true = F) {
        if (return_true) {
                tb_true_all = bind_rows(map(eval_tb_list, "true"))
                tb_true_all = mutate(tb_true_all, is_resolved = as.numeric(!is.na(node_gr)))
                tb_true_all = process_node_size(tb_true_all)
                tb_true_all = process_split_ratio(tb_true_all)
                tb_true_all = mutate(tb_true_all, log2_node_size = log2(node_size))
                tb_true_all = mutate(tb_true_all, gr_time_trans_error = gr_time_trans - node_time)
                tb_true_all = mutate(tb_true_all, gr_time_est_error = gr_time_est - node_time)
                tb_true_all = mutate(tb_true_all, gr_node_size_logfc = log2(gr_node_size_in / node_size))
                tb_true_all = mutate(tb_true_all, log2_node_sampled = log2(node_size_sampled / node_size))
                return(tb_true_all)
        } else {
                tb_recon_all = bind_rows(map(eval_tb_list, "recon"))
                tb_recon_all = mutate(tb_recon_all, is_resolved = !is.na(node))
                tb_recon_all = process_node_size(tb_recon_all)
                tb_recon_all = process_split_ratio(tb_recon_all)
                tb_recon_all = mutate(tb_recon_all, log2_node_size = log2(node_size))
                tb_recon_all = mutate(tb_recon_all, gr_time_trans_error = gr_time_trans - node_time)
                tb_recon_all = mutate(tb_recon_all, gr_time_est_error = gr_time_est - node_time)
                tb_recon_all = mutate(tb_recon_all, gr_node_size_logfc = log2(gr_node_size_in / node_size))
                tb_recon_all = mutate(tb_recon_all, log2_node_sampled = log2(node_size_sampled / node_size))
                return(tb_recon_all)
        }
}

eval_tb_wrap <- function(eval_tb, exp_params) {
        eval_tb_all = process_res(eval_tb, return_true = F)
        eval_tb_all = left_join(eval_tb_all, select(exp_params, j, big_graph_id, num_tip, sampling, bsum), by = "j")
        eval_tb_all = mutate(eval_tb_all, node_collect_size = map_dbl(node_size_collect, function(x) {
                if (is.null(x)){
                        return(NA)
                } else  {
                        return(x)
                }
        }))
        eval_tb_all = mutate(eval_tb_all, gr_collect_size_used = gr_node_size_in  / node_collect_size)
        eval_tb_all$suff_sampled = factor(eval_tb_all$log2_node_sampled > -2)
        eval_tb_all$high_cov = factor((1/eval_tb_all$gr_collect_size_used) > 2.5)
        eval_tb_all = eval_tb_all %>%
                mutate(log2_node_size = log2(node_size),
                       log2_gr_node_size_in = log2(gr_node_size_in))
        eval_tb_all
}

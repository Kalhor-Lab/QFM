sps_default = function(lin_tr, tip_id, type_func = get_type_from_id) {
        if (is.null(lin_tr)) {
                return(NULL)
        }
        lin_tr$edge.length[lin_tr$edge.length < 0] = 0
        node_tips_list = list_dd_and_tips(lin_tr)$tips
        if (lin_tr$Nnode == 1) {
                return(NULL)
        }
        node_type_list = map(node_tips_list[2:length(node_tips_list)],
                             function(x) table(type_func(x))[tip_id])
        sim_mat = matrix(0,
                         nrow = length(tip_id),
                         ncol = length(tip_id))
        rownames(sim_mat) = colnames(sim_mat) = tip_id
        for (x in node_type_list) {
                x = x[!is.na(x)]
                p_score = 1/2^(sum(x > 0)-1)
                # p_score = 1/sum(x > 0)
                if (sum(x>0) > 1 & sum(x>0) < length(tip_id)) {
                        all_comb = combn(names(x)[x>0], m = 2)
                        for (i in 1:ncol(all_comb)) {
                                sim_mat[all_comb[1, i], all_comb[2, i]] = sim_mat[all_comb[1, i], all_comb[2, i]] + p_score
                                sim_mat[all_comb[2, i], all_comb[1, i]] = sim_mat[all_comb[2, i], all_comb[1, i]] + p_score
                        }
                }
                if (sum(x>0) == 1) {
                        rr = names(x)[x>0]
                        sim_mat[rr, rr] = sim_mat[rr, rr] + 1
                }
                # else {
                #         y = names(x)[x>0]
                #         sim_mat[y, y] = sim_mat[y, y] +  2 * p_score
                # }
        }
        return(sim_mat)
}
process_sps <- function(tr_r) {
        tip_id = sort(unique(get_type_from_id(tr_r$tip.label)))
        sps_mat = sps_default(name_nodes(tr_r), tip_id)
        dmat = 1 - sps_mat/max(sps_mat)
        gr = phangorn::upgma(dmat)
        list(gr = gr, sps = sps_mat)
}

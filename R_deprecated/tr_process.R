process_parallel <- function(pipeline_func, data, type_graph, bp_param, shuffle_input = T) {
        if (shuffle_input) {
                shuffle_indices = sample(length(data), replace = F)
        } else {
                shuffle_indices = 1:length(data)
        }
        out = bplapply(data[shuffle_indices], function(x, run_func, type_graph) {
                run_func(x, type_graph)
        }, BPPARAM = bp_param,
        run_func = pipeline_func,
        type_graph = type_graph)
        return(out[order(shuffle_indices)])
}
tr2sps <- function(tr, type_graph) {
        devtools::load_all("./")
        # if ((length(tr$tip.label) %% 10) != 0) {
        #         return(NA)
        # }
        d = sps_ver1(tr, type_graph$tip_id)
        d
}
plot_sps <- function(sps_mat, tr) {
        Heatmap(sps_mat[tr$tip.label, tr$tip.label],
                cluster_rows = as.dendrogram(tr), cluster_columns = as.dendrogram(tr))
}
sc2tr <- function(sc_mat, type_graph) {
        devtools::load_all("./")
        reconstruct_lineage(sc_mat)
}
sc_pseudo_bulk <- function(x) {
        x_bulk = purrr::map(1:ncol(x), function(j) {
                x_j = barcode_allele_onehot_new(x[, j, drop = F])
                if (is.null(x_j)) {
                        return(NULL)
                }
                out = ghelper::aveMatFac(x_j,
                                         get_type_from_id(rownames(x)))
                out
        })
        x_bulk
}
simulate_bulk_and_sc <- function(type_graph, mut_p, sample_size) {
        require(extraDistr)
        require(matrixStats)
        require(data.table)
        devtools::load_all(".")
        produce_data <- function(type_graph, mut_p, sample_size) {
                gens0 = make_gens(type_graph = type_graph)
                gens1 = sample_size_gens(gens0, type_graph, sample_size = sample_size)
                mut_param = mut_param = do.call(make_mut_param_by_rate_rvec, mut_p)
                message("simluating sc..")
                gens2 = simulate_sc_muts(gens1,
                                         type_graph = type_graph,
                                         mut_param = mut_param)
                message("constructing true tree..")
                tr = construct_true_lineage(gens2, type_graph$tip_id)
                x = get_barcodes(gens2, type_graph$tip_id)
                x = remove_allele_ver(x)

                # message("simulating bulk..")
                # gens3 = simulate_muts(gens1, type_graph, mut_param = mut_param)
                # y = extract_allele_matrix(gens3[type_graph$tip_id])
                # y = lapply(y, collapse_reucr_ver)

                message("pooling pseudo-bulk...")
                x_bulk = purrr::map(1:ncol(x), function(j) {
                        x_j = barcode_allele_onehot_new(x[, j, drop = F])
                        if (is.null(x_j)) {
                                return(NULL)
                        }
                        out = ghelper::sumMatFac(x_j,
                                                 get_type_from_id(rownames(x)))
                        out
                        out
                })
                out = list(tr = tr,
                           sc = x,
                           bulk = x_bulk)
                return(out)
        }
        return(produce_data(type_graph = type_graph,
                            mut_p = mut_p,
                            sample_size = sample_size))
}
process_tr <- function(x, total_time, ...) {
        restrict_tb = estimation_pairwise_restrictions(x)
        gr = name_nodes(restrction_to_graph(restrict_tb, total_time, ...))
        list(restriction_times = restrict_tb,
             gr = gr)
}

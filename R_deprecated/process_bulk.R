# tip_mat_collapsed_fil = lapply(tip_mat_collapsed, filter_allele_matrix, frac_cut = 0.01)
# tip_mat_scaled_all = do.call(cbind,
#                              lapply(tip_mat_collapsed_fil, function(x) clr(normalize_row(x+100, m = 5000), add_one = F)))
# tip_mat_all = do.call(cbind, tip_mat_collapsed_fil)
# allele_drop = paste0("R-", 1:length(mut_p$recur_vec))[which(mut_p$recur_vec > 0.05)]
# tip_mat_dropped_all = do.call(cbind, lapply(tip_mat_collapsed_fil, function(x) {
#         clr(normalize_row(x[, !colnames(x) %in% allele_drop, drop = F]+100, m = 5000))
#         # clr(x[, !colnames(x) %in% allele_drop, drop = F])
# }))


variable_alleles_indices <- function(allele_mat, cutoff) {
        means = log2(1+colMeans(allele_mat))
        vars = log2(1+colVars(allele_mat))
        lowess_fit = lowess(means, vars)
        # col_vec = rep("black", length(vars))
        # col_vec[which((vars - approxfun(lowess_fit)(means)) > 1)] = "red"
        plot(means, vars, col = col_vec)
        lines(lowess_fit, col = "red")
        idx = which((vars - approxfun(lowess_fit)(means)) > cutoff)
        idx
}

bulk_pipeline <- function(allele_matrix_list, params) {
        # allele_matrix_list = dat0_phy[[1]]
        # params = params_list[[2]]
        allele_matrix_list = lapply(1:length(allele_matrix_list), function(j) {
                x = allele_matrix_list[[j]]
                colnames(x) = paste0("element_", j, "_", colnames(x))
                x
        })
        # filter
        allele_matrix_list = lapply(allele_matrix_list, filter_allele_matrix, frac_cut = params$frac_cut)
        element_sel = which(sapply(allele_matrix_list, ncol) > 1)
        if (length(element_sel) == 0) {
                return(NULL)
        }
        allele_matrix_list = allele_matrix_list[element_sel]
        allele_matrix_list = lapply(allele_matrix_list, function(x) normalize_row(x+params$normalize_pseudo_count,
                                                                                  m = 1.))
        if (params$var_sel) {
                temp_combined = do.call(cbind, allele_matrix_list)
                allele_sel = colnames(temp_combined)[variable_alleles_indices(temp_combined, cutoff = params$var_cut)]
                }
        if (params$clr) {
                allele_matrix_list = lapply(allele_matrix_list, function(x) compositions::clr(x))
                # allele_matrix_list = lapply(allele_matrix_list, function(x) compositions::ilr(x))
                }
        allele_matrix_combined = do.call(cbind, allele_matrix_list)
        if (params$var_sel) {
                allele_matrix_combined = allele_matrix_combined[, allele_sel]
                }
        if (params$scale_data) {
                allele_matrix_combined = scale(allele_matrix_combined)
        }
        # if (!is.na(params$pc_ndim)) {
        #         pr_out = prcomp(allele_matrix_combined)
        #         allele_matrix_combined = allele_matrix_combined %*% pr_out$rotation[, 1:params$pc_ndim]
        # }
        allele_matrix_combined
}
params_list = list("raw" = list(frac_cut = 0.01,
                                normalize_pseudo_count = 0,
                                clr = F,
                                var_sel = F,
                                var_cut = NA,
                                scale_data = F,
                                dist_method = "manhattan"),
                   "raw_scale" = list(frac_cut = 0.01,
                                      normalize_pseudo_count = 0,
                                      clr = F,
                                      var_sel = F,
                                      var_cut = NA,
                                      scale_data = T,
                                      dist_method = "manhattan"),
                   "clr_250" = list(frac_cut = 0.01,
                                    normalize_pseudo_count = 250,
                                    clr = T,
                                    var_sel = F,
                                    var_cut = NA,
                                    scale_data = F,
                                    dist_method = "euclidean"),
                   "clr_scale_250" = list(frac_cut = 0.01,
                                          normalize_pseudo_count = 250,
                                          clr = T,
                                          var_sel = F,
                                          var_cut = NA,
                                          scale_data = T,
                                          dist_method = "euclidean"))
# proc_list = lapply(params_list, function(x) bulk_pipeline(dat0_g2[[1]], x))
# bulk_pipeline(dat0_g2, params_list[[2]])
# phy1 = upgma(dist(proc_mat))
# Heatmap(proc_mat,
#         cluster_rows = as.dendrogram(phy1))

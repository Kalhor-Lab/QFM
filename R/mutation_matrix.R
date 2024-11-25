# function for processing mutation matrix
# imputation
impute_characters <- function (mat, nrounds = 50, max_depth = 4) {
        op = options(na.action = "na.pass")
        on.exit(options(op))
        col_div = apply(mat, 2, function(x) length(unique(x[!is.na(x)])))
        im_indices = which(col_div > 1)
        keep_indices = which(col_div == 1)
        mat_keep = rlang::duplicate(mat[, keep_indices, drop = F])
        if (ncol(mat_keep) > 0) {
                for (i in 1:ncol(mat_keep)) {
                        char_vec = mat_keep[, i]
                        char_vec_unique = unique(char_vec[!is.na(char_vec)])
                        assertthat::assert_that(length(char_vec_unique) ==
                                                        1)
                        mat_keep[is.na(char_vec), i] = char_vec_unique
                }
        }
        assertthat::assert_that(all(!is.na(mat_keep)))
        mat_im = as.data.frame(rlang::duplicate(mat[, im_indices,
                                                    drop = F]))
        element_id = 1:ncol(mat_im)
        na_frac = map_dbl(element_id, function(j) {
                mean(is.na(mat_im[, j]))
        })
        element_id = element_id[order(na_frac)]
        for (i in element_id) {
                # message(i)
                if (na_frac[i] == 0) {
                        next
                }
                design_mat = model.matrix(~., mat_im[, -i])[, -1]
                # reduce design mat
                design_mat = design_mat[, colSums(design_mat, na.rm = T) > 1]
                y_fac = factor(as.matrix(mat_im[i])[, 1])
                y_vec = as.numeric(y_fac)
                bst <- xgboost(data = as.matrix(design_mat[!is.na(y_vec), , drop = F]),
                               label = y_vec[!is.na(y_vec)] - 1, max_depth = max_depth,
                               eta = 1, nthread = 12, nrounds = nrounds, objective = "multi:softmax",
                               eval_metric = "mlogloss", num_class = max(y_vec[!is.na(y_vec)]),
                               verbose = 0)
                # suppressWarnings()
                y_pred = levels(y_fac)[predict(bst, xgb.DMatrix(as.matrix(design_mat[is.na(y_vec), , drop = F]))) + 1]
                mat_im[is.na(y_vec), i] = y_pred
        }
        mat_im = as.matrix(mat_im)
        rownames(mat_im) = rownames(mat)
        mat_out = cbind(mat_im, mat_keep)
        mat_out = mat_out[, colnames(mat)]
        assertthat::assert_that(!anyNA(mat_out))
        mat_out
}
remove_uniform_id <- function(mat, abund_thres = 0) {
        # count diversity ignoring single allele
        id_diversity = apply(mat, 2, function(x) {
                x_tab = table(x[!is.na(x)])
                x_tab = x_tab[x_tab > abund_thres]
                return(length(x_tab))
        })
        mat[, id_diversity > 1]
}
filter_cells <- function(mat, cutoff = 2, abund_thres = 1) {
        # make a list of informative spacers
        spacer_ind_mat = do.call(cbind, map(1:ncol(mat), function(j) {
                x = mat[, j]
                x_tab = table(x[!is.na(x)])
                x_tab = x_tab[x_tab > abund_thres]
                if (length(x_tab) <= 1) {
                        return(rep(F, nrow(mat)))
                }
                allele_freq = sort(table(mat[, j]))
                out_spacers = names(allele_freq)[allele_freq > 1]
                out_spacers = out_spacers[out_spacers != "0"]
                mat[, j] %in% out_spacers
        }))
        mat[rowSums(spacer_ind_mat) > cutoff, ]
}
filter_noise_spacers <- function(tb, noise_factor = 4) {
        if (nrow(tb) > 1) {
                tb = arrange(tb, desc(count))
                tb_count = tb$count
                if (tb_count[1] > sum(tb_count[2:length(tb_count)]) * noise_factor) {
                        return(tb[1, ])
                } else {
                        return(tb)
                }
        } else {
                return(tb)
        }
}
make_recur_vec_list <- function(allele_prob_tb) {
        allele_prob_nest = allele_prob_tb %>%
                nest(data = -site)
        map(allele_prob_nest$data, function(tb) {
                setNames(tb$prob, tb$mutation)
        })
}

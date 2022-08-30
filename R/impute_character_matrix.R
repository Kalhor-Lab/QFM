im_chr = function (mat, nrounds = 20, max_depth = 4)
{
        if (is.null(colnames(mat))) {
                colnames(mat) = paste0("el_", 1:ncol(mat))
        }
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
                mean(is.na(mat[, j]))
        })
        element_ord = order(na_frac)
        element_id = element_id[element_ord]
        for (i in element_id) {
                print(na_frac[i])
                if (na_frac[i] > 0) {
                        message(i)
                        design_mat = model.matrix(~., mat_im[, -i])[, -1]
                        y_fac = factor(as.matrix(mat_im[i])[, 1])
                        y_vec = as.numeric(y_fac)
                        suppressWarnings(bst <- xgboost(data = design_mat[!is.na(y_vec),
                        ], label = y_vec[!is.na(y_vec)] - 1, max_depth = max_depth,
                        eta = 1, nthread = 12, nrounds = nrounds, objective = "multi:softmax",
                        eval_metric = "mlogloss", num_class = max(y_vec[!is.na(y_vec)]),
                        verbose = 0))
                        y_pred = levels(y_fac)[predict(bst, design_mat[is.na(y_vec), , drop = F
                        ]) + 1]
                        mat_im[is.na(y_vec), i] = y_pred
                } else {
                        print('skiped')
                }
        }
        mat_im = as.matrix(mat_im)
        rownames(mat_im) = rownames(mat)
        mat_out = cbind(mat_im, mat_keep)
        mat_out = mat_out[, colnames(mat)]
        mat_out
}

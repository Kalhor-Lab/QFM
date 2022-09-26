loglike_func_case1 <- function(t, lambda, a, t_total) {
        e_lt = exp(-lambda * t)
        e_lt_total = exp(-lambda * t_total)
        e_ratio = e_lt_total / e_lt
        log(
                e_ratio * ((1 - e_lt) * a)^2 + (1 - e_ratio) * a
        )
}
loglike_func_case2 <- function(t, i_total, lambda) {
        e_lt = exp(-lambda * t)
        log(e_lt) * (1 - i_total) + log(1 - e_lt) * i_total
}
est_func_case1 <- function(t, lambda, a, t_total) {
        e_lt = exp(lambda * t)
        e_lt_total = exp(lambda * t_total)
        lambda * ((a - 1) * e_lt^2 - a) /
                ((e_lt_total/e_lt + a - 1) * e_lt^2 - 2 * a * e_lt + a)
}
est_func_case1a <- function(t, lambda, t_total) {
        # version of case 1 that ignores independent occurrences
        e_lt = exp(lambda * t)
        e_lt_total = exp(lambda * t_total)
        lambda / (1 - e_lt_total/e_lt)
}
est_func_case2 <- function(t, i_total, lambda) {
        e_lt_inv = exp(lambda * t)
        lambda * ((i_total - 1) * e_lt_inv + 1) / (e_lt_inv - 1)
}
#' Estimate coalescence time (not to be exported)
estimate_coal_time <- function(b1, b2, mut_p, total_coal_time, min_alele_prob = 0.) {
        avail_indices = !is.na(b1) & !is.na(b2)
        if (!any(avail_indices)) return(total_coal_time)
        b1 = b1[avail_indices]
        b2 = b2[avail_indices]
        mut_p = list(mut_rate = mut_p$mut_rate[avail_indices],
                     recur_vec_list = mut_p$recur_vec_list[avail_indices])
        case_ind = (b1 == b2 & b1 != "0" & b2 != "0")
        a_vec = map2_dbl(mut_p$recur_vec_list[case_ind], b1[case_ind], function(r_vec, m) {
                r_vec[m]
        })
        a_ind = a_vec > min_alele_prob
        est_func <- function(t_val) {
                c1 = est_func_case1(t = t_val,
                                    lambda = mut_p$mut_rate[case_ind][a_ind],
                                    a = a_vec[a_ind],
                                    t_total = total_coal_time)
                c1a = est_func_case1a(t = t_val,
                                      lambda = mut_p$mut_rate[case_ind][!a_ind],
                                      t_total = total_coal_time)
                c2 = est_func_case2(t = t_val,
                                    i_total = (b1[!case_ind] != "0") + (b2[!case_ind] != "0"),
                                    lambda = mut_p$mut_rate[!case_ind])
                sum(c(c1, c1a, c2))
        }
        if (est_func(total_coal_time) > 0 & est_func(0.1) > 0) {
                return(total_coal_time)
        }
        if (est_func(total_coal_time) < 0 & est_func(0.1) < 0) {
                return(0)
        }
        t_est = uniroot(est_func, interval = c(0.1, total_coal_time))$root
        t_est
}
#' Convert pairwise distanes in a data.frame format to a matrix format
dist_df2mat <- function(dist_df) {
        tips = unique(c(dist_df$N1, dist_df$N2))
        dist_mat = matrix(NA, length(tips), length(tips))
        colnames(dist_mat) = rownames(dist_mat) = tips
        dist_mat[cbind(dist_df$Var1, dist_df$Var2)] = dist_df$dist
        dist_mat[cbind(dist_df$Var2, dist_df$Var1)] = dist_df$dist
        dist_mat
}
#' Estimate mutation rate from single cell character matrix
estimate_mut_rate <- function(sc_mat, t_total) {
        p_est = map_dbl(1:ncol(sc_mat), function(j) {
                p = mean(sc_mat[, j] != "0", na.rm = T)
                p
        })
        lambda_est = map_dbl(p_est, function(p) {
                -log(1 - p + 1e-10)/t_total
        })
        lambda_est
}
#' Generate uniform allele emergence probabilities for Phylotime inference.
default_allele_prob <- function(sc_mat) {
        purrr::map(1:ncol(sc_mat), function(k) {
                muts = sort(unique(sc_mat[, k]))
                muts = muts[muts != "0"]
                obs_prob = rep(1, length(muts))/length(muts)
                names(obs_prob) = muts
                obs_prob
        })
}
#' Estimate allele emergence probabilities from single cell character matrix
estimate_allele_prob <- function(sc_mat) {
        purrr::map(1:ncol(sc_mat), function(k) {
                muts = sort(unique(sc_mat[, k]))
                muts = muts[muts != "0"]
                obs_prob = table(sc_mat[, k])[muts]/sum(sc_mat[, k] != "0", na.rm = T)
                obs_prob
        })
}
# Estimate a `mut_p` object from single cell character matrix
estimate_mut_p <- function(sc_mat, t_total) {
        list(mut_rate = estimate_mut_rate(sc_mat, t_total),
             recur_vec_list = default_allele_prob(sc_mat))
}
#' Wrapper of UPGMA in base R
upgma <- function(D, ...) {
        DD <- as.dist(D)
        hc <- hclust(DD, method = "average", ...)
        result <- as.phylo(hc)
        result <- reorder(result, "postorder")
        result
}

#' Run Phylotime for reconstructing time-scaled phylogeny based on lineage barcodes
#' @param sc_mat character matrix cell x barcoding sites
#' @param t_total total amount of time since barcode activation
#' @param return_dist whether to return pairwise distance between cells instead, by default, a `phylo` object is returned
#' @return a `phylo` object
phylotime <- function(sc_mat, t_total = 1., mut_p = NULL, return_dist = F, parallel = T) {
        assertthat::assert_that(!is.null(rownames(sc_mat)))
        # plan(multisession, workers = 12)
        if (is.null(mut_p)) {
                mut_p = estimate_mut_p(sc_mat, t_total)
        }
        dist_df = as_tibble(expand.grid(1:nrow(sc_mat), 1:nrow(sc_mat)))
        dist_df = dist_df[dist_df$Var1 < dist_df$Var2, ]
        # p <- progressor(along = 1:nrow(dist_df))
        # tt = proc.time()
        if (parallel == T) {
                dist_df$coal_time = future_map_dbl(1:nrow(dist_df), function(k) {
                        i = dist_df$Var1[k]
                        j = dist_df$Var2[k]
                        estimate_coal_time(sc_mat[i, ],
                                           sc_mat[j, ],
                                           mut_p,
                                           total_coal_time = t_total,
                                           min_alele_prob = 0.0)
                }, .progress = T)
        } else {
                dist_df$coal_time = map_dbl(1:nrow(dist_df), function(k) {
                        i = dist_df$Var1[k]
                        j = dist_df$Var2[k]
                        estimate_coal_time(sc_mat[i, ],
                                           sc_mat[j, ],
                                           mut_p,
                                           total_coal_time = t_total,
                                           min_alele_prob = 0.0)
                })
        }
        # print(proc.time() - tt)
        dist_df$dist = dist_df$coal_time * 2
        dist_df$N1 = rownames(sc_mat)[dist_df$Var1]
        dist_df$N2 = rownames(sc_mat)[dist_df$Var2]
        if (return_dist) {
                return(dist_df)
        } else {
                dmat = dist_df2mat(dist_df)
                ord_indices = sample(nrow(dmat), replace = F)
                tr = phangorn::upgma(dmat[ord_indices, ord_indices])
                # tr = upgma(dmat[ord_indices, ord_indices])
                tr = name_nodes(tr)
                return(tr)
        }
}

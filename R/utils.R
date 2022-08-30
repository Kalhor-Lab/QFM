#' deprecated
#' Assign counts based on total and fraction, rounded to integer
#' Ensure that each category gets some minimum of counts
# round_frac <- function(total, prob, min_val = rep(0, length(prob))) {
#         assertthat::assert_that(length(min_val) == length(prob))
#         assertthat::assert_that(sum(min_val) <= total)
#         n_ind = which.max(prob)[1]
#         rounded = round(total * prob[-n_ind])
#         rounded = pmax(rounded, min_val[-n_ind])
#         out = numeric(length(prob))
#         out[-n_ind] = rounded
#         out[n_ind] = total - sum(rounded)
#         assertthat::assert_that(out[n_ind] >= min_val[n_ind])
#         names(out) = names(prob)
#         out
# }

#' Scale counts based on a factor, rounded to integer, counts must be at least 1
scale_counts <- function(counts, z) {
        assertthat::assert_that(length(counts) == length(z))
        if (length(counts) == 1) {
                return(counts)
        }
        # TODO: find correct way to scale when counts are close to zero
        # Important sample size must be at least 1, exclude from scaling for now
        val = counts * z
        # stochastic rounding
        rounded = sapply(val[-1], rstoc)
        n_trial = 0
        # force a valid rounding to happen, temporary fix
        while (!(sum(rounded) <= (sum(counts) - 1))) {
                rounded = sapply(val[-1], rstoc)
                n_trial = n_trial + 1
                if (n_trial > 10) {
                        print(counts)
                        print(z)
                        stop('trouble finding scaling result')
                }
        }
        return(c(sum(counts)-sum(rounded), rounded))
}
#' Stochastic roundind, results all at least 1
rstoc <- function(x){
        q <- abs(x - trunc(x))
        adj <- sample(0:1, size = 1, prob = c(1 - q, q))
        if(x < 0) adj <- adj * -1
        # at least 1
        max(1, trunc(x) + adj)
}
#' Check if a vector is all the same value
is_zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
        if (length(x) == 1) return(TRUE)
        x <- range(x) / mean(x)
        isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
#' Check if a vector is non-decreasing
is_non_decs <- function(vec) {
        all(diff(vec) >= 0)
}
normalize_row <- function(x, m = 1000) {
        total = rowSums(x)
        total[total==0] = 1
        x / total * m
}
clr <- function(x, add_one = T) {
        if (add_one) {
                x <- log(x+1)
        } else {
                x <- log(x)
        }
        x <- sweep(x, 1, rowMeans(x), "-")
        x
}
subset_mut_p <- function(mut_p, indices) {
        list(mut_rate = mut_p$mut_rate[indices],
             active_time = mut_p$active_time,
             recur_prob = mut_p$recur_prob[indices],
             recur_vec_list = mut_p$recur_vec_list[indices])
}
concat_mut_p <- function(mut_p1, mut_p2) {
        assertthat::assert_that(all.equal(mut_p1$active_time, mut_p2$active_time))
        list(mut_rate = c(mut_p1$mut_rate, mut_p2$mut_rate),
             active_time = mut_p1$active_time,
             recur_prob = c(mut_p1$recur_prob, mut_p2$recur_prob),
             recur_vec_list = c(mut_p1$recur_vec_list, mut_p2$recur_vec_list))
}
concat_mut_p_list <- function(mut_p_list) {
        purrr::reduce(mut_p_list, concat_mut_p)
}

push_pdf <- function(g, file_name, width = 4, height = 4, ps = 12, dir = "../LTModelPlots/", open_file = T) {
        pdf(paste0(dir, file_name, ".pdf"), width = width, height = height, pointsize = ps)
        print(g)
        dev.off()
        if (open_file) {
                file_path = gsub("/", "\\\\", paste0(dir, file_name, ".pdf"))
                shell.exec(file_path)
        }
}
push_png <- function(g, file_name, width = 4, height = 4, ps = 12, res = 300, dir = "../LTModelPlots/", open_file = T) {
        png(paste0(dir, file_name, ".png"), res = res, units = "in", width = width, height = height, pointsize = ps)
        print(g)
        dev.off()
        if (open_file) {
                file_path = gsub("/", "\\\\", paste0(dir, file_name, ".png"))
                shell.exec(file_path)
        }
}
make_heatmap_col <- function(max_val, mid_val = max_val/2, min_val = 0, col = c("blue", "white", "red")) {
        circlize::colorRamp2(breaks = c(min_val, mid_val, max_val), colors = col)
}
plot_col_vec <- function(col_vec) {
        tibble(name = names(method_col),
               val = 1.) %>%
                ggbarplot(x = "name", y = "val",
                          color = NA, fill = "name") + scale_fill_manual(values = col_vec)
}
compute_total_ranking <- function(gr) {
        list_out = list_dd_and_tips(name_nodes(gr))
        gr_dd = list_out$dd
        gr_tip_lists = list_out$tips
        prod(map_dbl(gr_dd, function(x) {
                s1 = pmax(1, length(gr_tip_lists[[x[1]]])) - 1
                s2 = pmax(1, length(gr_tip_lists[[x[2]]])) - 1
                factorial(s1 + s2) / factorial(s1) / factorial(s2)
        }))
}
# deprecated
# sample_node_ranking <- function(gr) {
#         list_out = list_dd_and_tips(name_nodes(gr))
#         gr_dd = list_out$dd
#         gr_tip_lists = list_out$tips
#         # genrate the draw for each node first
#         node_shuffle = map(gr_dd, function(x) {
#                 s1 = pmax(1, length(gr_tip_lists[[x[1]]])) - 1
#                 s2 = pmax(1, length(gr_tip_lists[[x[2]]])) - 1
#                 sample(c(rep(x[1], s1), rep(x[2], s2)))
#         })
#         root_node = names(gr_dd)[1]
#         ranking = rep(root_node, gr$Nnode)
#         assign_dd <- function(dd) {
#                 if (dd %in% gr$tip.label) {
#                         return()
#                 } else {
#                         dd_indices = which(ranking == dd)
#                         if (length(dd_indices) <= 1) {
#                                 return()
#                         } else {
#                                 ranking[dd_indices[-1]] <<- node_shuffle[[dd]]
#                                 assign_dd(gr_dd[[dd]][1])
#                                 assign_dd(gr_dd[[dd]][2])
#                         }
#                 }
#         }
#         assign_dd(root_node)
#         ranking
# }
plot_legend <- function(g) {
        legend <- cowplot::get_legend(g)
        grid.newpage()
        grid.draw(legend)
}
#' compare if two strings are identical other than N characters
check_n_identity <- function(str_n, str_ref) {
        str_n_vec = strsplit(str_n, "")[[1]]
        str_ref_vec = strsplit(str_ref, "")[[1]]
        if (length(str_n_vec) != length(str_ref_vec)){
                return(F)
        }
        N_ind = map_lgl(str_n_vec, function(x) x == "N")
        if (all(str_n_vec[!N_ind] == str_ref_vec[!N_ind])) {
                return(T)
        } else {
                return(F)
        }
}
white_col_mapper <- function(types) {
        out = rep("white", length(types))
        names(out) = types
        return(out)
}
calc_bsum <- function(phy) {
        sum(apply(balance(phy), 1, function(x) abs(x[1] - x[2])))
}

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
        # force a valid rounding to happen, to be improved
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
#' Round a number stochasticly, returns integer greater than or equal to one.
#' @param x integer
rstoc <- function(x){
        q <- abs(x - trunc(x))
        adj <- sample(0:1, size = 1, prob = c(1 - q, q))
        if(x < 0) adj <- adj * -1
        # at least 1
        max(1, trunc(x) + adj)
}
#' Check if a vector is all the same value
#' @param x An input vector.
#' @param tol tolerance when comparing elements of the vector
is_zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
        if (length(x) == 1) return(TRUE)
        x <- range(x) / mean(x)
        isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
#' Check if a vector is non-decreasing
#' @param vec An input vector.
is_non_decs <- function(vec) {
        all(diff(vec) >= 0)
}
#' Normalize rows of a matrix to be the same sum `m`
normalize_row <- function(x, m = 1000) {
        total = rowSums(x)
        total[total==0] = 1
        x / total * m
}
# Centered log ratio transform
clr <- function(x, add_one = T) {
        if (add_one) {
                x <- log(x+1)
        } else {
                x <- log(x)
        }
        x <- sweep(x, 1, rowMeans(x), "-")
        x
}
#' Subset a `mut_p` object to a number of its elements
subset_mut_p <- function(mut_p, indices) {
        list(mut_rate = mut_p$mut_rate[indices],
             active_time = mut_p$active_time,
             recur_prob = mut_p$recur_prob[indices],
             recur_vec_list = mut_p$recur_vec_list[indices])
}
#' Concatenate two `mut_p` (mutagenesis parameters) objects.
concat_mut_p <- function(mut_p1, mut_p2) {
        assertthat::assert_that(all.equal(mut_p1$active_time, mut_p2$active_time))
        list(mut_rate = c(mut_p1$mut_rate, mut_p2$mut_rate),
             active_time = mut_p1$active_time,
             recur_prob = c(mut_p1$recur_prob, mut_p2$recur_prob),
             recur_vec_list = c(mut_p1$recur_vec_list, mut_p2$recur_vec_list))
}
#' Concatenate a list of `mut_p` (mutagenesis parameters) objects.
concat_mut_p_list <- function(mut_p_list) {
        purrr::reduce(mut_p_list, concat_mut_p)
}
#' Push a ggplot object to a PDF file
push_pdf <- function(g, file_name, width = 4, height = 4, ps = 12, dir = "./", open_file = T) {
        pdf(paste0(dir, file_name, ".pdf"), width = width, height = height, pointsize = ps)
        print(g)
        dev.off()
        if (open_file) {
                file_path = gsub("/", "\\\\", paste0(dir, file_name, ".pdf"))
                shell.exec(file_path)
        }
}
#' Push a ggplot object to a PNG file
push_png <- function(g, file_name, width = 4, height = 4, ps = 12, res = 300, dir = "../LTModelPlots/", open_file = T) {
        png(paste0(dir, file_name, ".png"), res = res, units = "in", width = width, height = height, pointsize = ps)
        print(g)
        dev.off()
        if (open_file) {
                file_path = gsub("/", "\\\\", paste0(dir, file_name, ".png"))
                shell.exec(file_path)
        }
}
#' Make a color palette for Heatmap visualization
#' @importFrom circlize colorRamp2
make_heatmap_col <- function(max_val, mid_val = max_val/2, min_val = 0, col = c("blue", "white", "red")) {
        circlize::colorRamp2(breaks = c(min_val, mid_val, max_val), colors = col)
}
#' Plot a color vector for visualization
plot_col_vec <- function(col_vec) {
        tibble(name = names(method_col),
               val = 1.) %>%
                ggbarplot(x = "name", y = "val",
                          color = NA, fill = "name") + scale_fill_manual(values = col_vec)
}
# compute_total_ranking <- function(gr) {
#         list_out = list_dd_and_tips(name_nodes(gr))
#         gr_dd = list_out$dd
#         gr_tip_lists = list_out$tips
#         prod(map_dbl(gr_dd, function(x) {
#                 s1 = pmax(1, length(gr_tip_lists[[x[1]]])) - 1
#                 s2 = pmax(1, length(gr_tip_lists[[x[2]]])) - 1
#                 factorial(s1 + s2) / factorial(s1) / factorial(s2)
#         }))
# }
#' Plot legend on a new canvas given a ggplot object
#' @importFrom cowplot get_legend
plot_legend <- function(g) {
        legend <- cowplot::get_legend(g)
        grid.newpage()
        grid.draw(legend)
}
#' Compre if two strings are identical other than 'N' characters
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
#' Generate a 'white' color vector for plotting
white_col_mapper <- function(types) {
        out = rep("white", length(types))
        names(out) = types
        return(out)
}
#' Compute BSUM or Colless index value
#' @param phy An `phylo` object.
#' @importFrom ape balance
calc_bsum <- function(phy) {
        sum(apply(balance(phy), 1, function(x) abs(x[1] - x[2])))
}

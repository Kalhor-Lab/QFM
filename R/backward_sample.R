# log choose
log_choose <- function(n, x) {
        return(lgamma(n+1) - lgamma(n-x+1) - lgamma(x+1))
}
#' the stochastic coal function
#' sample the number of doublet branches when dividing n cells into 2*n cell
sample_doublet_branch <- function(n, k) {
        # TODO: fix for very large n
        if (n > 1e8) {
                return(0)
        }
        # formula
        # (choose(n, z) * choose(n - z, k - 2*z) * 2^(k - 2*z))/choose(2*n, k)
        z = 0:floor(k/2)
        log_p = log_choose(n, z) +
                log_choose(n-z, n-k+z) +
                log(2)*(k-2*z) -
                log_choose(2*n, k)
        err = abs(1-sum(exp(log_p)))
        assertthat::assert_that(err < 1e-3, msg = paste0("error: ", err, ", n: ", n, ", k: ", k))
        return(sample(z, size = 1, prob = exp(log_p)))
}
#' the asymmetric coal function
sample_doublet_twotype <- function(n, k1, k2) {
        z = 0:pmin(k1, k2)
        # formula:
        # choose(n, z) * choose(n-z, k1 - z) * choose(n - k1, k2 - z) / (choose(n, k1) * choose(n, k2))
        log_p  = log_choose(n, z) + log_choose(n - z, k1 - z) + log_choose(n - k1, k2 - z) - log_choose(n, k1) - log_choose(n, k2)
        assertthat::assert_that(abs(1-sum(exp(log_p))) < 1e-3)
        return(sample(z, size = 1, prob = exp(log_p)))
}

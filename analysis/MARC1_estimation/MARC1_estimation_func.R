library(parallel)
library(optimParallel)
margin_func <- function(theta, n_vec, k_vec) {
        local_logprob <- function(n, k, theta) {
                p0 = dbinom(0, size = n, theta)
                return((1-k)*log(p0) + k*log(1-p0))
        }
        all_neglogprob <- function(n_vec, k_vec, theta) {
                assertthat::assert_that(length(n_vec) == length(k_vec))
                - sum(sapply(1:length(n_vec), function(j) {
                        local_logprob(n_vec[j], k_vec[j], theta)
                }))
        }
        all_neglogprob(n_vec, k_vec, theta)
}
logsumexp <- function (x) {
        y = max(x)
        y + log(sum(exp(x - y)))
}
softmax <- function (x) {
        exp(x - logsumexp(x))
}
joint_func <- function(theta_vec, n_vec, k_mat) {
        assertthat::assert_that(ncol(k_mat) == length(theta_vec))
        logsumexp <- function (x) {
                y = max(x)
                y + log(sum(exp(x - y)))
        }
        softmax <- function (x) {
                exp(x - logsumexp(x))
        }
        local_logprob <- function(n, k, theta) {
                p0 = dbinom(0, size = n, theta)
                return((1-k)*log(p0) + k*log(1-p0))
        }
        all_neglogprob <- function(n_vec, k_vec, theta) {
                assertthat::assert_that(length(n_vec) == length(k_vec))
                - sum(sapply(1:length(n_vec), function(j) {
                        local_logprob(n_vec[j], k_vec[j], theta)
                }))
        }
        prob_vec = softmax(theta_vec)
        sum(sapply(1:ncol(k_mat), function(j) {
                all_neglogprob(n_vec, k_mat[, j], prob_vec[j])
        }))
}
calc_obs_prob <- function(n_vec, p_vec) {
        1 - sapply(1:length(p_vec), function(k) {
                logprob0 = dbinom(0, size = n_vec, prob = p_vec[k], log = T)
                prob0 = exp(sum(logprob0))
        })
}
calc_recur_prob <- function(n_vec, p_vec) {
        1 - sapply(1:length(p_vec), function(k) {
                logprob0 = dbinom(0, size = n_vec, prob = p_vec[k], log = T)
                logprob1 = dbinom(1, size = n_vec, prob = p_vec[k], log = T)
                prob0 = exp(sum(logprob0))
                prob1 = sum(exp(sapply(1:length(n_vec),
                                       function(j) sum(c(logprob0[-j], logprob1[j])))))
                prob0 + prob1
        })
}

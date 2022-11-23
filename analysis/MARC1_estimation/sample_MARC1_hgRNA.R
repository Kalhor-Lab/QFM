source("analysis/MARC1_estimation//load_MARC1_data.R")
# sample fast and mid hgRNAs
id_prob = readRDS("./intermediate_data/MARC1_id_prob.rds")
# indelphi predictied allele probabilities
indelphi_tb = readRDS("./intermediate_data/MARC1_indelphi_predicted.rds")
indelphi_tb = indelphi_tb %>% mutate(recur_vec = purrr::map(indelphi, function(x) {
        p = x$prob
        p = sort(p/sum(p), decreasing = T)
        names(p) = paste0("R-", 1:length(p))
        p
}))
id_mean_mr = id_prob %>% nest(mr_data = -id) %>%
        transmute(id = id,
                  mutation_rate =  map_dbl(mr_data, function(x) {
                          p = exp(x$scaled_log_prob)
                          sum(x$mutation_rate*p/sum(p))
                  }))
# sample_mutp_mid_or_fast <- function(n, active_time) {
#         id_sampled = sample(table_s3$identifier[table_s3$Class %in% c("mid", "fast")], size = n, replace = T)
#         id_sampled_params = tibble(id = id_sampled) %>% left_join(id_mean_mr, by = "id")
#         id_sampled_params = id_sampled_params %>%
#                 left_join(select(indelphi_tb, id, recur_vec))
#         mut_p = list(mut_rate = id_sampled_params$mutation_rate,
#                      active_time = list(active_time),
#                      recur_prob = rep(1., nrow(id_sampled_params)),
#                      recur_vec_list = id_sampled_params$recur_vec)
#         mut_p
# }
id_params = tibble(id = table_s3$identifier[table_s3$Class != "inactive"]) %>% left_join(id_mean_mr, by = "id")
id_params = id_params %>%
        left_join(select(indelphi_tb, id, recur_vec))
id_params$hgRNA_class = table_s3$Class[match(id_params$id, table_s3$identifier)]
hgRNA_pool = list(mut_rate = id_params$mutation_rate,
                  recur_prob = rep(1., nrow(id_params)),
                  recur_vec_list = id_params$recur_vec,
                  hgRNA_class = id_params$hgRNA_class)

sample_mutp_mid_or_fast <- function(n, active_time) {
        hgrna_indices = sample(which(hgRNA_pool$hgRNA_class  %in% c("mid", "fast")), size = n, replace = T)
        out = subset_mut_p(hgRNA_pool, indices = hgrna_indices)
        out$active_time = list(active_time)
        out
}
sample_mutp <- function(n, active_time, hgrna_class) {
        hgrna_indices = sample(which(hgRNA_pool$hgRNA_class  %in% hgrna_class), size = n, replace = T)
        out = subset_mut_p(hgRNA_pool, indices = hgrna_indices)
        out$active_time = list(active_time)
        out
}




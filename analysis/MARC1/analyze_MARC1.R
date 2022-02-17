library(ComplexHeatmap)
library(ggpubr)
library(patchwork)

# process mouse embryo time course and compare it to simulated results
source("./analysis/MARC1/load_MARC1_data.R")

sim_obs_col = c("Simulated" = "#bd00ff",
                "Unobserved" = "#a0937d",
                "Observed" = "#00B8FF")

m_rate_seq = 10^(seq(-3.2, 0.5, length = 50))
frac_cut = 0.01
get_sim_diversity <- function(x, frac = 0.01) {
        z = x[1, ]
        sum((z>(sum(z)*frac_cut)) & names(z) != "0")
}
logit <- function(p) {
        log(p/(1-p+exp(-20))+exp(-20))
}
check_equal <- function(x) {
        assertthat::assert_that(all(x == x[1]))
        return(x[1])
}

dir0 = './data/MARC1/filtered_pairs/200401-FilteredPairs-PB7/'
dir1 = './data/MARC1/filtered_pairs/200401-FilteredPairs-PB3/'
files0 = list.files(dir0, pattern = 'filteredpairs*\\.txt')
files1 = list.files(dir1, pattern = 'filteredpairs*\\.txt')
# identify parent files
parent_idx0 = grep("parent", files0)
assertthat::assert_that(length(parent_idx0) == 1)
parent_idx1 = grep("parent", files1)
assertthat::assert_that(length(parent_idx1) == 1)
parent_data0 = load_count_file(paste0(dir0, files0[parent_idx0])) %>% rename(parent_spacer = spacer)
parent_data1 = load_count_file(paste0(dir1, files1[parent_idx1])) %>% rename(parent_spacer = spacer)

files0 = files0[-parent_idx0]
files1 = files1[-parent_idx1]
parent_data = rbind(parent_data0, parent_data1)

sample_match0 = str_match(files0, '(E\\d+p\\d|P21)-(((?:t|c)\\d+)(((?:emb)|m)\\d+)|[cj]\\d+m\\d+)(?:-(.*))*_filteredpairs\\.txt')
sample_match1 = str_match(files1, '(E\\d+p\\d|P21)-(((?:t|c)\\d+)(((?:emb)|m)\\d+)|[cj]\\d+m\\d+)-(.*)_filteredpairs\\.txt')
sample_match0[, 1] = paste0(dir0, sample_match0[, 1])
sample_match1[, 1] = paste0(dir1, sample_match1[, 1])
sample_match_all = rbind(sample_match0, sample_match1)

# Group by unique embryos
all_samples = tibble(file = sample_match_all[, 1],
                     emb_name = paste0(sample_match_all[, 2], "_", sample_match_all[, 4], "_", sample_match_all[, 5]),
                     is_embyro = grepl("E", sample_match_all[, 2]),
                     time = ifelse(grepl("E", sample_match_all[, 2]), as.numeric(str_replace(
                             str_replace(sample_match_all[, 2], "E", ""),
                             "p", ".")), 20 + 21),
                     mother = sample_match_all[, 4],
                     emb_id = sample_match_all[, 5],
                     tissue = sample_match_all[, 7])
all_samples = all_samples %>% mutate(data = map(file, load_count_file)) %>% unnest(data)
# add unmutated label
all_samples = all_samples %>% left_join(parent_data, by = "id") %>%
        mutate(spacer = map2(spacer, parent_spacer, function(s, p) left_join(s, p, by ="spacer", suffix = c("", "_parent")) %>%
                                     mutate(unmutated = !is.na(count_parent)) %>%
                                     select(-count_parent))) %>% select(-parent_spacer)

# pooling all spacers for each embryo
all_emb = all_samples %>% group_by(emb_name, id) %>% summarise(time = check_equal(time),
                                                               mother = check_equal(mother),
                                                               is_embyro = check_equal(is_embyro),
                                                               spacer = list(reduce(spacer, rbind) %>%
                                                                                     group_by(spacer) %>%
                                                                                     summarise(count = sum(count),
                                                                                               unmutated = check_equal(unmutated))))

# mutated fraction and observed diversity
all_emb = all_emb %>% mutate(mutated_fraction = map_dbl(spacer, function(x) {
        sum(x$count * (1 - x$unmutated))/sum(x$count)
        }),
        observed_diversity = map_int(spacer, function(x) {
                sum((x$count > (sum(x$count) * 0.01)) & !x$unmutated)
                }))
all_emb = subset(all_emb, is_embyro)
assertthat::assert_that(all(active_id %in% all_emb$id))

# naive mutation rate estimates
all_emb = all_emb %>% ungroup() %>% mutate(all_emb, mutrate_estimate = -log(1. - mutated_fraction)/(time - 0.6))
id_mutrate = all_emb %>% group_by(id) %>% summarise(mutrate_estimate = mean(pmin(pmax(10^(-3.2), mutrate_estimate), 10^0.5)))

left_join(id_mutrate, table_s3, by = c("id" = "identifier")) %>% group_by(Class) %>% summarize(mean(log10(mutrate_estimate)))

# function to simulate the MARC1 time course
sim_wrap <- function(sample_time, activation_time) {
        double_time0 = 0.6
        double_time1 = 0.35
        n_gen_stage = 7
        max_gen = 30
        sample_size = 2000
        frac_cut = 0.01
        mp0 = make_mut_param_by_rate(m_rate_seq,
                                     recur_prob = 0.,
                                     recur_vec = NA,
                                     active_time = list(c(activation_time, sample_time)))
        n_gen =  ifelse(sample_time < double_time0 * n_gen_stage,
                        floor(sample_time/double_time0),
                        n_gen_stage + pmin(max_gen - n_gen_stage,
                                           floor((sample_time - double_time0 * n_gen_stage)/double_time1)))
        n_gen0 = pmin(n_gen, n_gen_stage)
        dum_gen0 = list(cell_type = NA,
                        start_time = 0,
                        start_count = 1,
                        start_mode_counts = NA,
                        cell_counts = 2^(0:n_gen0),
                        num_gen = n_gen0,
                        double_time = double_time0,
                        end_time = double_time0 * n_gen0,
                        end_count = 2^n_gen0,
                        end_mode_counts = NA,
                        active = F)
        if (n_gen <= n_gen_stage) {
                dum_gen0 = generate_sample_size(dum_gen0, sample_size = pmin(dum_gen0$end_count, sample_size))
                dum_gen0$start_mut_counts = init_mut_counts_list(mp0, num_cells = 1)
                dum_gen0 = sample_mut_history_gens(dum_gen0,
                                                   mut_param = mp0)
                allele_mat_list = extract_allele_matrix(list(dum_gen0))
        } else {
                dum_gen1 = list(cell_type = NA,
                                start_time = dum_gen0$end_time,
                                start_count = dum_gen0$end_count,
                                start_mode_counts = NA,
                                num_gen = n_gen - n_gen_stage,
                                double_time = double_time1,
                                end_time = n_gen_stage * double_time0 + double_time1 * (n_gen - n_gen_stage),
                                end_count = dum_gen0$end_count * 2^(n_gen - n_gen_stage),
                                end_mode_counts = NA,
                                active = F)
                dum_gen1 = generate_sample_size(dum_gen1, sample_size = pmin(dum_gen1$end_count, sample_size))
                dum_gen0 = generate_sample_size(dum_gen0, sample_size = dum_gen1$sample_size[1])
                dum_gen0$start_mut_counts = init_mut_counts_list(mp0, num_cells = 1)
                dum_gen0 = sample_mut_history_gens(dum_gen0,
                                                   mut_param = mp0)
                dum_gen1$start_mut_counts = dum_gen0$end_mut_counts
                dum_gen1 = sample_mut_history_gens(dum_gen1,
                                                   mut_param = mp0,
                                                   target_time = sample_time)
                allele_mat_list = extract_allele_matrix(list(dum_gen1))
        }
        return(tibble(mutation_rate = m_rate_seq,
                      allele = allele_mat_list))
}
# sim_out = (1:100) %>% map_dfr(function(i) {
#         message(i)
#         tibble(expand.grid(sample_time = sort(unique(all_emb$time)),
#                            activation_time = c(0.6))) %>% mutate(data = pmap(., sim_wrap))
# })
# saveRDS(sim_out, file = "./intermediate_data/sim_MARC1_time_course.rds")
sim_out = readRDS("./intermediate_data/sim_MARC1_time_course.rds")
# sim_out = subset(sim_out, activation_time == 0.6)
sim_out = sim_out %>% unnest(data) %>% group_by(sample_time, activation_time, mutation_rate) %>% summarise(data_list = list(c(allele)))
sim_out = sim_out %>% mutate(mutated_fraction = map(data_list, function(x) x %>% map_dbl(get_mutated_frac)),
                             total_diversity = map(data_list, function(x) x %>% map_int(get_sim_diversity, frac = 0.005)))
sim_out = sim_out %>% mutate(mf_den = map(mutated_fraction, function(x) density(logit(x))),
                             div_den = map(total_diversity, function(x) {
                                     freq = table(x)/length(x)
                                     den = list(x = as.numeric(names(freq)), y = freq)
                             })) %>% ungroup()
all_emb_like = all_emb %>% left_join(sim_out %>% select(sample_time, activation_time, mutation_rate, mf_den) %>%
                                     nest(fit_data = -sample_time), by = c("time" =  "sample_time")) %>%
        mutate(fit_data = map2(mutated_fraction, fit_data, function(mf, s_data) {
                mutate(s_data, prob = map_dbl(mf_den, function(den) {
                        prob = approxfun(den)(logit(mf))
                        ifelse(is.na(prob), 0, prob)
                })) %>% select(-mf_den)
                }))
id_prob = all_emb_like %>% unnest(cols = fit_data) %>%
        group_by(id, mutation_rate) %>%
        summarise(log_prob = sum(log(1e-10+prob))) %>% mutate(scaled_log_prob = log_prob - max(log_prob))
# saveRDS(id_prob, "intermediate_data/MARC1_id_prob.rds")
id_prob_max = id_prob %>% top_n(n = 1, wt = scaled_log_prob) %>% arrange(mutation_rate) %>% left_join(table_s3, by = c("id" = "identifier"))
col_func = circlize::colorRamp2(c(-3.2, -1.35, 0.5), colors = c("green", "yellow", "red"))
# id_prob_max %>% mutate(color = map_chr(log10(mutation_rate), col_func))

# do high density interval estiamation
# These packages are not required by need to be installed to reproduce this part of the analysis
library(HDInterval)
library(sfsmisc)
id_prob_est = id_prob %>% nest(data = -id) %>% subset(id %in% active_id) %>% mutate(mr_smooth_den = map(data, function(mr_data) {
        prob_val = exp(mr_data$scaled_log_prob)
        b = ksmooth(log10(mr_data$mutation_rate), prob_val, bandwidth = 0.3)
        b$y = b$y/integrate.xy(b)
        b_den = rlang::duplicate(b)
        class(b_den)  = "density"
        b_den
}))
id_prob_est = id_prob_est %>% mutate(mr_estimates = map(mr_smooth_den, function(b_den) {
        interval = hdi(b_den, allowSplit = F, credMass = 0.9)
        tibble(mr_est = 10^(b_den$x[which.max(b_den$y)]),
               mr_lo90 =  10^(interval[1]),
               mr_hi90 = 10^(interval[2]))
}))
id_prob_est = left_join(id_prob_est, rename(id_mutrate, mr_naive_estimate = mutrate_estimate))
id_prob_est_sum = id_prob_est %>% unnest(cols = mr_estimates) %>% select(-c(data, mr_smooth_den))
write_tsv(id_prob_est_sum, path = "./supplementary_data/supplementary_data_2_mutation_rates/MARC1_all_active_mutation_rates.tsv")

# plotting individual density
mr_glist = purrr::map(id_prob_est$id, function(id_name) {
        den = id_prob_est$mr_smooth_den[[which(id_prob_est$id == id_name)]]
        est = id_prob_est$mr_estimates[[which(id_prob_est$id == id_name)]]
        class(den) = "list"
        ggline(data.frame(den), x = "x", y = "y", numeric.x.axis = T, plot_type = "l", xlab = "log10(MR)", ylab = "Density") +
                geom_vline(xintercept = log10(est$mr_lo90), linetype = 2) + geom_vline(xintercept = log10(est$mr_hi90), linetype = 2) +
                # adding naive estimates
                geom_vline(xintercept = log10(id_prob_est$mr_naive_estimate[which(id_prob_est$id == id_name)]), color = "red")
})
mr_glist[[25]]
id_den_plot = id_prob_est %>% transmute(data = map(mr_smooth_den, function(x) {
        class(x) = "list"
        as_tibble(x)
})) %>% unnest(cols = data) %>% left_join(table_s3, by = c("id" = "identifier"))
ggplot(id_den_plot, aes(x, y, color = Class, group = id)) + geom_line() + theme_bw()
# ggplot(id_den_plot %>% filter(id == "AAGCCGCGCG"), aes(x, y, color = Class)) + geom_line() + theme_bw()

id_mr_prob = t(as.matrix(spread(select(id_prob, -log_prob), key = id, value = scaled_log_prob)[, -1]))
id_mr_prob = t(scale(t(id_mr_prob)))
# id_ord = split(1:length(active_id), table_s3$Class[match(active_id, table_s3$identifier)]) %>%
#         map(function(idx) rev(idx[hclust(dist(id_mr_prob[active_id[idx], ]))$order])) %>% flatten_int()
id_ord = split(1:length(active_id), table_s3$Class[match(active_id, table_s3$identifier)]) %>%
        map(function(idx) rev(idx[order(apply(id_mr_prob[active_id[idx], ], 1, which.max))])) %>% flatten_int()

# Figure S3A
h_mr = Heatmap(rbind(log10(m_rate_seq)), cluster_columns = F, col = col_func)
h_est = Heatmap(id_mr_prob[active_id, ], cluster_columns = F, row_order = id_ord, heatmap_legend_param = list(direction = "horizontal"),
        left_annotation = rowAnnotation(Class = table_s3$Class[match(active_id, table_s3$identifier)],
                                        Founder = table_s3$Founder[match(active_id, table_s3$identifier)],
                                        col = list(Class = c("fast" = "#fc8d62",
                                                             "mid" = "#66c2a5",
                                                             "slow" = "#8da0cb"),
                                                   Founder = c("PB3" = "#2f4b7c",
                                                               "PB7" = "#a05195"))),
        show_row_names = F, name = "Scaled_log_prob")
h_mr %v% h_est


# id_mr_cdf = id_prob %>% nest(mr_data = -id) %>% mutate(mr_cdf = map(mr_data, function(x) {
#         b = ksmooth(log10(x$mutation_rate), exp(x$log_prob - max(x$log_prob)), bandwidth = 0.6)
#         b$y = b$y/sfsmisc::integrate.xy(b)
#         pdf_to_cdf(b)
# }))
expit <- function(x){
        exp(x)/(1+exp(x))
}
num_sim = 100
id_simulated = tibble(id = active_id) %>%
        left_join(id_prob %>% nest(mr_data = -id), by = "id") %>%
        mutate(simulated = map(mr_data, function(x) tibble(time = rep(unique(all_emb$time), num_sim),
                                                           mutation_rate =  sample(x$mutation_rate,
                                                                                   size = length(unique(all_emb$time)) * num_sim,
                                                                                   prob = exp(x$scaled_log_prob),
                                                                                   replace = T))))
id_simulated_data = id_simulated %>%
        select(-mr_data) %>%
        unnest(cols = simulated) %>%
        group_by(time, mutation_rate) %>% nest(data = id) %>% ungroup()
id_simulated_data = id_simulated_data %>%
        left_join(select(subset(sim_out, activation_time == 0.6), sample_time, mutation_rate, mf_den, div_den),
                  by = c("time" = "sample_time", "mutation_rate")) %>%
        mutate(data = map2(data, mf_den, function(data, den) {
                mutate(data, mutated_fraction = expit(sample(den$x, size = nrow(data), prob = den$y, replace = T)))
                })) %>%
        mutate(data = map2(data, div_den, function(data, den) {
                mutate(data, observed_diversity = sample(den$x, size = nrow(data), prob = den$y, replace = T))
        })) %>%
        select(time, mutation_rate, data) %>% unnest(cols = data) %>% group_by(id, time)
combined_data = rbind(id_simulated_data %>% select(id, time, mutated_fraction, observed_diversity) %>% mutate(type = "simulated") %>% ungroup(),
                      subset(all_emb, id %in% active_id) %>% ungroup() %>% select(id, time, mutated_fraction, observed_diversity) %>% mutate(type = "observed"))
# combined_data = left_join(combined_data, id_prob_max, by = "id")
library(patchwork)
pdf("./plots/MARC1/mut_rate.pdf", onefile = T)
combined_data %>% group_by(id) %>% group_walk(
        ~ print(ggline(.x, x = "time", y = "mutated_fraction", color = "type", add = c("mean_se"), numeric.x.axis = T, ylim = c(0, 1)) +
                ggboxplot(.x, x = "time", y = "observed_diversity", color = "type", add = "density", numeric.x.axis = T, ylim = c(0, 30)) +
                        plot_annotation(title = paste0(.y, ", class: ", table_s3$Class[match(.y, table_s3$identifier)]))
                )
        )
# dev.off()

# g_list0 = c("fast", "mid", "slow") %>% map(function(x) {
#         combined_data %>% subset(id %in% table_s3$identifier[table_s3$Class == x]) %>%
#                 ggline(x = "time", y = "mutated_fraction", color = "id", facet.by = "type", add = "mean_se",
#                        xlab = "Emb Time", ylab = "Mutated Fraction",
#                        # palette = (group_by(., id) %>% summarize(color = first(color)))$color,
#                        ylim = c(0, 1), numeric.x.axis = T) +
#                 theme(legend.position = "none",
#                       strip.background = element_blank(),
#                       text = element_text(size = 12))
# })
# g_list1 = c("fast", "mid", "slow") %>% map(function(x) {
#         combined_data %>% subset(id %in% table_s3$identifier[table_s3$Class == x]) %>%
#                 ggboxplot(x = "time", y = "observed_diversity", color = "id", facet.by = "type",
#                           xlab = "Emb Time", ylab = "Allele Diversity",
#                           # palette = (group_by(., id) %>% summarize(color = first(color)))$color,
#                           numeric.x.axis = T, ylim = c(0, 30)) +
#                 theme(legend.position = "none",
#                       strip.background = element_blank(),
#                       text = element_text(size = 12))
# })
# g_simulation_all = (g_list0[[1]] + g_list0[[2]] + g_list0[[3]]) /
#         (g_list1[[1]] + g_list1[[2]] + g_list1[[3]])
# g_simulation_all
# push_pdf(g_simulation_all, file_name = "mutation_sim_vs_obs", width = 7.0, height = 4.0)

combined_data %>% left_join(select(table_s3, identifier, Class), by = c(id = "identifier")) %>%
        ggline(x = "time", y = "mutated_fraction", color = "id", facet.by = c("Class", "type"), add = "mean_se",
               xlab = "Emb Time", ylab = "Mutated Fraction",
               # palette = (group_by(., id) %>% summarize(color = first(color)))$color,
               ylim = c(0, 1), numeric.x.axis = T) +
        theme(legend.position = "none",
              strip.background = element_blank(),
              text = element_text(size = 12))

combined_data %>% left_join(select(table_s3, identifier, Class), by = c(id = "identifier")) %>%
        ggboxplot(x = "time", y = "observed_diversity", fill = "id", facet.by = c("Class", "type"), color = "id",
                  xlab = "Emb Time", ylab = "Allele Diversity",
                  numeric.x.axis = T, ylim = c(0, 30)) +
        theme(legend.position = "none",
              strip.background = element_blank(),
              text = element_text(size = 12))

(combined_data %>% subset(id %in% table_s3$identifier[table_s3$Class == "mid"]) %>%
        ggboxplot(x = "time", y = "observed_diversity", color = "id", facet.by = "type", numeric.x.axis = T, ylim = c(0, 30)) + theme(legend.position = "none")) /
(combined_data %>% subset(id %in% table_s3$identifier[table_s3$Class == "slow"]) %>%
        ggboxplot(x = "time", y = "observed_diversity", color = "id", facet.by = "type", numeric.x.axis = T, ylim = c(0, 30)) + theme(legend.position = "none"))

indelphi_tb = readRDS("./intermediate_data/MARC1_indelphi_predicted.rds")
all_emb = all_emb %>% nest(data = -id) %>% subset(id %in% active_id) %>% left_join(indelphi_tb)

# source("./analysis/MARC1_estimation_func.R")
# summarize data for each id
# all_emb_id = all_emb %>% group_by(id) %>% nest(data = -id)
# all_emb_id = all_emb_id %>% mutate(spacer_ind = map(data, function(x) {
#         spacer_counts = spread(select(subset(x %>% unnest(spacer), unmutated == F),
#                                       emb_name, spacer, count),
#                                key = emb_name, value = count)
#         out = t(as.matrix(spacer_counts[, -1]))
#         colnames(out) = spacer_counts$spacer
#         out = (!is.na(out)) * 1
#         out
# }))
# all_emb_id = all_emb_id %>% mutate(p_vec = map(spacer_ind, function(y) {
#         n_vec = rowSums(y)
#         if (ncol(y) > 0) {
#                 p_mar = sapply(1:ncol(y), function(j) {
#                         message(j)
#                         optim_out = optimParallel(par = 0.5, lower = 0.0+1e-3, upper = 1.0-1e-3,
#                                                   fn = margin_func,
#                                                   n_vec = n_vec,
#                                                   k_vec = y[, j])
#                         return(optim_out$par)
#                 })
#                 p_vec = p_mar / sum(p_mar)
#                 return(p_vec)
#         } else {
#                 return(NULL)
#         }
# }))
# # reorder based on chance of occurring
# all_emb_id = all_emb_id %>% mutate(temp_spacer_ord = map(p_vec, function(p_vec) {
#         if (length(p_vec) >= 1) {
#                 return(order(p_vec, decreasing = T))
#         } else {
#                 return(NULL)
#         }
# }))
# all_emb_id = all_emb_id %>%
#         mutate(p_vec = map2(p_vec, temp_spacer_ord,
#                             function(x, ord) {
#                                     if (length(x) >= 1) {
#                                             return(x[ord])
#                                             } else {
#                                                     return(NULL)
#                                             }
#                                     }
#                             ),
#                spacer_ind = map2(spacer_ind, temp_spacer_ord,
#                                  function(spacer_ind ,ord) {
#                                          if (ncol(spacer_ind) >= 1) {
#                                                  return(spacer_ind[, ord, drop = F])
#                                          } else {
#                                                  return(spacer_ind)
#                                          }
#                                  }))
# saveRDS(all_emb_id, file = "./intermediate_data/MARC1_hgRNA_estimates.rds")

all_emb_id = readRDS("./intermediate_data/MARC1_hgRNA_estimates.rds")

all_emb_id = all_emb_id %>% subset(id %in% active_id) %>% left_join(indelphi_tb)
all_emb_id = all_emb_id %>% left_join(transmute(table_s3, id = identifier, class = Class), by = "id")
all_emb_id = all_emb_id %>% mutate(spacer_pred = map2(all_spacer, indelphi, function(x, y) {
        y = mutate(y, obs = str_sub(obs, 8))
        out = full_join(x, y, by = c("spacer" = "obs"))
        out = mutate(out, log_prob_norm = log(prob / sum(prob, na.rm = T)))
        out
}))

all_emb_id = all_emb_id %>% subset(id %in% active_id)
p_vec_df = map2(all_emb_id$id, all_emb_id$p_vec, function(x, y) {
        if (!is.null(y)) {
                return(
                        tibble(id = x, prob = log(y), rank = seq_along(y))
                )
        }
}) %>% reduce(bind_rows)
p_vec_df = p_vec_df %>% left_join(table_s3, by = c("id" = "identifier")) %>% select(id, prob, rank, Class)
p_vec_dir_df = (1:20) %>% map(function(i) {
        p_vec_dir = sort(extraDistr::rdirichlet(n = 1, alpha = rep(0.1, 1000))[1, ], decreasing = T)
        tibble(id = paste0("dir_", i),
               prob = log(p_vec_dir),
               rank = seq_along(p_vec_dir),
               Class = "Dir")
}) %>% reduce(bind_rows)
p_vec_df = rbind(p_vec_df, p_vec_dir_df)
# ggplot(p_vec_df, aes(x = rank, y = prob, color = id)) + geom_line() + xlim(c(0, 500)) + facet_wrap(~Class) + ylim(c(-20, 0)) + theme_bw()+ theme(legend.position = "none")
ggplot(p_vec_df, aes(x = rank, y = prob, group = id, color = Class)) + geom_line() + xlim(c(0, 500)) + ylim(c(-20, 0)) + theme_bw()

# below are simulations of emb data using indelphi probabilities and posterior mutation rates
id_mrden = id_prob %>% nest(mrden = -id)
id_mrden = id_mrden %>% filter(id %in% active_id)
id_mrden = left_join(id_mrden, indelphi_tb) %>% mutate(r_vec = map(indelphi, function(x) {
        r_vec = x$prob
        r_vec = sort(r_vec/sum(r_vec), decreasing = T)
        r_vec
})) %>% select(-indelphi)
sample_mrden <- function(mrden, n) {
        sample(mrden$mutation_rate,
               size = n,
               prob = exp(mrden$scaled_log_prob),
               replace = T)
}

# mut_p = readRDS("../LTModelData/temp_mut_p.rds")
# make a scatter plot displaying in division rate in mouse embyro
# Figure S1
double_time0 = 0.6
double_time1 = 0.35
n_gen_stage = 7
max_gen = 30
n_gen = 30
n_gen0 = pmin(n_gen, n_gen_stage)
time_seq = c(double_time0 * (0:n_gen_stage))
stage0_end_time = time_seq[length(time_seq)]
time_seq = c(time_seq, stage0_end_time + double_time1 * (1:(max_gen - n_gen_stage)))
count_seq = c(2^(0:n_gen0))
stage0_end_count = count_seq[length(count_seq)]
count_seq = c(count_seq, stage0_end_count * 2^(1:(max_gen - n_gen_stage)))
sim_tb = tibble(time = time_seq, log2count = log2(count_seq))
obs_tb = tibble(time = c(0, 2., 2.5, 3.5, 4.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0),
                log2count = log2(c(c(1, 4, 12, 32, 64), c(120, 250, 660, 4500, 15000, 75000) * 2)))
sim_tb %>%
        ggline(x = "time", y = "log2count",
               plot_type  = "l",
               numeric.x.axis = T,
               xlim = c(0, 10),
               xlab = "Time (days)",
               ylab = "log2(Cell count)") +
        geom_point(data = obs_tb, aes(x = time, y = log2count))


simulate_data_pvec <- function(mr, r_vec, sample_time, activation_time = 0.6) {
        # dummy param
        mr = mut_p$mut_rate[1]
        r_vec = mut_p$recur_vec_list[[1]]
        sample_time = 16.5
        activation_time = 0.6
        # end dummy param

        # built-in MARC1 parameters
        double_time0 = 0.6
        double_time1 = 0.35
        n_gen_stage = 7
        max_gen = 30
        sample_size = 2000
        frac_cut = 0.01
        # end built-in MARC1 parameters

        mp0 = make_mut_param_by_rate_rvec(mr,
                                          recur_prob_vec = 1.,
                                          recur_vec_list = list(r_vec),
                                          active_time = list(c(activation_time, sample_time)))

        n_gen =  ifelse(sample_time < double_time0 * n_gen_stage,
                        floor(sample_time/double_time0),
                        n_gen_stage + pmin(max_gen - n_gen_stage,
                                           floor((sample_time - double_time0 * n_gen_stage)/double_time1)))
        n_gen0 = pmin(n_gen, n_gen_stage)
        dum_gen0 = list(cell_type = NA,
                        start_time = 0,
                        start_count = 1,
                        start_mode_counts = NA,
                        cell_counts = 2^(0:n_gen0),
                        num_gen = n_gen0,
                        double_time = double_time0,
                        end_time = double_time0 * n_gen0,
                        end_count = 2^n_gen0,
                        end_mode_counts = NA,
                        active = F)
        if (n_gen <= n_gen_stage) {
                dum_gen0 = generate_sample_size(dum_gen0, sample_size = pmin(dum_gen0$end_count, sample_size))
                dum_gen0$start_mut_counts = init_mut_counts_list(mp0, num_cells = 1)
                dum_gen0 = sample_mut_history_gens(dum_gen0,
                                                   mut_param = mp0)
                allele_mat_list = extract_allele_matrix(list(dum_gen0))
        } else {
                dum_gen1 = list(cell_type = NA,
                                start_time = dum_gen0$end_time,
                                start_count = dum_gen0$end_count,
                                start_mode_counts = NA,
                                cell_counts = dum_gen0$end_count * 2^(n_gen - n_gen_stage),
                                num_gen = n_gen - n_gen_stage,
                                double_time = double_time1,
                                end_time = n_gen_stage * double_time0 + double_time1 * (n_gen - n_gen_stage),
                                end_count = dum_gen0$end_count * 2^(n_gen - n_gen_stage),
                                end_mode_counts = NA,
                                active = F)
                dum_gen1 = generate_sample_size(dum_gen1, sample_size = pmin(dum_gen1$end_count, sample_size))
                dum_gen0 = generate_sample_size(dum_gen0, sample_size = dum_gen1$sample_size[1])
                dum_gen0$start_mut_counts = init_mut_counts_list(mp0, num_cells = 1)
                dum_gen0 = sample_mut_history_gens(dum_gen0,
                                                   mut_param = mp0)
                dum_gen1$start_mut_counts = dum_gen0$end_mut_counts
                dum_gen1 = sample_mut_history_gens(dum_gen1,
                                                   mut_param = mp0,
                                                   target_time = sample_time)
                allele_mat_list = extract_allele_matrix(list(dum_gen1))
        }
        return(allele_mat_list)
}

sim_pvec = as_tibble(expand.grid(time = sort(unique(all_emb$time)),
                                 id = id_mrden$id))
# library(furrr)
# plan(multisession, workers = 12)
# num_sim = 10
# sim_pvec = sim_pvec %>% left_join(id_mrden) %>% mutate(simulated_data = future_pmap(., function(mrden, r_vec, time, ...) {
#         map(1:num_sim, function(i) {
#                 simulate_data_pvec(sample_mrden(mrden, 1), r_vec, sample_time = time)
#         })
# }, .progress = T))
# saveRDS(sim_pvec, file = "./intermediate_data/sim_MARC1_mutation.rds")

sim_pvec = readRDS("./intermediate_data/sim_MARC1_mutation.rds")
sim_pvec_res = select(sim_pvec, time, id, simulated_data) %>% unnest(simulated_data)
sim_pvec_res = mutate(sim_pvec_res, simulated_data_proc = purrr::map(simulated_data, function(x) {
        filter_allele_matrix(collapse_reucr_ver(x[[1]]), frac_cut = 0.01)
}))

col_func = circlize::colorRamp2(breaks = c(-3.2, -2.2, -1, 0, 0.5), colors = c( "#325288", "#325288", "#ffdf6b", "#a20a0a", "#a20a0a"))
# illustrate the bad performance of naive allele praobility estimates
dummy_y_loc = -10
spacer_prob_example = filter(sim_pvec_res, id == "AAGCCGCGCG") %>% transmute(id, allele_prob_est = map(simulated_data_proc, function(x) {
        cts = x[1, colnames(x)!=0]
        tibble(spacer = names(cts), prob = cts / sum(cts))
})) %>% unnest(cols = allele_prob_est)
spacer_prob_example = spacer_prob_example %>% group_by(spacer) %>% summarize(prob = mean(prob)) %>% mutate(log_prob_norm = log(prob/sum(prob)))
r_vec = sim_pvec$r_vec[[which(sim_pvec$id ==  "AAGCCGCGCG")[1]]]
spacer_prob_example = full_join(spacer_prob_example, tibble(spacer = paste0("R-", 1:length(r_vec)), true_log_prob = log(r_vec)))
spacer_prob_example$is_obs = ifelse(!is.na(spacer_prob_example$log_prob_norm), "Observed", "Unobserved")
spacer_prob_example$log_prob_norm[is.na(spacer_prob_example$log_prob_norm)] = rnorm(sum(is.na(spacer_prob_example$log_prob_norm)), mean = dummy_y_loc, sd = 0.1)
# g1 = ggscatter(spacer_prob_example, x = "true_log_prob", y = "log_prob_norm", color = "is_obs", ylim = c(dummy_y_loc - 0.5, 0)) + geom_abline(intercept = 0, slope = 1, color = "red") + xlab("True Allele Probability") + ylab("Estimated Allele Probability")
# allele probability estimates using indepdent samples
spacer_prob_example1 = filter(sim_pvec_res, id == "AAGCCGCGCG") %>% transmute(id, allele_occurance = map(simulated_data_proc, function(x) {
        cts = x[1, colnames(x)!=0]
        if (length(cts) == 0) {
                return(NULL)
        }
        tibble(spacer = names(cts), occurance = 1.)
})) %>% unnest(cols = allele_occurance)

spacer_prob_example1 = spacer_prob_example1 %>% group_by(spacer) %>% summarise(prob = sum(occurance)/num_sim)
spacer_prob_example1 = mutate(spacer_prob_example1, log_prob_norm = log(prob/sum(prob)))
spacer_prob_example1 = full_join(spacer_prob_example1, tibble(spacer = paste0("R-", 1:length(r_vec)), true_log_prob = log(r_vec)))
spacer_prob_example1$is_obs = ifelse(!is.na(spacer_prob_example1$log_prob_norm), "Observed", "Unobserved")
spacer_prob_example1$log_prob_norm[is.na(spacer_prob_example1$log_prob_norm)] = rnorm(sum(is.na(spacer_prob_example1$log_prob_norm)), mean = dummy_y_loc, sd = 0.1)

# final figure for simulated allele prob estimates
bind_rows(mutate(spacer_prob_example, Method = "Within Sample Estimate"),
          mutate(spacer_prob_example1, Method = "Across Sample Estimate")) %>%
        ggscatter(x = "true_log_prob", y = "log_prob_norm", color = "is_obs",
                  ylim = c(dummy_y_loc - 0.5, 0),
                  xlab = "log(True Allele Probability)", ylab = "log(Estimated Allele Probability)") %>%
        facet(facet.by = "Method") + geom_abline(intercept = 0, slope = 1, color = "#afb9c8", size = 1) + scale_color_manual(values=  sim_obs_col)

sim_pvec_res = mutate(sim_pvec_res, observed_diversity = map_dbl(simulated_data_proc, function(x) {
        sum(colnames(x) != "0")
}))
sim_pvec_res = mutate(sim_pvec_res, mutated_fraction = map_dbl(simulated_data_proc, function(x) {
        sum(x[colnames(x) != "0"])/sum(x)
}))

combined_data = bind_rows(select(sim_pvec_res, id, time, mutated_fraction , observed_diversity) %>% mutate(type = "Simulated"),
                          filter(all_emb, id %in% active_id) %>% ungroup() %>% select(id, time, mutated_fraction, observed_diversity) %>% mutate(type = "Observed"))

# testing color
id_col = col_func(log10(id_prob_max$mutation_rate))
names(id_col) = id_prob_max$id
# names(id_col) = c(table_s3$identifier[table_s3$Class == "slow"],
#                   table_s3$identifier[table_s3$Class == "mid"],
#                   table_s3$identifier[table_s3$Class == "fast"])

# col_slow = col_func(seq(-3, -1, length = 32)); names(col_slow) = table_s3$identifier[table_s3$Class == "slow"]
# col_mid = col_func(seq(-1.1, 1, length = 23)); names(col_mid) = table_s3$identifier[table_s3$Class == "mid"]
# col_fast = col_func(seq(1.2, 3, length = 9)); names(col_fast) = table_s3$identifier[table_s3$Class == "fast"]

# col_slow = rand_color(32, hue = "red"); names(col_slow) = table_s3$identifier[table_s3$Class == "slow"]
# col_mid = rand_color(23, hue = "green"); names(col_mid) = table_s3$identifier[table_s3$Class == "mid"]
# col_fast = rand_color(9, hue = "blue"); names(col_fast) = table_s3$identifier[table_s3$Class == "fast"]
# id_col = c(col_slow, col_mid, col_fast)
# assertthat::assert_that(!anyDuplicated(c(col_fast, col_mid, col_slow)))
# plot(NULL, xlim = c(1, 32), ylim = c(1, 3), axes = FALSE, ann = FALSE)
# points(1:32, rep(1, 32), pch = 16, cex = 5,
#        col = col_slow)
# points(1:23, rep(2, 23), pch = 16, cex = 5,
#        col = col_mid)
# points(1:9, rep(3, 9), pch = 16, cex = 5,
#        col = col_fast)
# table_s3$color = id_col[match(table_s3$identifier, names(id_col))]
# end testing color

#  %>% left_join(select(table_s3, identifier, Class), by = c(id = "identifier"))

# g_mr = combined_data %>% group_by(id, time, type) %>% summarise(mean = mean(mutated_fraction), se = sqrt(mutated_fraction * (1 - mutated_fraction) / n())) %>%
#         left_join(select(table_s3, identifier, Class), by = c(id = "identifier")) %>%
#         ggplot(aes(x = time, y = mean, color = id)) + geom_line() + facet_grid(Class ~ type) + geom_errorbar(aes(ymin = pmax(0., mean - se), ymax = pmin(1., mean + se)), width = 0.1) +
#         scale_color_manual(values = id_col) + scale_y_continuous(limits = c(0, 1))+ xlab("Emb Time") + ylab("Mutated Fraction")+
#         theme_pubr() + theme(legend.position = "none",
#               strip.background = element_blank(),
#               text = element_text(size = 11))

g_mr = combined_data %>% left_join(select(table_s3, identifier, Class), by = c(id = "identifier")) %>%
        ggline(x = "time", y = "mutated_fraction", color = "id", facet.by = c("Class", "type"), add = "mean_se", size = 0.5, width = 0.5,
               xlab = "Emb Time", ylab = "Mutated Fraction", numeric.x.axis = T) + scale_color_manual(values = id_col) + scale_y_continuous(limits = c(0, 1))+
        theme(legend.position = "none",
              strip.background = element_blank(),
              text = element_text(size = 11))
g_div = combined_data %>% left_join(select(table_s3, identifier, Class), by = c(id = "identifier")) %>%
        ggboxplot(x = "time", y = "observed_diversity", fill = "id", facet.by = c("Class", "type"), color = "id",
                  xlab = "Emb Time", ylab = "Mutant Allele Diversity", outlier.size = 0.1, outlier.stroke = 0.1, width = 0.2, lwd = 0.1,
                  numeric.x.axis = T, ylim = c(0, 30)) + scale_fill_manual(values = id_col) + scale_color_manual(values = id_col) +
        theme(legend.position = "none",
              strip.background = element_blank(),
              text = element_text(size = 11))

g_mr$layers[[3]]$aes_params$size = 0.25
g_all = g_mr + g_div
g_all
push_pdf(g_all, "MARC1_simulation_compare", width = 7.2, height= 4, dir = "./plots/MARC1/")

g_mr_div  = combined_data %>% left_join(select(id_prob_max, id, mutation_rate)) %>% group_by(type, id) %>% summarise(log_mr = log10(mutation_rate), mean_obs = mean(observed_diversity)) %>%
        ggline(x = "log_mr", y = "mean_obs", add = "mean_sd", color = "type", numeric.x.axis = T, ylab = "Average Diversity", xlab = "log10(Mutation Rate)") +
        scale_color_manual(values = sim_obs_col) + theme(legend.position = "right", text = element_text(size = 11), legend.text = element_text(size = 10))
push_pdf(g_mr_div, "MARC1_mutation_rate_diversity", width = 4, height= 2.8, dir = "./plots/MARC1/")

# g_list0 = c("fast", "mid", "slow") %>% map(function(x) {
#         combined_data %>% subset(id %in% table_s3$identifier[table_s3$Class == x]) %>%
#                 ggline(x = "time", y = "mutated_fraction", color = "id", facet.by = "type", add = "mean_sd",
#                        # palette = (group_by(., id) %>% summarize(color = first(color)))$color,
#                        ylim = c(0, 1), numeric.x.axis = T) +
#                 theme(legend.position = "none", strip.background = element_blank())
# })
# g_list1 = c("fast", "mid", "slow") %>% map(function(x) {
#         combined_data %>% subset(id %in% table_s3$identifier[table_s3$Class == x]) %>%
#                 ggboxplot(x = "time", y = "observed_diversity", color = "id", facet.by = "type",
#                           # palette = (group_by(., id) %>% summarize(color = first(color)))$color,
#                           numeric.x.axis = T, ylim = c(0, 30)) +
#                 theme(legend.position = "none", strip.background = element_blank())
# })
# (g_list0[[1]] + g_list0[[2]] + g_list0[[3]]) /
#         (g_list1[[1]] + g_list1[[2]] + g_list1[[3]])

sim_pvec_res = mutate(sim_pvec_res, allele_frequency = map(simulated_data_proc, function(r_count) {
        r_count = sort(r_count[1, colnames(r_count) != " 0"], decreasing = T)
        mutate(tibble(allele_rank = 1:length(r_count), count = r_count), allele_fraction = count/sum(count))
}))
sim_pvec_allele = sim_pvec_res %>% unnest(cols = allele_frequency)
all_emb = all_emb %>% mutate(allele_frequency = map(spacer, function(x) {
        if (any(!x$unmutated & (x$count > (sum(x$count) * 0.01)))) {
                r_count = x$count
                names(r_count) = x$spacer
                r_count = r_count[(r_count > (sum(r_count) * 0.01)) & !x$unmutated]
                r_count = sort(r_count, decreasing = T)
                return(
                        mutate(tibble(spacer = names(r_count),
                                      allele_rank = 1:length(r_count),
                                      count = r_count,
                                      total = sum(count)),
                               allele_fraction = (count)/sum(count))
                )
        } else {
                return(NULL)
        }
}))
# demonstrate bad allele probability estimates using averaged per sample proportions
all_emb_allele = select(all_emb, id, time, allele_frequency) %>% filter(id %in% active_id) %>% unnest(cols = allele_frequency)
all_emb_allele_prob = all_emb_allele %>% group_by(id, spacer) %>% summarise(allele_prob = mean(allele_fraction))
all_emb_allele_prob = all_emb_allele_prob %>% group_by(id) %>% mutate(allele_prob_norm = allele_prob / sum(allele_prob))

dummy_y_loc = -8
cominbed_allele = full_join(all_emb_allele_prob,
                            indelphi_tb %>% unnest(indelphi) %>% transmute(id, pred, spacer = str_sub(obs, 8), gen, prob, log_prob))
# normalize prediction so that observed sum to one
cominbed_allele = mutate(cominbed_allele, prob_norm = prob / sum(prob, na.rm = T))
cominbed_allele = cominbed_allele %>% mutate(log_prob_obs = log(allele_prob), log_prob_pred = log(prob_norm))
cominbed_allele$is_obs = ifelse(!is.na(cominbed_allele$log_prob_obs), "Observed", "Unobserved")
cominbed_allele$log_prob_obs_fill = cominbed_allele$log_prob_obs
cominbed_allele$log_prob_obs_fill[is.na(cominbed_allele$log_prob_obs_fill)] = rnorm(n = sum(is.na(cominbed_allele$log_prob_obs_fill)), mean = dummy_y_loc, sd = 0.1)
spacer_prob_example3 = filter(cominbed_allele, id == "AAGCCGCGCG") %>% transmute(log_p_vec = log_prob_obs_fill, log_prob_norm = log_prob_pred, is_obs)
spacer_prob_example2 = filter(all_emb_id, id == "AAGCCGCGCG")$spacer_pred[[1]] %>% mutate(log_p_vec = log(p_vec)) %>%
        mutate(is_pred = !is.na(gen), is_obs = ifelse(!is.na(p_vec), "Observed", "Unobserved"))
spacer_prob_example2$log_p_vec[is.na(spacer_prob_example2$log_p_vec)] = rnorm(n = sum(is.na(spacer_prob_example2$log_p_vec)), mean = dummy_y_loc, sd = 0.1)
spacer_prob_example2 = select(spacer_prob_example2, log_prob_norm, log_p_vec, is_obs)

# final figure for comparing estimated with indelphi prediciton using MARC1 data
g_x = bind_rows(mutate(spacer_prob_example3, Method = "Within Sample Estimate"),
          mutate(spacer_prob_example2, Method = "Across Sample Estimate")) %>%
        ggscatter(x = "log_prob_norm", y = "log_p_vec", color = "is_obs",
                  # ylim = c(dummy_y_loc - 0.5, 0),
                  ylim = c(- 8.5, 0),
                  xlab = "log(True Allele Probability)", ylab = "log(Estimated Allele Probability)") %>%
        facet(facet.by = "Method") + geom_abline(intercept = 0, slope = 1, color = "#afb9c8", size = 1) + scale_color_manual(values = sim_obs_col) +
        theme(legend.position = "right",
              strip.background = element_blank(),
              text = element_text(size = 12))
g_y = bind_rows(mutate(spacer_prob_example, Method = "Within Sample Estimate"),
          mutate(spacer_prob_example1, Method = "Across Sample Estimate")) %>%
        ggscatter(x = "true_log_prob", y = "log_prob_norm", color = "is_obs",
                  # ylim = c(dummy_y_loc - 0.5, 0),
                  ylim = c(- 10.5, 0),
                  xlab = "log(Indelphi Predicted Allele Probability)", ylab = "log(Estimated Allele Probability)") %>%
        facet(facet.by = "Method") + geom_abline(intercept = 0, slope = 1, color = "#afb9c8", size = 1) + scale_color_manual(values=  sim_obs_col) +
        theme(legend.position = "right",
              strip.background = element_blank(),
              text = element_text(size = 12))
push_pdf(g_x / g_y, "allele_prob_acros_within_sample", width = 7.2, height = 5., dir = "./plots/MARC1/")

# Figure for allele fraction
combined_allele_data = bind_rows(select(sim_pvec_allele, id, time, allele_rank, allele_fraction) %>% mutate(type = "Simulated"),
                                 filter(all_emb, id %in% active_id) %>% ungroup() %>%
                                         select(id, time,  allele_frequency) %>% unnest(cols = allele_frequency) %>% select(-count)  %>% mutate(type = "Observed"))
g_allele_frac = ggline(combined_allele_data, y = "allele_fraction", x = "allele_rank", add = "mean_sd", color = "type", ylab = "Mutant Allele Fraction",
                       xlab = "Mutant Allele Abundance Rank", xlim = c(1, 16), numeric.x.axis = T) +
        theme(legend.position = "none", text = element_text(size = 11), legend.text = element_text(size = 10)) +
        scale_color_manual(values = sim_obs_col)
push_pdf(g_allele_frac, "MARC1_allele_fraction", width = 3.5, height= 2.8, dir = "./plots/MARC1/")





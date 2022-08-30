output_dir = "./intermediate_data/panel_mod2_v1/"
plot_dir = "./plots/panel_mod2_v1/"
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))

# simulating the experiments
library(furrr)
plan(multisession, workers = 6)
exp_params = expand.grid(big_graph_id = 1:nrow(tree_panel),
                         sample_size = c(100),
                         num_element = c(50),
                         sampling = c("fixed", "proportional"),
                         i_sim = 1) %>%
        as_tibble()
# exp_params = mutate(exp_params, mut_p = map(num_element, function(i) {
#         sample_mutp_mid_or_fast(i, c(0.6, 11.5))
# }))
# saveRDS(exp_params, file = paste0(output_dir, "exp_data_1rep.rds"))

exp_params = readRDS(paste0(output_dir, "exp_data_1rep.rds"))

rep_dir = paste0(output_dir, "exp_data_1rep/")
dir.create(rep_dir)
# future_walk(1:nrow(exp_params), function(i) {
# future_walk(1:nrow(exp_params), function(i) {
future_walk(1:nrow(exp_params), function(i) {
        big_graph_id = exp_params$big_graph_id[i]
        sample_size = exp_params$sample_size[i]
        mut_p = exp_params$mut_p[[i]]
        sampling = exp_params$sampling[i]

        if (sampling == "fixed") {
                ss = sample_size
        }
        if (sampling == "proportional") {
                fm = tree_panel$type_graph[[big_graph_id]]
                gens0 = make_gens_mod2(fm)
                tip_size = map_dbl(gens0[fm$tip_id], "end_count")
                tip_sample = extraDistr::rmvhyper(nn = 1, k = length(fm$tip_id) * sample_size, n = tip_size)[1, ]
                tip_sample = pmax(tip_sample, 5)
                names(tip_sample) = fm$tip_id
                ss = tip_sample
        }
        out_file = paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        if (!file.exists(out_file)) {
                try({
                        out = simulate_sc_data_mod2(tree_panel$type_graph[[big_graph_id]], mut_p, sample_size = ss)
                        saveRDS(out, file = out_file)
                })
        }
}, .progress = T, .options = furrr_options(seed = T))

# exp_params = exp_params[exp_params$sampling == "fixed", ]
exp_params$data = map(1:nrow(exp_params), function(i) {
        out_file = paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        if (file.exists(out_file)) {
                out = readRDS(out_file)
        } else {
                out = NULL
        }
        out
})
exp_params = exp_params %>% mutate(sc = map(data, "sc"))
exp_params = exp_params %>% mutate(tr = map(data, "tr"))
saveRDS(exp_params, file = paste0(output_dir, "exp_data_1rep_sim.rds"))

plan(multisession, workers = 12)
exp_params$res = future_map(exp_params$tr, function(tr) {
        res = ice_fase_mod(tr,
                           get_type_from_id(tr$tip.label),
                           total_time = 11.5 - 0.6,
                           root_time = 0.6,
                           theta = 0.)
        res
}, .progress = T, .options = furrr_options(seed = T))
exp_params$res = map2(exp_params$big_graph_id, exp_params$res, function(i, res) {
        res = evaluate_gr_res(res, tree_panel$gr[[i]])
        res
})
exp_params$kc0 = map_dbl(exp_params$res, function(x) {
        x$gr_eval$kc0
})
exp_params$bsum = tree_panel$bsum[exp_params$big_graph_id]
exp_params$num_tip = tree_panel$num_tip[exp_params$big_graph_id]
exp_params$qfm_eval = map(1:nrow(exp_params), function(i) {
        message(i)
        eval_out = evaluate_qfm(exp_params$res[[i]],
                                tree_panel$type_graph[[exp_params$big_graph_id[i]]],
                                exp_params$data[[i]])
        eval_out
})
exp_params = saveRDS(exp_params, file = paste0(output_dir, "exp_data_1rep_proc.rds"))

# generate control
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))
exp_params = readRDS(paste0(output_dir, "exp_data_1rep_proc.rds"))

exp_params = generate_control(exp_params, tree_panel)
expm1_trans <-  function() trans_new("expm1", "expm1", "log1p")
method_col_v1 = c("sps" = "#377eb8",
                  "ice_fase" = "#e41a1c",
                  "random_control" = "#949494")

(exp_params %>% transmute(sampling = sampling,
                          num_tip = num_tip,
                          bsum = bsum,
                          sps = map_dbl(exp_params$sps_gr_eval, function(x) {
                                  if (is.null(x)) {
                                          return(NA)
                                  } else {
                                          return(x$kc0)
                                  }
                          }),
                          random_control = map_dbl(exp_params$gr_control_eval, function(x) {
                                  if (is.null(x)) {
                                          return(NA)
                                  } else {
                                          return(x$kc0)
                                  }
                          }),
                          # ice_fase_phylidite = map_dbl(exp_params$gr3_eval, function(x) {
                          #   if (is.null(x)) {
                          #     return(NA)
                          #     } else {
                          #       return(x$kc0)
                          #       }
                          #   }),
                          # ice_fase_hamming = map_dbl(exp_params$gr5_eval, function(x) {
                          #   if (is.null(x)) {
                          #     return(NA)
                          #     } else {
                          #       return(x$kc0)
                          #       }
                          #   }),
                          ice_fase = map_dbl(exp_params$res, function(x) {
                                  if (is.null(x)) {
                                          return(NA)
                                  } else {
                                          return(x$gr_eval$kc0)
                                  }
                          })
) %>%
                gather(key = "method", value = "rf", -c("sampling", "num_tip", "bsum")) %>%
                ggscatter(x = "bsum", y = "rf", color = "method",
                          facet.by = c("sampling", "num_tip"),
                          ylab = "KC0",
                          xlab = "BSUM",
                          xlim = c(0, NA),
                          scales = "free",
                          size = 0.1,
                          alpha  = 0.2) +
                scale_color_manual(values = method_col_v1) +
                geom_smooth(aes(color = method, group = method),
                            size = 0.5, se = F,
                            span = 1, method = "loess") +
                scale_y_continuous(trans=log1p_trans()) +
                coord_trans(y = expm1_trans()) +
                theme(text = element_text(size = 12),
                      legend.position = "bottom",
                      axis.text.x = element_text(size = 9, angle = 65),
                      strip.background = element_blank())) %>%
        push_pdf(file_name = "eval_topology", width = 5.85, height = 4.55, ps = 12, dir = plot_dir)

# (ggscatter(exp_params, x = "bsum", y = "kc0", facet.by = c("sampling", "num_tip"), scale = "free_x")  +
#         geom_smooth()) %>%
#         push_pdf("temp_kc0", dir = plot_dir)
#
# (ggscatter(exp_params, x = "bsum", y = "kc0", color = "sampling", facet.by = c("num_tip"), scale = "free_x")  +
#                 geom_smooth(aes(color = sampling))) %>%
#         push_pdf("temp_kc0_col", dir = plot_dir)

#### SPS part ####
exp_params$sps = future_map(exp_params$tr, function(tr_r) {
        if (is.null(tr_r)) {
                return(NULL)
        }
        shared_progenitor_score(tr_r,
                                sort(unique(get_type_from_id(tr_r$tip.label))))
}, .progress = T)
exp_params$sps_dist = map(exp_params$sps, function(sps_mat) {
        if (is.null(sps_mat)) {
                return(NULL)
        }
        dmat = 1 - sps_mat/max(sps_mat)
        dmat
})
total_time = 11.5 - 0.6
exp_params$sps_gr = map(exp_params$sps_dist, function(x) {
        if (is.null(x)) {
                return(NULL)
        }
        gr = phangorn::upgma(as.dist(x))
        total_depth = max(node.depth.edgelength(gr))
        gr$edge.length = gr$edge.length / total_depth * total_time #TODO: what to use for sps as max dist?
        gr = name_nodes(gr)
        gr
})
# exp_params$sps_nj_gr = map(exp_params$sps_dist, function(x) {
#         if (is.null(x)) {
#                 return(NULL)
#         }
#         phytools::midpoint.root(phangorn::NJ(as.dist(x)))
# })
exp_params$sps_gr_eval = future_map2(exp_params$sps_gr, exp_params$big_graph_id,
                                     function(x, i) {
                                             if (is.null(x)) {
                                                     return(NULL)
                                             }
                                             evalute_gr(x, tree_panel$gr[[i]])
                                     }, .progress = T, .options = furrr_options(seed = T))
#### end sps part ####
wrap_ccs_topology <- function(sc_mat) {
        sc_mat_mod = do.call(cbind, purrr::map(1:ncol(sc_mat), function(j) {
                paste0("site", j, "_", sc_mat[, j])
        }))
        rownames(sc_mat_mod) = rownames(sc_mat)
        sc_mat_onehot = barcode_allele_onehot_new(sc_mat_mod)
        type_ind_mat = do.call(rbind, purrr::map(split(1:nrow(sc_mat_onehot), get_type_from_id(rownames(sc_mat_mod))), function(x) {
                apply(sc_mat_onehot[x, ], 2, any)
        }))
        score_tb = as_tibble(t(combn(rownames(type_ind_mat), 2)))
        score_tb = bind_rows(score_tb, tibble(V1 = rownames(type_ind_mat),
                                              V2 = rownames(type_ind_mat)))
        score_tb$dist = map2_dbl(score_tb$V1, score_tb$V2, function(x, y) {
                sum(type_ind_mat[x, ] * type_ind_mat[y, ])
        })
        score_tb = mutate(score_tb, N1 = V1, N2 = V2, Var1 = V1, Var2 = V2)
        ccs_mat = dist_df2mat(score_tb)
        ccs_rows = rowSums(ccs_mat)
        ccs_cols = colSums(ccs_mat)
        ccs_exp = matrix(1, nrow = nrow(ccs_mat) , ncol = ncol(ccs_mat))
        for (i in 1:nrow(ccs_exp)) {
                ccs_exp[i, ] = ccs_exp[i, ] * ccs_rows[i]
        }
        for (j in 1:nrow(ccs_exp)) {
                ccs_exp[, j] = ccs_exp[, j] * ccs_cols[j]
        }
        ccs_exp = ccs_exp / sum(ccs_mat)
        ccs_norm = ccs_mat / ccs_exp
        hcl = hclust(dist(ccs_norm), method = "average")
        gr = dendextend:::as.phylo.dendrogram(as.dendrogram(hcl))
        gr
}
sc_dump_dir = paste0(output_dir, "1rep_sc/")
dir.create(sc_dump_dir)
for (j in 1:nrow(exp_params)) {
        sc_mat = exp_params$data[[j]]$sc
        saveRDS(sc_mat, file = paste0(sc_dump_dir, stringr::str_pad(j, pad = "0", width = 4), ".rds"))
}
ccs_dump_dir = paste0(output_dir, "1rep_ccs/")
dir.create(ccs_dump_dir)
future_walk(1:nrow(exp_params), function(j) {
        sc_mat = readRDS(paste0(sc_dump_dir, stringr::str_pad(j, pad = "0", width = 4), ".rds"))
        gr = wrap_ccs_topology(sc_mat)
        saveRDS(gr, file = paste0(ccs_dump_dir, stringr::str_pad(j, pad = "0", width = 4), ".rds"))
}, .progress = T)

csc_kc0 = map_dbl(1:nrow(exp_params), function(j) {
        data_file = paste0(ccs_dump_dir, stringr::str_pad(j, pad = "0", width = 4), ".rds")
        if (file.exists(data_file)) {
                gr = readRDS(data_file)
                evalute_gr(gr, tree_panel$gr[[exp_params$big_graph_id[j]]])$kc0
        } else {
                return(NA)
        }
})

bind_rows(tibble(method = "Clonal coupling score",
                 kc0 = csc_kc0[!is.na(csc_kc0)]),
          tibble(method = "ICE_FASE",
                 kc0 = exp_params$kc0[!is.na(csc_kc0)])) %>%
        ggboxplot(x = "method", y = "kc0", xlab = "", ylab = "KC0")



library(phylogram)
library(dendextend)
evalute_gr(gr, tree_panel$gr[[1]])


# exp_params = run_fase_eval(exp_params, tree_panel)

# when absolute size is small, mean fase can be wrong by chance
fase_eval_all = bind_rows(map(1:nrow(exp_params), function(i) {
        mutate(exp_params$fase_eval[[i]], exp_id = i)
}))
saveRDS(fase_eval_all, file = "./intermediate_data/panel_mod2_v1/fase_eval_all.rds")
# fase_large_error = fase_eval_all[fase_eval_all$log2_state_sample_frac == 0 & fase_eval_all$error > 4, ]
# fase_large_error$data[[1]]

fase_eval_all = readRDS("./intermediate_data/panel_mod2_v1/fase_eval_all.rds")
(fase_eval_all %>%
                mutate(log2_n_fase = log2(n_fase)) %>%
                sample_n(1e4) %>%
                ggscatter(x = "log2_state_sample_size",
                          y = "log2_n_fase") + geom_smooth(method = "lm")) %>%
        push_pdf("n_fase_psize", dir = "./plots/panel_mod2_v1/")

fase_eval_all %>% group_by(exp_id) %>%
        summarise(mean_n_fase = mean(n_fase)) %>% summarise(median(mean_n_fase))

fase_eval_all$num_tip = exp_params$num_tip[fase_eval_all$exp_id]
(fase_eval_all %>% group_by(exp_id) %>%
                summarise(mean_n_fase = mean(n_fase)) %>%
                gghistogram(x = "mean_n_fase", fill = "grey", xlab = "Average number of FASEs per type")) %>%
        push_pdf("n_fase_hist", dir = "./plots/panel_mod2_v1/")


exp_params %>% group_by(sampling, num_tip) %>%
        summarise(mean(kc0 == 0))

eval_tb$reslv = ifelse(eval_tb$kc0 == 0, "resolved", "non-reseolved")
eval_tb %>%
        group_by(exp_id, num_tip, reslv, sampling) %>%
        summarise(min_log2_sample_frac = min(log2_node_sampled, na.rm = T)) %>%
        ggboxplot(x = "reslv", y = "min_log2_sample_frac", facet.by = c("num_tip", "sampling"))

eval_tb$log2_gr_node_size_in = eval_tb$gr_node_size_in
eval_tb$log2_node_size_in = eval_tb$log2_node_size

produce_strat_time_sum(eval_tb, group_by = c("suff_sampled", "sampling"))
produce_strat_size_sum(eval_tb, group_by = c("suff_sampled", "sampling"))
produce_strat_split_sum(eval_tb, group_by = c("suff_sampled", "sampling"))

produce_strat_time_sum(eval_tb, group_by = c("sampling"))
produce_strat_size_sum(eval_tb, group_by = c("sampling"))
produce_strat_split_sum(eval_tb, group_by = c("sampling"))


fase_eval_all %>%
        sample_n(5e4) %>%
        ggscatterhist(x = "log2_state_sample_frac", y = "error", alpha = 0.1) +
        geom_smooth()

fase_eval_all %>%
        sample_n(1e5) %>%
        ggplot(aes(x = log2_state_sample_frac, y = error)) +
        # geom_point(alpha = 0.5) +
        geom_density_2d_filled(h = c(1, 1))
# ggscatterhist(x = "log2_state_sample_frac", y = "error", alpha = 0.1) +
# geom_smooth()

fase_eval_all$log2_state_sample_frac_bin = cut(fase_eval_all$log2_state_sample_frac,
                                               breaks = c(-Inf, seq(-10, 0, by = 1.))
)
(fase_eval_all %>%
                sample_n(1e5) %>%
                ggboxplot(x = "log2_state_sample_frac_bin", y = "error",
                          xlab = "log2(Progenitor state sampling fraction",
                          ylab = "FASE time error",
                          size = 0.5,
                          outlier.size = 0.01) +
                theme(text = element_text(size = 8),
                      axis.text.x = element_text(angle = 90, size = 6))) %>%
        push_pdf(file_name = "fase_time_error", w = 2., h = 2, dir = plot_dir)

(fase_eval_all %>%
                group_by(log2_state_sample_frac_bin) %>%
                sample_n(500, replace = T) %>%
                ggline(x = "log2_state_sample_frac_bin", y = "error",
                       xlab = "log2(Progenitor state sampling fraction",
                       ylab = "FASE time error",
                       size = 0.5,
                       add = c("mean_sd", "jitter"),
                       add.params = list(alpha = 0.02),
                       ylim = c(0, 7.5)) +
                theme(text = element_text(size = 8),
                      axis.text.x = element_text(angle = 90, size = 6))) %>%
        push_pdf(file_name = "fase_time_error_line_dot", w = 4., h = 4, dir = plot_dir)

bind_rows(exp_params$fase_eval) %>%
        filter(log2_state_size > 8) %>%
        sample_n(1e4) %>%
        ggscatter(x = "log2_state_sample_frac",
                  y = "error", alpha = 0.05,
                  xlab = "log2(Progenitor state sampling fraction",
                  ylab = "FASE time error") +
        geom_smooth(method = "loess")

bind_rows(exp_params$fase_eval) %>%
        filter(log2_state_size > 4) %>%
        sample_n(1e4) %>%
        ggscatter(x = "log2_state_sample_frac",
                  y = "error", alpha = 0.05,
                  xlab = "log2(Progenitor state sampling fraction",
                  ylab = "FASE time error") + geom_smooth()

lm(error ~ log2_state_size + log2_state_sample_frac,
   data = bind_rows(exp_params$fase_eval)) %>% summary


tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))
exp_params = readRDS(paste0(output_dir, "exp_data_1rep_proc.rds"))
sampling_col = RColorBrewer::brewer.pal(9, "Spectral")
eval_tb = bind_rows(purrr::map(1:nrow(exp_params), function(i) {
        mutate(exp_params$qfm_eval[[i]]$recon,
               exp_id = i,
               sampling = exp_params$sampling[i],
               num_tip = exp_params$num_tip[i])
}))
eval_tb$suff_sampled = eval_tb$log2_node_sampled > -2

eval_tb_true = bind_rows(map(1:nrow(exp_params), function(i) {
        mutate(exp_params$qfm_eval[[i]]$true,
               exp_id = i,
               sampling = exp_params$sampling[i],
               num_tip = exp_params$num_tip[i])
}))
eval_tb_true %>% group_by(sampling, num_tip) %>%
        summarise(accuracy = mean(!is.na(node_gr)))

eval_tb_true$log2_ps_cov = log2(map_dbl(eval_tb_true$node_size_collect, function(x) {
        if (is.null(x)) {
                return(NA)
        } else {
                return(x)
        }
}) / eval_tb_true$gr_node_size_in)
table(eval_tb_true$log2_node_sampled > -2,
      eval_tb_true$log2_ps_cov > log2(5))

eval_tb_true %>%
        sample_n(1e4) %>%
        ggscatter(x = "node_time", y = "is_resolved",
                  facet.by = "sampling") +
        geom_smooth()
# eval_tb_true %>%
#         filter(sampling == "fixed" & num_tip == 16) %>%
#         ggscatter(x = "log2_node_sampled", y = "is_resolved") +
#         geom_smooth()
eval_tb_true %>%
        sample_n(1e4) %>%
        ggscatter(y = "log2_node_sampled", x = "node_time") +
        geom_smooth()

eval_tb_true %>%
        ggscatter(x = "node_time", y = "is_resolved") +
        geom_smooth()

eval_tb_true$time_bin = cut(eval_tb_true$node_time, breaks = seq(2, 11, by = 0.5))
eval_tb_true %>%
        sample_n(1e4) %>%
        ggboxplot(x = "time_bin", y = "log2_node_sampled",
                  xlab = "Commitment time",
                  ylab = "log2(Sampling fraction)") +
        theme(axis.text.x = element_text(angle = 90))


hist(eval_tb$log2_node_sampled, breaks = 40)
mean(eval_tb$gr_time_trans_error^2)



g1 = (eval_tb %>%
              filter(!is.na(node)) %>%
              sample_n(10000) %>%
              ggscatter(x = "node_time", y=  "gr_time_trans",
                        color = "log2_node_sampled",
                        size = 0.15,
                        alpha = 0.5,
                        # facet.by = c("sampling", "num_tip"),

                        # xlab = "True commitment time",
                        # ylab = "Inferred commitment time",
                        xlab = "", ylab = ""
              ) %>% facet(facet.by = c("sampling"), ncol = 1) +
              geom_abline(size = 0.25)  +
              geom_smooth(aes(x =  node_time, y = gr_time_trans),
                          color = sampling_col[9],
                          method = "loess",
                          se = F,
                          span = 1,
                          data = filter(eval_tb, suff_sampled == "TRUE")) +
              geom_smooth(aes(x =  node_time, y = gr_time_trans),
                          color = sampling_col[1],
                          method = "loess",
                          span = 1,
                          se = F,
                          data = filter(eval_tb, suff_sampled == "FALSE")) +
              stat_cor(aes(group = suff_sampled,
                           label =  ..r.label..),
                       method = "spearman",
                       cor.coef.name = "rho",
                       digits = 3) +
              scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                                     breaks = (-4):(-1), na.value = sampling_col[1]) +
              theme(text = element_text(size = 10),
                    panel.background = element_rect(fill = "#bebebe"),
                    # legend.position = "top",
                    legend.position = "none",
                    strip.background = element_blank()))

eval_tb = eval_tb %>%
        # filter(suff_sampled) %>%
        mutate(log2_node_size = log2(node_size),
               log2_gr_node_size_in = log2(gr_node_size_in))

# saveRDS(exp_params, file = "./intermediate_data/panel_mod2/exp_data_1rep_proc.rds")
# eval_tb$log2_node_size = eval_tb$log2_node_size

g2 = (eval_tb %>%
              ggscatter(x = "log2_node_size",
                        y = "log2_gr_node_size_in",
                        size = 0.15,
                        alpha = 0.5,
                        # facet.by = c("sampling", "num_tip"),
                        color = "log2_node_sampled")  %>%
              facet(facet.by = c("sampling"), ncol = 1) +
              geom_abline() +
              # geom_smooth() +
              geom_smooth(aes(x =  log2_node_size, y = log2_gr_node_size_in),
                          color = sampling_col[9],
                          method = "loess",
                          se = F,
                          span = 1,
                          data = filter(eval_tb, suff_sampled == "TRUE")) +
              geom_smooth(aes(x =  log2_node_size, y = log2_gr_node_size_in),
                          color = sampling_col[1],
                          method = "loess",
                          span = 1,
                          se = F,
                          data = filter(eval_tb, suff_sampled == "FALSE")) +
              xlab("") +
              ylab("") +
              # xlab("log2(Progenitor population size)") +
              # ylab("log2(Inferred progenitor population size)") +
              stat_cor(aes(group = suff_sampled,
                           label =  ..r.label..),
                       method = "spearman",
                       cor.coef.name = "rho",
                       digits = 3) +
              scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                                     breaks = (-4):(-1), na.value = sampling_col[1]) +
              theme(text = element_text(size = 10),
                    panel.background = element_rect(fill = "#bebebe"),
                    legend.position = "none",
                    strip.background = element_blank()))


g3 = (eval_tb %>%
              # filter(suff_sampled) %>%
              ggscatter(x = "node_split_order", y = "gr_node_split_order",
                        # facet.by = c("sampling", "num_tip"),
                        color = "log2_node_sampled",
                        alpha = 0.5,
                        size = 0.15,
                        # xlab = "True commitment ratio",
                        # ylab = "Inferred commitment ratio"
                        xlab = "",
                        ylab = ""
              ) %>%
              facet(facet.by = c("sampling"), ncol = 1) +
              geom_abline() +
              geom_smooth(aes(x =  node_split_order, y = gr_node_split_order),
                          color = sampling_col[9],
                          method = "loess",
                          se = F,
                          span = 1,
                          data = filter(eval_tb, suff_sampled == "TRUE")) +
              geom_smooth(aes(x =  node_split_order, y = gr_node_split_order),
                          color = sampling_col[1],
                          method = "loess",
                          span = 1,
                          se = F,
                          data = filter(eval_tb, suff_sampled == "FALSE")) +
              stat_cor(aes(group = suff_sampled,
                           label =  ..r.label..),
                       method = "spearman",
                       cor.coef.name = "rho",
                       digits = 3) +
              scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                                     breaks = (-4):(-1), na.value = sampling_col[1]) +
              theme(text = element_text(size = 10),
                    panel.background = element_rect(fill = "#bebebe"),
                    legend.position = "none",
                    strip.background = element_blank()))

push_pdf(g1, file_name = "time_error", w = 3, h = 5, dir = plot_dir)
push_pdf(g2, file_name = "size_error", w = 3, h = 5, dir = plot_dir)
push_pdf(g3, file_name = "split_error", w = 3, h = 5, dir = plot_dir)

eval_tb$log2_ps_cov = log2(map_dbl(eval_tb$node_size_collect, function(x) {
        if (is.null(x)) {
                return(NA)
        } else {
                return(x)
        }
}) / eval_tb$gr_node_size_in)
eval_tb$abs_time_trans_error = abs(eval_tb$gr_time_trans_error)
eval_tb$abs_gr_node_size_logfc = abs(eval_tb$gr_node_size_logfc)
eval_tb$abs_gr_node_split_error = abs(eval_tb$gr_node_split_order - eval_tb$node_split_order)

# auROC
library(pROC)
pROC_obj <- roc(eval_tb$log2_node_sampled > -2,
                eval_tb$log2_ps_cov,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

plot(pROC_obj, print.auc=TRUE) %>%
        push_pdf("sampling_coverage_auc", w = 2.5, h = 2.5, ps = 10, dir = plot_dir)

# psCOV
(eval_tb %>%
                sample_n(1e4) %>%
                ggscatter(x = "log2_ps_cov",
                          y = "log2_node_sampled",
                          xlab = "log2(PScov)",
                          ylab = "log2(Progenitor state sampling fraction)",
                          alpha = 0.2,
                          size = 0.05) + geom_smooth(size = 0.5, se = F) +
                geom_rect(xmin = log2(5.), xmax = Inf, ymin = -2, ymax = 0.2, fill = NA, color = "red") +
                theme(text = element_text(size = 10))) %>%
        push_pdf("sampling_pscov_scatter", w = 3.0, h = 3.0, ps = 12, dir = plot_dir)


g_s = sampling_plots[[1]] + geom_vline(xintercept = -2)
g_hist = eval_tb %>% gghistogram(x = "log2_node_sampled", fill = "#6d6d6d", color = NA, bins = 50,
                                 xlim = c(-10, 0),
                                 xlab = "log2(Sampling fraction)") +
        geom_vline(xintercept = -2)
((g_s / g_hist) + plot_layout(heights = c(3, 1))) %>%
        push_pdf("time_error_sample_frac_whist", w = 3.5, h = 4.5, dir = plot_dir)

g_hist_alt = eval_tb %>% gghistogram(x = "log2_node_sampled", fill = "#6d6d6d", color = NA, bins = 50,
                                     xlim = c(-10, 0),
                                     scales = "free_y",
                                     xlab = "log2(Sampling fraction)") %>%
        facet(facet.by = "sampling", ncol = 1) +
        geom_vline(xintercept = -2)
(g_hist_alt) %>%
        push_pdf("time_error_sample_frac_alt", w = 3.5, h = 3.5, dir = plot_dir)

sampling_plots[[1]] %>%
        push_pdf("time_error_sample_frac", w = 3.5, h = 3.375, dir = plot_dir)
sampling_plots[[3]] %>%
        push_pdf("size_error_sample_frac", w = 3.5, h = 3.375, dir = plot_dir)
sampling_plots[[5]] %>%
        push_pdf("split_error_sample_frac", w = 3.5, h = 3.375, dir = plot_dir)

sampling_plots = list(eval_tb %>%
                              filter(!is.na(node)) %>%
                              sample_n(size = 5000) %>%
                              ggscatter(x = "log2_node_sampled",
                                        y = "abs_time_trans_error",
                                        xlab = "log2(Sampling fraction)",
                                        ylab = "Absolute commitment time error",
                                        xlim = c(-10, 0),
                                        alpha = 0.2,
                                        size = 0.15) + geom_smooth(method = "loess"),
                      eval_tb %>%
                              filter(!is.na(node)) %>%
                              sample_n(size = 5000) %>% ggscatter(x = "log2_ps_cov",
                                                                  y = "abs_time_trans_error",
                                                                  xlab = "log2(PS_cov)",
                                                                  ylab = "Absolute commitment time error",
                                                                  xlim = c(0, 11),
                                                                  alpha = 0.05,
                                                                  size = 0.15) + geom_smooth(method = "loess"),
                      eval_tb %>%
                              filter(!is.na(node)) %>%
                              sample_n(size = 5000)%>% ggscatter(x = "log2_node_sampled",
                                                                 y = "abs_gr_node_size_logfc",
                                                                 xlab = "log2(Sampling fraction)",
                                                                 ylab = "Absolute log2(Progenitor population size fold change)",
                                                                 xlim = c(-10, 0),
                                                                 alpha = 0.2,
                                                                 size = 0.15) + geom_smooth(method = "loess"),
                      eval_tb %>%
                              filter(!is.na(node)) %>%
                              sample_n(size = 5000)%>% ggscatter(x = "log2_ps_cov",
                                                                 y = "abs_gr_node_size_logfc",
                                                                 xlab = "log2(PS_cov)",
                                                                 ylab = "Absolute log2(Progenitor population size fold change)",
                                                                 xlim = c(0, 11),
                                                                 alpha = 0.2,
                                                                 size = 0.15) + geom_smooth(method = "loess"),
                      eval_tb %>%
                              filter(!is.na(node)) %>%
                              sample_n(size = 5000)%>% ggscatter(x = "log2_node_sampled",
                                                                 y = "abs_gr_node_split_error",
                                                                 xlab = "log2(Sampling fraction)",
                                                                 ylab = "Absolute commitment bias absolute error",
                                                                 xlim = c(-10, 0),
                                                                 alpha = 0.2,
                                                                 size = 0.15) + geom_smooth(method = "loess"),
                      eval_tb %>%
                              filter(!is.na(node)) %>%
                              sample_n(size = 5000) %>% ggscatter(x = "log2_ps_cov",
                                                                  y = "abs_gr_node_split_error",
                                                                  xlab = "log2(PS_cov)",
                                                                  ylab = "Absolute commitment bias absolute error",
                                                                  xlim = c(0, 11),
                                                                  alpha = 0.2,
                                                                  size = 0.15) + geom_smooth(method = "loess")
)

eval_tb %>%
        filter(!is.na(node)) %>%
        sample_n(size = 5000) %>% ggscatter(x = "log2_ps_cov",
                                            y = "abs_gr_node_size_logfc",
                                            color = "temp",
                                            xlab = "log2(PS_cov)",
                                            ylab = "Absolute log2(Progenitor population size fold change)",
                                            xlim = c(0, 11),
                                            alpha = 0.2,
                                            size = 0.15) +
        geom_smooth(method = "loess")


sampling_plots[[2]] %>%
        push_pdf("time_error_pscov", w = 3.5, h = 3.375, dir = plot_dir)
sampling_plots[[4]] %>%
        push_pdf("size_error_pscov", w = 3.5, h = 3.375, dir = plot_dir)
sampling_plots[[6]] %>%
        push_pdf("split_error_pscov", w = 3.5, h = 3.375, dir = plot_dir)


eval_tb_true = bind_rows(map(1:nrow(exp_params), function(i) {
        mutate(exp_params$qfm_eval[[i]]$true,
               exp_id = i,
               sampling = exp_params$sampling[i],
               num_tip = exp_params$num_tip[i])
}))
eval_tb_true$bsum = tree_panel$bsum[exp_params$big_graph_id[eval_tb_true$exp_id]]
eval_tb_true$bsum_bin = tree_panel$bsum_bin[exp_params$big_graph_id[eval_tb_true$exp_id]]

eval_tb_true_sum = eval_tb_true %>% group_by(exp_id) %>%
        summarise(mean_log2_sampling = mean(log2_node_sampled),
                  mean_undersample = mean(log2_node_sampled < - 3.),
                  kc0 = kc0[1],
                  sampling = sampling[1],
                  num_tip = num_tip[1],
                  bsum = bsum[1])

x1 = eval_tb_true_sum %>%
        ggscatter(x = "bsum", y = "mean_undersample", facet.by = c("sampling", "num_tip"),
                  scales = "free_x") +
        geom_smooth()

x2 = eval_tb_true_sum %>% ggscatter(x = "mean_log2_sampling", y = "kc0", facet.by = c("sampling", "num_tip")) +
        geom_smooth()

x1 / x2




# total_time = 11.5 - 0.6
# for (j in 1:nrow(exp_params)) {
#         message(j)
#         print(paste0("Ntip: ", length(all_graphs[[exp_params$big_graph_id[j]]]$tip_id)))
#         tr3 = name_nodes(phylotime(exp_params$sc[[j]], total_time, mut_p = exp_params1$mut_p[[j]]))
#         saveRDS(tr3,
#                 file = paste0(output_dir, "temp_tr3/", stringr::str_pad(j, width = 3, side = "left", "0")))
# }
# exp_params = mutate(exp_params, tr3 = map2(sc, mut_p, function(x, m) {
#         name_nodes(phylotime(x, total_time, mut_p = m))
# }))
# saveRDS(exp_params, file = paste0(output_dir, "exp_data_proc.rds"))


xx[xx$from_type == "Node-1" & xx$to_type != "Node-1", ]

# test edge depth function
tr = gr_tr_data$tr
tr_time_alt = gr_tr_data$tr_edges_tb$to_time
names(tr_time_alt) = gr_tr_data$tr_edges_tb$to

tr_node_depth = node.depth.edgelength(tr)
names(tr_node_depth) = c(tr$tip.label, tr$node.label)

tr_node_depth[["type_15_gen_0_1"]]

gr_tr_data$tr_edges_tb[gr_tr_data$tr_edges_tb$from == "type_15_gen_0_1", ]
gr_tr_data$tr_edges_tb[gr_tr_data$tr_edges_tb$from == "type_15_gen_1_2", ]

plot(tr_time_alt[names(tr_node_depth)], tr_node_depth)
abline(0, 1)

qplot(sort(node.depth.edgelength(as.phylo_mod2.type_graph(type_graph))) + 2.4,
      sort(get_true_diff_time_mod3(type_graph)) + 0.6) +
        geom_abline(slope = 1, intercept = 0.)

qplot(sort(exp_obj$gr_trans_time), sort(get_true_diff_time_mod3(type_graph))) +
        geom_abline(slope = 1, intercept = 0.)


# trying to resolve the small timing error here
# clarify definitions and check each function
# First, each node is one-to-one with an edge that lead to it, the edge has "to" field being the node
# the potency of the node is determined by its subtree
# the time of a node is defined as the end time of the edge

# For state "15", fate becomes pure at time 3.0, true commitment time is 2.4
# the true diff time is for type "15" is 2.4
# get_true_diff_time_mod3(type_graph)

# now are the ICE nodes correctly identified? should be gen3. Seems correct



# for example: type_15_gen_3_1 is an ICE node, its node time is the commitment time
#

# output data for batch phylotime
exp_params$target_time = map_dbl(exp_params$big_graph_id, function(graph_id) tree_panel$type_graph[[graph_id]]$target_time)
push_sc_tb <- function(exp_params, dir) {
        out = select(exp_params,
                     sc, `target_time`, mut_p,
                     big_graph_id, sample_size, num_element, sampling, i_sim,
                     bsum, num_tip)
        saveRDS(out, file = paste0(dir, "phylotime_input.rds"))
}
push_sc_tb(exp_params, output_dir)
run_batch_phylotime <- function(exp_params, i, output_dir ,root_edge = 0.6) {
        if (!dir.exists(output_dir)) {
                dir.create(output_dir)
        }
        tr = phylotime(exp_params$sc[[i]],
                       t_total = exp_params$target_time[i] - root_edge,
                       mut_p = exp_params$mut_p[[i]],
                       parallel = T)
        out_file = paste0(output_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        saveRDS(tr, file = out_file)
}
run_batch_phylotime(exp_params, 1, paste0(output_dir, "tr3/"))


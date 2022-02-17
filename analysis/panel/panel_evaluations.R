method_col = c("sps x truth" = "#11468F",
               "sps" = "#11468F",
               "ice_fase x hamming" = "#FC28FB",
               "ice_fase_hamming" = "#FC28FB",
               "ice_fase x phylidite" = "#FFC900",
               "ice_fase_phylidite" = "#FFC900",
               "ice_fase x truth" = "#DA1212",
               "ice_fase" = "#DA1212")
plot_col_vec(method_col)

output_dir = "./intermediate_data/panel/"
exp_params = readRDS(paste0(output_dir, "exp_data_rep1_2_proc2.rds"))
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))

# FIGURE 2: panel of all graphs
# type_graph = all_graphs[[1]]
# node_label_letter <- function(x) {
#   t_ind = grepl("-", x)
#   out = x
#   out[t_ind] = stringr::str_replace_all(out[t_ind], "-", "T")
#   out[!t_ind] = paste0("P", out[!t_ind])
#   out
# }
# g_all = map(all_graphs, function(type_graph) {
#   out = plot_type_graph(type_graph, function(x) {
#     gr_color(type_graph)[x]
#   }, node_label_mapper = node_label_letter,
#   show_node_text = F)
#   if (length(type_graph$tip_id) > 16) {
#     out = out  + theme(axis.title = element_blank(),
#                        axis.line.y = element_blank(),
#                        axis.text.y = element_blank(),
#                        axis.ticks.y = element_blank())
#   }
#   out
# })

#### Figure #2 Plotting example graphs from the panel ####
g_all_simple = map(all_graphs, function(type_graph) {
  out = plot_type_graph_clean(type_graph,
                              node_col_mapper = function(x) {
                                gr_color(type_graph)[x]
                              },
                              node_size = 2.5)
  if (length(type_graph$tip_id) > 16) {
    out = out  + theme(axis.title = element_blank(),
                       axis.line.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())
  }
  out
})
all_graph_num_tips = map_dbl(all_graphs, function(x) length(x$tip_id))
plot_indices = c(c(1, 9, 17),
                 c(18, 34, 50),
                 c(51, 89, 127))
plot_indices = c(matrix(plot_indices, ncol = 3, byrow = T))
(reduce(g_all_simple[plot_indices], `+`) + plot_layout(nrow = 3, width = c(1, 2, 4))) %>%
  push_pdf(file_name = "panel", dir = "./plots/panel/",
           width = 10.5, height = 6.,
           ps = 12)
#### End Figure 2 ####

# BEGIN Evaluations Below #
expm1_trans <-  function() trans_new("expm1", "expm1", "log1p")
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
                          ice_fase_phylidite = map_dbl(exp_params$gr3_eval, function(x) {
                            if (is.null(x)) {
                              return(NA)
                              } else {
                                return(x$rf)
                                }
                            }),
                          ice_fase_hamming = map_dbl(exp_params$gr5_eval, function(x) {
                            if (is.null(x)) {
                              return(NA)
                              } else {
                                return(x$rf)
                                }
                            }),
                          ice_fase = map_dbl(exp_params$gr_eval, function(x) {
                            if (is.null(x)) {
                              return(NA)
                              } else {
                                return(x$rf)
                              }
                            })) %>%
        gather(key = "method", value = "rf", -c("sampling", "num_tip", "bsum")) %>%
  ggscatter(x = "bsum", y = "rf", color = "method",
            facet.by = c("sampling", "num_tip"),
            ylab = "KC0",
            xlab = "BSUM",
            xlim = c(0, NA),
            scales = "free", size = 0.1) +
    scale_color_manual(values = method_col) +
  geom_smooth(aes(color = method, group = method),
              size = 0.5, se = F,
              span = 1, method = "loess") +
  scale_y_continuous(trans=log1p_trans()) +
  coord_trans(y = expm1_trans()) +
  theme(text = element_text(size = 12),
        legend.position = "right",
        axis.text.x = element_text(size = 9, angle = 65),
        strip.background = element_blank())) %>%
  push_pdf(file_name = "eval_ice_fase_phylidite_kc0", width = 5.5, height = 2.5, ps = 12, dir = "./plots/panel/")

# Generating evaluations on the experiment
# evaluations with truth:
eval_tb_gr = generate_evaluate_tb(exp_params, all_graphs, tr_col = "tr", gr_col = "gr", gr_eval_col = "gr_eval")
# evaluations with phylotime:
eval_tb_gr3 = generate_evaluate_tb(exp_params, all_graphs, tr_col = "tr3", gr_col = "gr3", gr_eval_col = "gr3_eval")
save(eval_tb_gr, file = paste0(output_dir, "./intermediate_data/panel/exp_gr_eval.rda"))
save(eval_tb_gr3, file = paste0(output_dir, "./intermediate_data/panel/exp_gr3_eval.rda"))

load(paste0(output_dir, "exp_gr_eval.rda"))
source("./analysis/process_evaluation_results.R")

eval_tb_all = eval_tb_wrap(eval_tb_gr, exp_params)

# looking at node resolution vs sampling fraction
# eval_tb_true = process_res(eval_tb_gr, return_true = T)
# eval_tb_true$suff_sampled = eval_tb_true$log2_node_sampled > -2
# eval_tb_true %>%
#   group_by(suff_sampled) %>%
#   summarise(mean(!is.na(node_gr)))
# eval_tb_true$resolved = as.numeric(!is.na(eval_tb_true$node_gr))
# ggscatter(eval_tb_true,
#           x = "log2_node_sampled", y = "resolved") +
#   geom_smooth(method = "loess")

# set this tibble to be plotted
eval_tb_plot = eval_tb_all
plot_dir = "./plots/panel/"

# FIGURE S2: robustness plots
eval_tb_plot = eval_tb_plot %>%
  mutate(log2_est_coverage = log2(1/gr_collect_size_used))
(eval_tb_plot %>% ggscatter(x = "gr_node_mean_final", y = "log2_node_sampled") + geom_smooth()) +
        (eval_tb_plot %>% ggscatter(x = "gr_collect_size_used", y = "log2_node_sampled") + geom_smooth())
(eval_tb_plot %>%
         ggscatter(x = "log2_est_coverage",
                   y = "log2_node_sampled",
                   # xlab = "log2(Estimated progenitor state coverage)",
                   # ylab = "log2(Progenitor state sampling fraction)",
                   xlab = "",
                   ylab = "",
              size = 0.1) + geom_smooth(size = 0.5, se = F) +
    geom_rect(xmin = log2(2.5), xmax = Inf, ymin = -2, ymax = 0.5, fill =NA, color = "red")) %>%
 push_png("sampling_coverage_scatter", w = 2.5, h = 1.5, ps = 12, res = 600, dir = plot_dir)
table(eval_tb_plot$log2_est_coverage > log2(2.5), eval_tb_plot$log2_node_sampled > -2)

library(pROC)
# pROC_obj <- roc(eval_tb_plot$log2_node_sampled > -2,
#                 1 - eval_tb_plot$gr_node_mean_final,
#                 smoothed = TRUE,
#                 # arguments for ci
#                 ci=TRUE, ci.alpha=0.9, stratified=FALSE,
#                 # arguments for plot
#                 plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
#                 print.auc=TRUE, show.thres=TRUE)
pROC_obj <- roc(eval_tb_plot$log2_node_sampled > -2,
                eval_tb_plot$log2_est_coverage,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)
plot(pROC_obj) %>%
  push_pdf("sampling_coverage_auc", w = 2., h = 2.0, ps = 10, dir = plot_dir)


sampling_col = RColorBrewer::brewer.pal(9, "Spectral")

(eval_tb_plot %>%
        # filter(suff_sampled) %>%
        ggscatter(x = "node_time", y=  "gr_time_trans",
                  color = "log2_node_sampled",
                  size = 0.15,
                  alpha = 0.05,
                  # facet.by = c("sampling", "num_tip"),
                  facet.by = c("sampling"),
                  # xlab = "True commitment time",
                  # ylab = "Inferred commitment time",
                  xlab = "", ylab = ""
        ) +
    geom_abline(size = 0.25)  +
    geom_smooth(aes(x =  node_time, y = gr_time_trans),
                color = sampling_col[9],
                method = "loess",
                se = F,
                span = 1,
                data = filter(eval_tb_plot, suff_sampled == "TRUE")) +
    geom_smooth(aes(x =  node_time, y = gr_time_trans),
                color = sampling_col[1],
                method = "loess",
                span = 1,
                se = F,
                data = filter(eval_tb_plot, suff_sampled == "FALSE")) +
    # stat_cor(aes(group = suff_sampled,
    #              label =  ..r.label..),
    #          method = "spearman",
    #          cor.coef.name = "rho",
    #          digits = 3) +
    scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                           breaks = (-4):(-1), na.value = sampling_col[1]) +
                theme(text = element_text(size = 10),
                      # legend.position = "top",
                      legend.position = "none",
                      strip.background = element_blank())) %>%
  # plot_legend()
  # push_png(file_name = "node_time_est_bytip", res = 1200, w = 6.5, h = 2.8, dir = plot_dir)
  push_png(file_name = "node_time_est_merge_nolab", res = 1200, w = 4.0, h = 2.0, dir = plot_dir)

eval_tb_plot = eval_tb_plot %>%
  # filter(suff_sampled) %>%
  mutate(log2_node_size = log2(node_size),
         log2_gr_node_size_in = log2(gr_node_size_in))

(eval_tb_plot %>%
        ggscatter(x = "log2_node_size",
                  y = "log2_gr_node_size_in",
                  size = 0.15,
                  alpha = 0.05,
                  facet.by = c("sampling"),
                  # facet.by = c("sampling", "num_tip"),
                  color = "log2_node_sampled") +
    geom_abline() +
    # geom_smooth() +
    geom_smooth(aes(x =  log2_node_size, y = log2_gr_node_size_in),
                color = sampling_col[9],
                method = "loess",
                se = F,
                span = 1,
                data = filter(eval_tb_plot, suff_sampled == "TRUE")) +
    geom_smooth(aes(x =  log2_node_size, y = log2_gr_node_size_in),
                color = sampling_col[1],
                method = "loess",
                span = 1,
                se = F,
                data = filter(eval_tb_plot, suff_sampled == "FALSE")) +
    xlab("") +
    ylab("") +
    # xlab("log2(Progenitor population size)") +
    # ylab("log2(Inferred progenitor population size)") +
    # stat_cor(aes(group = suff_sampled,
    #              label =  ..r.label..),
    #          method = "spearman",
    #          cor.coef.name = "rho",
    #          digits = 3) +
    scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                           breaks = (-4):(-1), na.value = sampling_col[1]) +
    theme(text = element_text(size = 10), legend.position = "none", strip.background = element_blank())) %>%
  # push_png(file_name = "node_size_est_bytip", res = 1200, w = 6.5, h = 2.8, dir = plot_dir)
push_png(file_name = "node_size_est_merge_nolab", res = 1200, w = 4.0, h = 2.0, dir = plot_dir)

# this plot has not been updated to the latest
# (eval_tb_plot %>%
#         # filter(suff_sampled) %>%
#         mutate(log2_node_size_sampled = log2(node_size_sampled),
#                log2_gr_node_size_in = log2(gr_node_size_in)) %>%
#         ggscatter(x = "log2_node_size_sampled",
#                   y =  "log2_gr_node_size_in",
#                   size = 0.25,
#                   facet.by = c("sampling", "num_tip"),
#                   color = "log2_node_sampled") + geom_abline() + geom_smooth(method = "lm") + stat_cor(aes(label =  ..r.label..)) +
#         xlab("log2(Progenitor sample size)") + ylab("log2(Inferred progenitor population size)") +
#         scale_colour_gradientn(limits = c(-4, 0),
#                                colors = c("red", "white", "blue"),
#                                breaks = c(-4, -0.5, 0),
#                                na.value = "red") +
#         theme(text = element_text(size = 10), legend.position = "none")) %>%
#         push_png(file_name = "node_size_sampled_est_all", res = 1200, w = 4.5, h = 3., dir = plot_dir)

(eval_tb_plot %>%
        # filter(suff_sampled) %>%
        ggscatter(x = "node_split_order", y = "gr_node_split_order",
                  # facet.by = c("sampling", "num_tip"),
                  facet.by = "sampling",
                  color = "log2_node_sampled",
                  alpha = 0.05,
                  size = 0.15,
                  # xlab = "True commitment ratio",
                  # ylab = "Inferred commitment ratio"
                  xlab = "",
                  ylab = ""
                  ) +
    geom_abline() +
    geom_smooth(aes(x =  node_split_order, y = gr_node_split_order),
                color = sampling_col[9],
                method = "loess",
                se = F,
                span = 1,
                data = filter(eval_tb_plot, suff_sampled == "TRUE")) +
    geom_smooth(aes(x =  node_split_order, y = gr_node_split_order),
                color = sampling_col[1],
                method = "loess",
                span = 1,
                se = F,
                data = filter(eval_tb_plot, suff_sampled == "FALSE")) +
    # stat_cor(aes(group = suff_sampled,
    #              label =  ..r.label..),
    #          method = "spearman",
    #          cor.coef.name = "rho",
    #          digits = 3) +
    scale_colour_gradientn(limits = c(-4, 0), colors = sampling_col[c(2, 5, 8)],
                           breaks = (-4):(-1), na.value = sampling_col[1]) +
        theme(text = element_text(size = 10), legend.position = "none", strip.background = element_blank())) %>%
  # push_png(file_name = "node_split_est_bytip", res = 1200, w = 6.5, h = 2.8, dir = plot_dir)
  push_png(file_name = "node_split_est_merge_nolab", res = 1200, w = 4.0, h = 2.0, dir = plot_dir)


# summary plots
# (eval_tb_plot %>%
#                 filter(suff_sampled) %>%
#                 ggline(x = "bsum",
#                        y = "gr_node_size_logfc",
#                        facet.by = c("sampling", "num_tip"),
#                        add = "mean_se",
#                        scales = "free_x") +
#                 geom_hline(yintercept = 0)) #%>%
        # push_pdf(file_name = "node_size", w = 7.5, h = 2.5, dir = "../LTModelPlots/temp_panel_eval/")

# range of node sizes
# eval_tb_plot %>%
#         filter(suff_sampled) %>%
#         gghistogram(x = "log2_node_size", fill = "orange")

(eval_tb_plot %>%
                filter(suff_sampled) %>%
                ggline(x = "bsum",
                       y = "gr_time_trans_error",
                       facet.by = c("sampling", "num_tip"),
                       add = "mean_se",
                       scales = "free_x") +
                geom_hline(yintercept = 0)) #%>%
        #push_pdf(file_name = "time_error", w = 7.5, h = 2.5, dir = "../LTModelPlots/temp_panel_eval/")

(eval_tb_plot %>% group_by(j) %>%
                summarise(bsum = bsum[1],
                          sampling = sampling[1],
                          num_tip = num_tip[1],
                          frac_suff_sampled = mean(suff_sampled),
                          frac_resolved = mean(is_resolved)) %>%
                ggline(x = "bsum",  y = "frac_resolved", facet.by = c("sampling", "num_tip"), scales = "free_x",
                       add = "mean_se") +
                theme(text = element_text(size = 10),
                      axis.text.x = element_text(size = 9, angle = 65))
) %>%
        push_pdf(file_name = "frac_resolved_v1", w = 6.5, h = 3., dir = plot_dir)



(eval_tb_plot %>%
    group_by(j) %>%
    summarise(bsum = bsum[1],
              sampling = sampling[1],
              num_tip = num_tip[1],
              frac_suff_sampled = sum(suff_sampled, na.rm = T) / (num_tip - 1),
              frac_resolved = mean(is_resolved)) %>%
    ggline(x = "bsum",  y = "frac_suff_sampled", facet.by = c("sampling", "num_tip"), scales = "free_x",
           add = "mean_se", ylim = c(0, 1)) +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(size = 9, angle = 65))
) %>%
  push_pdf(file_name = "frac_suff_sampled_v1", w = 6.5, h = 3., dir = plot_dir)

(eval_tb_plot %>%
    group_by(j) %>%
    summarise(bsum = bsum[1],
              sampling = sampling[1],
              num_tip = num_tip[1],
              frac_resolved = mean(is_resolved),
              frac_suff_sampled = mean(suff_sampled)) %>%
    ggline(x = "bsum",
           y = "frac_suff_sampled",
           facet.by = c("sampling", "num_tip"),
           scales = "free_x",
           add = "mean_se", ylim = c(0, 1)) +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(size = 9, angle = 65))
  ) #%>%
   #     push_pdf(file_name = "frac_suff_sampled", w = 6.5, h = 3., dir = plot_dir)
# ggscatter(eval_tb_plot, x = "mean_cosar_error", y = "is_resolved")+ geom_smooth(method = "lm")


hist(log2(eval_tb_plot$node_size_sampled / eval_tb_plot$node_size), breaks = 100)
abline(v = -1)

ggscatter(eval_tb_plot, x = "")

qplot(eval_tb_plot$gr_frac_short_edge,
      eval_tb_plot$log2_node_sampled) + geom_hline(yintercept = -1)

node_collect_size = map_dbl(eval_tb_plot$node_size_collect, function(x) {
  if (is.null(x)){
    return(NA)
  } else  {
    return(x)
  }
})
qplot((eval_tb_plot$gr_node_size_in / node_collect_size),
      eval_tb_plot$gr_time_trans_error)

# histogram for node assignment accuracy per experiment
node_assign_raw = bind_rows(evaluate_node_assign(exp_params, all_graphs, tr_col = "tr", gr_col = "gr"))
saveRDS(node_assign_raw, file = "./intermediate_data/panel/exp_node_assign.rds")
node_assign_accracy = node_assign_raw %>%
        group_by(j) %>%
        summarise(accuracy = mean(correct)) %>%
        gghistogram(x = "accuracy", fill = "#1572A1", color = NA)

# further evaluations of node assignments and sampling fraction (unused)
# node_assign_tb = bind_rows(evaluate_node_assign(exp_params, all_graphs, tr_col = "tr", gr_col = "gr"))
# node_assign_tb = mutate(node_assign_tb, logit_sampling_frac = log(sampling_frac / (1 - sampling_frac)))
# node_assign_tb = mutate(node_assign_tb,
#                         logit_sampling_frac = log(sampling_frac / (1 - sampling_frac)),
#                         logit_frac_correct = log(frac_correct / (1 - frac_correct)))
# node_assign_tb %>% group_by(j) %>% summarize()
#
# cor(node_assign_tb$frac_correct, node_assign_tb$logit_sampling_frac, use = "pairwise.complete.obs")
# ggscatter(node_assign_tb, x = "logit_sampling_frac", y = "logit_frac_correct") + geom_smooth()


#### Second Section ####
# further panel evaluations involving Phylotime
load(paste0(output_dir, "exp_gr3_eval.rda"))
set.seed(73)
eval_sum = bind_rows(
        eval_tb_wrap(eval_tb_gr, exp_params) %>% group_by(sampling, num_tip, suff_sampled) %>%
                filter(!is.na(gr_time_trans) & !is.na(node_time)) %>%
                summarise(cc = cor(node_time, gr_time_trans, method = "pearson"),
                          cs = cor(node_time, gr_time_trans, method = "spearman"),
                          rmse = sqrt(mean((node_time - gr_time_trans)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x truth"),
        eval_tb_wrap(eval_tb_gr3, exp_params) %>% group_by(sampling, num_tip, suff_sampled) %>%
                filter(!is.na(gr_time_trans) & !is.na(node_time)) %>%
                summarise(cc = cor(node_time, gr_time_trans, method = "pearson"),
                          cs = cor(node_time, gr_time_trans, method = "spearman"),
                          rmse = sqrt(mean((node_time - gr_time_trans)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x phylidite"))
eval_sum$method = factor(eval_sum$method, c("ice_fase x phylidite", "ice_fase x truth"))

eval_sum_merge = bind_rows(
        eval_tb_wrap(eval_tb_gr, exp_params) %>% group_by(sampling, suff_sampled) %>%
                filter(!is.na(gr_time_trans) & !is.na(node_time)) %>%
                summarise(cc = cor(node_time, gr_time_trans, method = "pearson"),
                          cs = cor(node_time, gr_time_trans, method = "spearman"),
                          rmse = sqrt(mean((node_time - gr_time_trans)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x truth"),
        eval_tb_wrap(eval_tb_gr3, exp_params) %>% group_by(sampling, suff_sampled) %>%
                filter(!is.na(gr_time_trans) & !is.na(node_time)) %>%
                summarise(cc = cor(node_time, gr_time_trans, method = "pearson"),
                          cs = cor(node_time, gr_time_trans, method = "spearman"),
                          rmse = sqrt(mean((node_time - gr_time_trans)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x phylidite"))
eval_sum_merge
eval_sum_merge %>% group_by(sampling, suff_sampled) %>% summarise(diff = cs[1] - cs[2])

g1 = ggbarplot(eval_sum, x = "suff_sampled", y = "cs", fill = "method",
               color = NA,
               position = position_dodge(),
               facet.by = c("sampling", "num_tip"),
               ylab = "Correlation with Truth") +
        scale_fill_manual(values = method_col) +
        theme(text = element_text(size = 10),
              legend.position = "right")
g1m = ggbarplot(eval_sum_merge, x = "suff_sampled", y = "cs", fill = "method",
                color = NA,
                position = position_dodge(),
                facet.by = c("sampling"),
                title = "Commitment time",
                xlab = ">= 25% Sampled",
                ylab = "Correlation with Truth") +
        scale_fill_manual(values = method_col) +
        theme(text = element_text(size = 10),
              legend.position = "right",
              strip.background = element_blank())

eval_sum = bind_rows(
        eval_tb_wrap(eval_tb_gr) %>% group_by(sampling, num_tip, suff_sampled) %>%
                filter(!is.na(log2_gr_node_size_in) & !is.na(log2_node_size)) %>%
                summarise(cc = cor(log2_node_size, log2_gr_node_size_in, method = "pearson"),
                          cs = cor(log2_node_size, log2_gr_node_size_in, method = "spearman"),
                          ma_log2fc = mean(abs(log2_node_size - log2_gr_node_size_in))),
        total = n()) %>%
        mutate(method = "ice_fase x truth"),
eval_tb_wrap(eval_tb_gr3, exp_params) %>% group_by(sampling, num_tip, suff_sampled) %>%
        filter(!is.na(log2_gr_node_size_in) & !is.na(log2_node_size)) %>%
        summarise(cc = cor(log2_node_size, log2_gr_node_size_in, method = "pearson"),
                  cs = cor(log2_node_size, log2_gr_node_size_in, method = "spearman"),
                  ma_log2fc = mean(abs(log2_node_size - log2_gr_node_size_in)),
                  total = n()) %>%
        mutate(method = "ice_fase x phylidite")
eval_sum$method = factor(eval_sum$method, c("ice_fase x phylidite", "ice_fase x truth"))

eval_sum_merge = bind_rows(
        eval_tb_wrap(eval_tb_gr, exp_params) %>% group_by(sampling, suff_sampled) %>%
                filter(!is.na(log2_gr_node_size_in) & !is.na(log2_node_size)) %>%
                summarise(cc = cor(log2_node_size, log2_gr_node_size_in, method = "pearson"),
                          cs = cor(log2_node_size, log2_gr_node_size_in, method = "spearman"),
                          rmse = sqrt(mean((log2_node_size - log2_gr_node_size_in)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x truth"),
        eval_tb_wrap(eval_tb_gr3, exp_params) %>% group_by(sampling, suff_sampled) %>%
                filter(!is.na(log2_gr_node_size_in) & !is.na(log2_node_size)) %>%
                summarise(cc = cor(log2_node_size, log2_gr_node_size_in, method = "pearson"),
                          cs = cor(log2_node_size, log2_gr_node_size_in, method = "spearman"),
                          rmse = sqrt(mean((log2_node_size - log2_gr_node_size_in)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x phylidite"))
eval_sum_merge
eval_sum_merge %>% group_by(sampling, suff_sampled) %>% summarise(diff = cs[1] - cs[2])

g2 = ggbarplot(eval_sum, x = "suff_sampled", y = "cs", fill = "method",
               color = NA,
               position = position_dodge(),
               facet.by = c("sampling", "num_tip"),
               xlab = ">= 25% Sampled",
               ylab = "Correlation with Truth") +
        scale_fill_manual(values = method_col) +
        theme(text = element_text(size = 10),
              legend.position = "right")

g2m = ggbarplot(eval_sum_merge, x = "suff_sampled", y = "cs", fill = "method",
                color = NA,
                position = position_dodge(),
                facet.by = c("sampling"),
                title = "Progenitor population size",
                xlab = ">= 25% Sampled",
                ylab = "Correlation with Truth") +
        scale_fill_manual(values = method_col) +
        theme(text = element_text(size = 10),
              legend.position = "right",
              strip.background = element_blank())

eval_sum = bind_rows(
        eval_tb_wrap(eval_tb_gr, exp_params) %>% group_by(sampling, num_tip, suff_sampled) %>%
                filter(!is.na(gr_node_split_order) & !is.na(node_split_order)) %>%
                summarise(cc = cor(node_split_order, gr_node_split_order)) %>% mutate(method = "ice_fase"),
        eval_tb_wrap(eval_tb_gr3, exp_params) %>% group_by(sampling, num_tip, suff_sampled) %>%
                filter(!is.na(gr_node_split_order) & !is.na(node_split_order)) %>%
                summarise(cc = cor(node_split_order, gr_node_split_order)) %>% mutate(method = "ice_fase_phylidite"))
eval_sum$method = factor(eval_sum$method, c("ice_fase_phylidite", "ice_fase"))

eval_sum = bind_rows(
        eval_tb_wrap(eval_tb_gr, exp_params) %>% group_by(sampling, num_tip, suff_sampled) %>%
                filter(!is.na(gr_node_split_order) & !is.na(node_split_order)) %>%
                summarise(cc = cor(node_split_order, gr_node_split_order, method = "pearson"),
                          cs = cor(node_split_order, gr_node_split_order, method = "spearman"),
                          rmse = sqrt(mean((node_split_order - gr_node_split_order)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x truth"),
        eval_tb_wrap(eval_tb_gr3, exp_params) %>% group_by(sampling, num_tip, suff_sampled) %>%
                filter(!is.na(gr_node_split_order) & !is.na(node_split_order)) %>%
                summarise(cc = cor(node_split_order, gr_node_split_order, method = "pearson"),
                          cs = cor(node_split_order, gr_node_split_order, method = "spearman"),
                          rmse = sqrt(mean((node_split_order - gr_node_split_order)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x phylidite"))
eval_sum$method = factor(eval_sum$method, c("ice_fase x phylidite", "ice_fase x truth"))

eval_sum_merge = bind_rows(
        eval_tb_wrap(eval_tb_gr, exp_params) %>% group_by(sampling, suff_sampled) %>%
                filter(!is.na(gr_node_split_order) & !is.na(node_split_order)) %>%
                summarise(cc = cor(node_split_order, gr_node_split_order, method = "pearson"),
                          cs = cor(node_split_order, gr_node_split_order, method = "spearman"),
                          rmse = sqrt(mean((node_split_order - gr_node_split_order)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x truth"),
        eval_tb_wrap(eval_tb_gr3, exp_params) %>% group_by(sampling, suff_sampled) %>%
                filter(!is.na(gr_node_split_order) & !is.na(node_split_order)) %>%
                summarise(cc = cor(node_split_order, gr_node_split_order, method = "pearson"),
                          cs = cor(node_split_order, gr_node_split_order, method = "spearman"),
                          rmse = sqrt(mean((node_split_order - gr_node_split_order)^2)),
                          total = n()) %>%
                mutate(method = "ice_fase x phylidite"))
eval_sum_merge
eval_sum_merge %>% group_by(sampling, suff_sampled) %>% summarise(diff = cs[1] - cs[2])

g3 = ggbarplot(eval_sum, x = "suff_sampled", y = "cs", fill = "method",
               color = NA,
               position = position_dodge(),
               facet.by = c("sampling", "num_tip"),
               xlab = ">= 25% Sampled",
               ylab = "Correlation with Truth") + scale_fill_manual(values = method_col) +
        theme(text = element_text(size = 10),
              legend.position = "right")

g3m = ggbarplot(eval_sum_merge, x = "suff_sampled", y = "cc", fill = "method",
                color = NA,
                position = position_dodge(),
                facet.by = c("sampling"),
                title = "Commitment bias",
                xlab = ">= 25% Sampled",
                ylab = "Correlation with Truth") + scale_fill_manual(values = method_col) +
        theme(text = element_text(size = 10),
              legend.position = "right",
              strip.background = element_blank())

method_col = c("sps x truth" = "#11468F",
               "sps" = "#11468F",
               "ice_fase x hamming" = "#FC28FB",
               "ice_fase_hamming" = "#FC28FB",
               "ice_fase x phylidite" = "#FFC900",
               "ice_fase_phylidite" = "#FFC900",
               "ice_fase x truth" = "#DA1212",
               "ice_fase" = "#DA1212")

(g1m / g2m / g3m) %>%
        push_pdf(file_name = "eval_pp_params_merge_cov",
                 width = 6.5, height = 6, ps = 10,
                 dir = plot_dir)
#### End Second Section ####








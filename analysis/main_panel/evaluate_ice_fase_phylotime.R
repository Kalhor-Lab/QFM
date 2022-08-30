library(qfm)
plot_dir = "./plots/panel_mod2_v1/"
output_dir = "./intermediate_data/panel_mod2_v1/"
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))

exp_params = readRDS(paste0(output_dir, "exp_data_5rep_mut_all.rds"))
exp_params$exp_id = 1:nrow(exp_params)
data_dir = paste0(output_dir, "exp_data_5rep_mut_all/")
# results with true phylogeny
res_dir = paste0(data_dir, "res/")
exp_params = gather_ice_fase(exp_params, tree_panel, res_dir)

# evaluate node assignment accuracy
tr_node_assign_long = do.call(c, map(exp_params$res, function(res) {
        res$tr_node_assign
}))
tr_node_assign_long = tr_node_assign_long[sample(length(tr_node_assign_long), 1e5)]
1 - mean(tr_node_assign == get_type_from_id(names(tr_node_assign)))

mean(map_dbl(exp_params$res[exp_params$i_sim == 1], function(res) {
        tr_node_assign = res$tr_node_assign
        mean(tr_node_assign == get_type_from_id(names(tr_node_assign)))
}))
# end evalaute node assignment accuracy

# load experiment data
load_exp_tr <- function(exp_params, rep_dir) {
        data_list = map(1:nrow(exp_params), function(i) {
                message(i)
                out_file = paste0(rep_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
                try({
                        out = readRDS(out_file)
                        out
                })
        })
        # exp_params$sc = map(data_list, "sc")
        exp_params$tr = map(data_list, "tr")
        exp_params$true_sampled_sizes = map(data_list, "true_sampled_sizes")
        exp_params
}
exp_params = load_exp_tr(exp_params, data_dir)

# evaluations with true phylogeny
eval_tb_list_true = readRDS("./intermediate_data/panel_mod2_v1/5rep_suff_eval_tb.rds")

# exp_params = generate_control(exp_params, tree_panel)
expm1_trans <-  function() trans_new("expm1", "expm1", "log1p")
method_col_v1 = c("sps" = "#377eb8",
                  "ice_fase" = "#e41a1c",
                  "random_control" = "#949494")
exp_params$kc0 = readRDS("./intermediate_data/panel_mod2_v1/5rep_eval_kc0.rds")
exp_params$kc0_control = map_dbl(exp_params$gr_control_eval, "kc0")
write_csv(select(exp_params, -c(res, gr_control_eval)),
          file = "./intermediate_data/panel_mod2_v1/5rep_exp_stats.csv")

#### SPS part ####
run_sps <- function(exp_params, data_dir, i, res_dir) {
        data_obj = readRDS(paste0(data_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds"))
        tr_r = data_obj$tr
        sps_mat = shared_progenitor_score(tr_r,
                                          sort(unique(get_type_from_id(tr_r$tip.label))))
        dmat = 1 - sps_mat/max(sps_mat)
        gr = phangorn::upgma(as.dist(dmat))
        total_depth = max(node.depth.edgelength(gr))
        gr$edge.length = gr$edge.length / total_depth * total_time #TODO: what to use for sps as max dist?
        gr = name_nodes(gr)
        out_file = paste0(res_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        saveRDS(gr, file = out_file)
}
sps_dir = paste0(data_dir, "sps/")
# dir.create(sps_dir)
future_walk(1:nrow(exp_params), function(i) {
        run_sps(exp_params, data_dir, i, sps_dir)
}, .progress = T)

exp_params$sps_gr = map(1:nrow(exp_params), function(i) {
        message(i)
        out_file = paste0(sps_dir, stringr::str_pad(i, width = 4, pad = "0"), ".rds")
        readRDS(out_file)
})
exp_params$sps_gr_eval = future_map2(exp_params$sps_gr, exp_params$big_graph_id,
                                     function(x, i) {
                                             if (is.null(x)) {
                                                     return(NULL)
                                             }
                                             evalute_gr(x, tree_panel$gr[[i]])
                                     }, .progress = T, .options = furrr_options(seed = T))
# END SPS part

# exp_params = readRDS("./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all_proc.rds")
# true_kc0 = exp_params$kc0
# saveRDS(true_kc0, file = "./intermediate_data/panel_mod2_v1/5rep_eval_kc0.rds")

exp_params = read_csv("./intermediate_data/panel_mod2_v1/5rep_exp_stats.csv")
# collecting results from mixed 50site
exp_params$tr3_kc0 = readRDS("./intermediate_data/panel_mod2_v1/5rep_50site_eval_kc0.rds")
exp_params$num_tip = tree_panel$num_tip[exp_params$big_graph_id]
exp_params$bsum = tree_panel$bsum[exp_params$big_graph_id]

exp_params %>% group_by(num_tip, sampling) %>%
        summarize(mean_kc0 = mean(kc0),
                  mean_kc0_tr3 = mean(tr3_kc0))

method_col = c("sps x truth" = "#11468F",
               "sps" = "#11468F",
               # "ice_fase x hamming" = "#FC28FB",
               # "ice_fase_hamming" = "#FC28FB",
               "ice_fase_phylotime" = "#6870a4",
               "ice_fase_truth" = "#99b9e1",
               "random" = "#d3d3d3")
expm1_trans <-  function() trans_new("expm1", "expm1", "log1p")
(exp_params %>% transmute(sampling = sampling,
                          num_tip = num_tip,
                          bsum = bsum,
                          random = kc0_control,
                          ice_fase_phylotime = tr3_kc0,
                          ice_fase_truth = kc0
                          ) %>%
                gather(key = "method", value = "kc0", -c("sampling", "num_tip", "bsum")) %>%
                ggscatter(x = "bsum", y = "kc0", color = "method",
                          facet.by = c("sampling", "num_tip"),
                          ylab = "KC0",
                          xlab = "Colless index",
                          xlim = c(0, NA),
                          scales = "free",
                          size = 0.05,
                          alpha  = 0.1) +
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
        push_pdf(file_name = "eval_res_kc0_wctrl", width = 5.5, height = 2.5, ps = 12, dir = plot_dir)

eval_tb_list = readRDS("./intermediate_data/panel_mod2_v1/5rep_50site_suff_eval_tb.rds")
eval_tb_list_true = readRDS("./intermediate_data/panel_mod2_v1/5rep_suff_eval_tb.rds")

# evaluate average correlations using true vs phylotime phylogeny
combine_type <- function(type) {
        eval_sum = bind_rows(map(eval_tb_list, type))
        eval_sum = dplyr::rename(eval_sum) %>% left_join(exp_params)
        eval_sum_true = bind_rows(map(eval_tb_list_true, type))
        eval_sum_true = dplyr::rename(eval_sum_true) %>% left_join(exp_params)
        bind_rows(mutate(eval_sum_true, method = "ice_fase_truth"),
                  mutate(eval_sum, method = "ice_fase_phylotime")) %>%
                mutate(type = type)
}
temp_sum = bind_rows(combine_type("time"),
                     combine_type("size"),
                     combine_type("split"))
temp_sum$cs[is.na(temp_sum$cs)] = 0

temp_sum %>%
        group_by(method, type) %>%
        filter(suff_sampled) %>%
        summarize(cs = mean(cs)) %>%
        group_by(type, sampling) %>%
        summarize(pct = cs[2] - cs[1])

(temp_sum %>%
        group_by(sampling, method, type, suff_sampled) %>%
        ggbarplot(x = "type", y = "cs",
                  color = NA,
                  fill = "method",
                  add = "mean_se",
                  position = position_dodge(),
                  facet.by = c("sampling", "suff_sampled")) +
        scale_fill_manual(values = method_col)) %>%
        push_pdf(file_name = "5rep_eval_param_cor_v1", width = 4, height = 3, ps = 12, dir = plot_dir)
# end evaluate average correlations

# results from different sets of hgRNAs (slow, mid, fast, 25, 50, 100)
mut_flat_all = readRDS("./intermediate_data/panel_mod2_v1/mut_flat_all.rds")
fast_ind = mut_flat_all$class == "fast"
slow_ind = mut_flat_all$class == "slow"
mut_flat_all$class[fast_ind] = "slow"
mut_flat_all$class[slow_ind] = "fast"
eval_gr_list = purrr::reduce(map(1:20, function(i) {
        readRDS(paste0("./intermediate_data/panel_mod2_v1/mut_flat_gr_eval_part", i, ".rds"))
}), c)
mut_flat_all$kc0 = map_dbl(eval_gr_list, "kc0")

all_mut_indices = readRDS(paste0(output_dir, "all_mut_indices.rds"))
fast_ind = all_mut_indices$class == "fast"
slow_ind = all_mut_indices$class == "slow"
all_mut_indices$class[fast_ind] = "slow"
all_mut_indices$class[slow_ind] = "fast"

eval_tb_list = readRDS("./intermediate_data/panel_mod2_v1/mut_flat_eval_gr_tb.rds")
eval_tb_list_true = readRDS("./intermediate_data/panel_mod2_v1/5rep_eval_gr_tb.rds")
# end results different sets of hgRNAs

plot_sum_template <- function(type = c("time", "size", "split")) {
        eval_sum = bind_rows(map(eval_tb_list, type))
        eval_sum_true = bind_rows(map(eval_tb_list_true, type))
        eval_sum = dplyr::rename(eval_sum, mut_id = exp_id) %>% left_join(mut_flat_all)
        eval_sum_true = left_join(eval_sum_true, exp_params, by = "exp_id")
        bind_rows(eval_sum,
                  mutate(eval_sum_true, class = "fast"),
                  mutate(eval_sum_true, class = "mid"),
                  mutate(eval_sum_true, class = "slow")) %>%
                filter(class != "somatic") %>%
                ggline(x = "diversity", y = "rmse", add = "mean_se", color = "class")
}

g0 = bind_rows(mut_flat_all[mut_flat_all$class != "somatic", ],
          select(exp_params, -res) %>% mutate(class = "fast"),
          select(exp_params, -res) %>% mutate(class = "mid"),
          select(exp_params, -res) %>% mutate(class = "slow")) %>%
        ggline(x = "num_element", y = "kc0", color = "class", add = "mean_se")

g0 +
        plot_sum_template("time") +
        plot_sum_template("size") +
        plot_sum_template("split") +
        plot_layout(nrow = 1, guide = "collect") &
        theme(legend.position = "bottom")








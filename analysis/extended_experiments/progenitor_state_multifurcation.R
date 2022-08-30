library(qfm)
library(furrr)
plan(multisession, workers = 12)
plot_dir = "./plots/panel_mod2_v1_examples/"
#' delete a node and connect its daughter states to parent
delete_state <- function(type_graph, state = "14") {
        assertthat::assert_that(state %in% type_graph$node_id)
        assertthat::assert_that(state != type_graph$root_id)
        parent_node = type_graph$edges$in_node[type_graph$edges$out_node == state]

        parent_merge_old = type_graph$merge[[parent_node]]
        state_merge_old = type_graph$merge[[state]]
        parent_diff_mode_probs_old = type_graph$diff_mode_probs[[parent_node]]
        state_diff_mode_probs_old = type_graph$diff_mode_probs[[state]]
        parent_state_ind = parent_merge_old == state

        assertthat::assert_that(parent_diff_mode_probs_old[3] == 0)
        parent_merge_new = c(parent_merge_old[!parent_state_ind], state_merge_old)
        parent_diff_mode_probs_new = c(parent_diff_mode_probs_old[1:2][!parent_state_ind],
                                       parent_diff_mode_probs_old[1:2][parent_state_ind] * state_diff_mode_probs_old[1:2],
                                       rep(0, 3))
        type_graph$merge[[parent_node]] = parent_merge_new
        type_graph$diff_mode_probs[[parent_node]] = parent_diff_mode_probs_new

        type_graph$node_id = setdiff(type_graph$node_id, state)
        type_graph$all_id = c(type_graph$node_id, type_graph$tip_id)
        type_graph$merge = type_graph$merge[type_graph$node_id]
        type_graph$diff_time = type_graph$diff_time[type_graph$all_id]
        type_graph$diff_mode_probs = type_graph$diff_mode_probs[type_graph$node_id]
        type_graph$lifetime = type_graph$lifetime[type_graph$all_id]
        type_graph$prob_loss = type_graph$prob_loss[type_graph$all_id]
        type_graph$prob_pm = type_graph$prob_pm[type_graph$all_id]
        type_graph = generate_edge_df_mod2(type_graph)
        type_graph
}
tree_panel = readRDS(paste0("./intermediate_data/panel_mod2_v1/", "tree_panel.rds"))
type_graph = tree_panel$type_graph[[22]]

set.seed(73)
gr_col_vec = gr_color_v1(type_graph)
plot_type_graph_clean_mod2(type_graph, node_col_mapper = function(x) gr_col_vec[x], show_node_text = T) %>%
        push_pdf("g_multifurc_original.pdf", w = 2.5, h = 2.5, dir = plot_dir)

type_graph_mod1 = delete_state(type_graph, "8")
plot_type_graph_clean_mod2(type_graph_mod1, node_col_mapper = function(x) gr_col_vec[x], show_node_text = T) %>%
        push_pdf("g_multifurc_mod1.pdf", w = 2.5, h = 2.5, dir = plot_dir)

# type_graph_mod2 = delete_state(type_graph, "7")
# plot_type_graph_clean_mod2(type_graph_mod2, node_col_mapper = function(x) gr_col_vec[x], show_node_text = T)
# type_graph_mod3 = delete_state(type_graph, "4")
# plot_type_graph_clean_mod2(type_graph_mod3, node_col_mapper = function(x) gr_col_vec[x], show_node_text = T)

# tree_panel_mod = tibble(type_graph = list(type_graph_mod1,
#                                           type_graph_mod2,
#                                           type_graph_mod3))
# tree_panel_mod$gr = map(tree_panel_mod$type_graph, as.phylo_mod2.type_graph)
# exp_params = expand.grid(big_graph_id = 1:3,
#                          sample_size = c(100),
#                          num_element = c(50),
#                          sampling = c("fixed"),
#                          i_sim = 1:50) %>% as_tibble()

tree_panel_mod = tibble(type_graph = list(type_graph,
                                          type_graph_mod1))
tree_panel_mod$gr = map(tree_panel_mod$type_graph, as.phylo_mod2.type_graph)
exp_params = expand.grid(big_graph_id = 1:2,
                         sample_size = c(100),
                         num_element = c(50),
                         sampling = c("fixed"),
                         i_sim = 1:100) %>% as_tibble()
rep_dir = "./intermediate_data/panel_mod2_v1/multifurc_v2/"
# dir.create(rep_dir)
data_dir = paste0(rep_dir, "data/")
# dir.create(data_dir)
# saveRDS(tree_panel_mod, file = paste0(rep_dir, "tree_panel_mod.rds"))
# saveRDS(exp_params, file = paste0(rep_dir, "exp_params.rds"))
mut_p = readRDS("./intermediate_data/panel_mod2_v1/mut_p_marc1.rds")
future_walk(1:nrow(exp_params), function(i) {
        run_experiment(exp_params, tree_panel_mod, i, data_dir, mut_p = mut_p)
})
# exp_params = load_exp_data(exp_params, data_dir)

exp_params = readRDS(paste0(rep_dir, "exp_params.rds"))
tree_panel_mod = readRDS(paste0(rep_dir, "tree_panel_mod.rds"))

# dir.create(paste0(rep_dir, "tr3_res/"))
walk_ice_fase_data_tr3(exp_params,
                       input_dir = paste0(rep_dir, "tr3/"),
                       tree_panel_mod,
                       out_dir = paste0(rep_dir, "tr3_res/"))

exp_params = gather_ice_fase(exp_params, tree_panel_mod, out_dir = paste0(rep_dir, "tr3_res/"))
plot(exp_params$res[[1]]$gr)

# report:
# is the sub-topology topology resolved?
exp_params = exp_params[exp_params$big_graph_id == 1, ]

# Question #1: is the parent state potency correct
# Question #2: if the cloest immedaite descendent gets merged, is the trifurcation correct?
parent_node_potent = as.character((-1):(-5))
split_chr = list(as.character((-1)),
                 as.character((-2):(-3)),
                 as.character((-4):(-5))) %>%
        map_chr(function(x) paste0(sort(x), collapse = "_"))


# Here dummy is a node that is a pseudo bifurcation
node_edge_len_all = bind_rows(map(1:nrow(exp_params), function(i) {
        res = exp_params$res[[i]]
        gr_tr_data = res$gr_tr_data
        node_ind = which(map_lgl(gr_tr_data$gr_tip_list, function(x) {
                satet_vec = parent_node_potent
                all(satet_vec %in% x) & (length(x) == length(satet_vec))
        }))
        if (length(node_ind) == 0) {
                return(NULL)
        }
        gr_edges_time_tb = res$gr_tr_data$gr_edges_tb
        gr_edges_time_tb$from_time = res$gr_trans_time[gr_edges_time_tb$from]
        gr_edges_time_tb$to_time = res$gr_trans_time[gr_edges_time_tb$to]
        gr_edges_time_tb$length = gr_edges_time_tb$to_time - gr_edges_time_tb$from_time

        gr_node_mapped = get_node_mapping_mod2(gr_tr_data, tree_panel_mod$type_graph[[1]])
        gr_edges_node = gr_edges_time_tb[!gr_edges_time_tb$to %in% gr_tr_data$gr$tip.label, ]
        gr_edges_node$from_mapped = gr_node_mapped[gr_edges_node$from]
        gr_edges_node$to_mapped = gr_node_mapped[gr_edges_node$to]
        gr_edges_node$from_pot = map_chr(gr_tr_data$gr_tip_list[gr_edges_node$from], function(x) {
                paste(sort(x), collapse = "_")
        })
        gr_edges_node$to_pot = map_chr(gr_tr_data$gr_tip_list[gr_edges_node$to], function(x) {
                paste(sort(x), collapse = "_")
        })
        # take the shorter edge downstream
        gr_edges_sel = gr_edges_time_tb[gr_edges_time_tb$from == names(node_ind), ]
        node_dummy = gr_edges_sel$to[which.min(gr_edges_sel$length)]
        node_splits = c(gr_edges_sel$to[which.max(gr_edges_sel$length)],
                        gr_tr_data$gr_dd[[node_dummy]])
        node_split_tips = gr_tr_data$gr_tip_list[node_splits]
        if (all(split_chr %in% map_chr(node_split_tips, function(x) paste0(sort(x), collapse = "_")))) {
                # print(node_split_tips)
                gr_edges_node$dum_ind = gr_edges_node$to == node_dummy
                gr_edges_node$cond = exp_params$big_graph_id[i]
                return(gr_edges_node)
        } else {
                return(NULL)
        }
}))
node_edge_len_all = node_edge_len_all %>% mutate(to_merge = length < 0.45)
node_edge_len_all %>% group_by(cond, dum_ind) %>%
        summarise(total_merged = mean(to_merge))

node_edge_len_all[node_edge_len_all$to_merge & !node_edge_len_all$dum_ind, ] %>%
        count(from_mapped, from_pot, to_mapped, to_pot) %>% arrange(desc(n))
node_edge_len_all[node_edge_len_all$to_merge & !node_edge_len_all$dum_ind, ] %>%
        count(from_pot) %>% arrange(desc(n))

node_edge_len_all$

node_edge_len_all %>% ggboxplot(x = "dum_ind", y = "length", xlab = "Pseudo-edge",facet.by = "cond") %>%
        push_pdf("multifurc_pseudo_edge_len", w = 2.5, h = 2.5, dir = plot_dir)

node_edge_len_all %>% group_by(cond, dum_ind) %>%
        summarise(mean_merge = mean(length <= 0.45))

node_edge_len_all %>% group_by(dum_ind) %>%
        do(w = wilcox.test(length ~ cond, data=., paired=FALSE)) %>%
        summarise(dum_ind, Wilcox = w$p.value)

node_edge_len_all %>% group_by(cond) %>%
        do(w = wilcox.test(length ~ dum_ind, data=., paired=FALSE)) %>%
        summarise(cond, Wilcox = w$p.value)


res = exp_params$res[[1]]
gr_col_alt = gr_color_v1(res$gr)
gr_col_alt[res$gr$tip.label] = gr_col_vec[res$gr$tip.label]
gr_tip_time = rep(type_graph$target_time, length(res$gr$tip.label))
gr_node_lab = map_gr_label(c(res$gr$node.label, res$gr$tip.label))
plot_gr_clean(res$gr,
              gr_node_time = c(res$gr_trans_time, gr_tip_time),
              total_time = type_graph$target_time,
              type_col = gr_col_alt,
              show_node_label = T,
              node_label = gr_node_lab
              ) %>%
        push_pdf("example_gr_multifurc_original.pdf", w = 2.5, h = 2.5, dir = plot_dir)

res_correct = correct_gr_edge_len(res)
res_correct$gr = di2multi(res_correct$gr, tol = 0.5)

res_update = ice_fase_mutlifurc(res$tr,
                                res$sc_celltypes,
                                total_time = res$total_time,
                                root_time = res$root_time,
                                gr = res_correct$gr)
(plot_gr_clean(res_update$gr,
              gr_node_time = c(res_update$gr_trans_time, gr_tip_time),
              total_time = type_graph$target_time,
              type_col = gr_col_alt,
              show_node_label = T,
              node_label = gr_node_lab
) + scale_x_reverse()) %>%
        push_pdf("example_gr_multifurc_original_merged.pdf", w = 2.5, h = 2.5, dir = plot_dir)


res2 = exp_params$res[[10]]
res1_node_mapped = get_node_mapping_mod2(res$gr_tr_data, tree_panel_mod$type_graph[[1]])
res2_node_mapped = get_node_mapping_mod2(res2$gr_tr_data, tree_panel_mod$type_graph[[1]])
res1_node_mapped_ord = names(res1_node_mapped)[match(type_graph$node_id, res1_node_mapped)]
res2_node_mapped_ord = names(res2_node_mapped)[match(type_graph$node_id, res2_node_mapped)]
res2_node_mapped_ord
gr_col_alt_mapped = gr_col_alt[res1_node_mapped_ord]
names(gr_col_alt_mapped) = res2_node_mapped_ord
gr_col_alt_mapped = c(gr_col_alt_mapped, gr_col_vec[res$gr$tip.label])
gr_node_lab_mapped = gr_node_lab[res1_node_mapped_ord]
names(gr_node_lab_mapped) = res2_node_mapped_ord
gr_node_lab_mapped = c(gr_node_lab_mapped, gr_node_lab[res$gr$tip.label])

# gr_col_alt_mapped[is.na(names(gr_col_alt_mapped))] = c("#e7298a", "#ce1256")
# names(gr_col_alt_mapped)[is.na(names(gr_col_alt_mapped))] = setdiff(res2$gr$node.label, res2_node_mapped_ord)


res2_correct = correct_gr_edge_len(res2)
res2_correct$gr = di2multi(res2_correct$gr, tol = 0.5)
res2_update = ice_fase_mutlifurc(res2$tr,
                                 res2$sc_celltypes,
                                 total_time = res2$total_time,
                                 root_time = res2$root_time,
                                 gr = res2_correct$gr)

(plot_gr_clean(res2$gr,
              gr_node_time = c(res2$gr_trans_time, gr_tip_time),
              total_time = type_graph$target_time,
              type_col = gr_col_alt_mapped,
              show_node_label = T,
              node_label = gr_node_lab_mapped
              )) %>%
        push_pdf("example_gr_multifurc_mod.pdf", w = 2.5, h = 2.5, dir = plot_dir)

(plot_gr_clean(res2_update$gr,
              gr_node_time = c(res2_update$gr_trans_time, gr_tip_time),
              total_time = type_graph$target_time,
              type_col = gr_col_alt_mapped,
              show_node_label = T,
              node_label = gr_node_lab
) + scale_x_reverse()) %>%
        push_pdf("example_gr_multifurc_mod_merged.pdf", w = 2.5, h = 2.5, dir = plot_dir)

# use modified ICE_FASE to estimate multiple bias, assuming true topology
# res_list = future_map(rep_indices, function(i) {
#         tr3 = readRDS(paste0(rep_dir, "tr3/", stringr::str_pad(i, pad = "0", width = 4), ".rds"))
#         res = ice_fase_mutlifurc(tr3,
#                                  sc_celltypes = get_type_from_id(tr3$tip.label),
#                                  total_time = 11.5,
#                                  root_time = 0.6,
#                                  gr = tree_panel_mod$gr[[2]])
# }, .progress = T)
# time_all = map_dbl(res_list, function(x) x$gr_trans_time["11"])
# time_all_tb = tibble(time = time_all)
# (gghistogram(time_all_tb, x = "time", fill = "grey") +
#         geom_vline(xintercept = get_true_diff_time_mod3(tree_panel_mod$type_graph[[2]])["11"],
#                    color = "red",
#                    size = 2)) %>%
#         push_pdf("eg_multifurc_time", w = 2.5, h = 2.5, dir = plot_dir)
#
#
# size_all_tb = bind_rows(map(1:length(res_list), function(i) {
#         x = res_list[[i]]
#         size_vec = x$gr_node_sizes[["11"]]
#         tibble(exp_id = i,
#                node = names(size_vec),
#                size = size_vec)
# }))
# (size_all_tb %>% filter(node == "11") %>%
#         gghistogram(x = "size", fill = "grey") +
#         geom_vline(xintercept = get_true_size_mod2(type_graph_mod1)[["11"]][1],
#                    color = "red", size = 2)) %>%
#         push_pdf("eg_multifurc_size", w = 2.5, h = 2.5, dir = plot_dir)
#
# true_split = get_true_size_mod2(type_graph_mod1)[["11"]][2:4]
# true_split = true_split / sum(true_split)
# split_all_tb = size_all_tb %>% filter(node != "11") %>%
#         group_by(exp_id) %>%
#         mutate(frac = size / sum(size))
# (split_all_tb %>% ggbarplot(x = "exp_id", y = "frac", fill = "node", position =position_stack(), color = NA) +
#         scale_fill_manual(values = gr_col_vec) +
#         geom_hline(yintercept = c(true_split[3],
#                                   true_split[3] + true_split[2]),
#                    color = "red", size = 2)) %>%
#         push_pdf("eg_multifurc_split", w = 2.5, h = 2.5, dir = plot_dir)
#

# repeat above for birfucating example
res_list = future_map(1:nrow(exp_params), function(i) {
        tr3 = readRDS(paste0(rep_dir, "tr3/", stringr::str_pad(i, pad = "0", width = 4), ".rds"))
        graph_id = exp_params$big_graph_id[i]
        res = ice_fase_mutlifurc(tr3,
                                 sc_celltypes = get_type_from_id(tr3$tip.label),
                                 total_time = 11.5,
                                 root_time = 0.6,
                                 gr = tree_panel_mod$gr[[graph_id]])
}, .progress = T)
saveRDS(res_list, file = paste0(rep_dir, "res_list.rds"))

time_truth = get_true_diff_time_mod3(tree_panel_mod$type_graph[[1]])[["11"]]
res_list = readRDS(paste0(rep_dir, "res_list.rds"))
time_all = map_dbl(res_list, function(x) x$gr_trans_time["11"])
time_all_tb = tibble(time = time_all,
                     cond = factor(exp_params$big_graph_id))
time_all_tb %>% group_by(cond) %>%
        summarize(rmse = sqrt(mean(time - time_truth)^2))

(gghistogram(time_all_tb, x = "time", fill = "cond", add = "mean") %>%
                facet(facet.by = "cond", ncol = 1) +
                geom_vline(xintercept = get_true_diff_time_mod3(tree_panel_mod$type_graph[[2]])["11"],
                           color = "red",
                           size = 1)) %>%
        push_pdf("eg_multifurc_time", w = 2.5, h = 3.5, dir = plot_dir)

size_all_tb = bind_rows(map(1:length(res_list), function(i) {
        x = res_list[[i]]
        size_vec = x$gr_node_sizes[["11"]]
        tibble(exp_id = i,
               cond = factor(exp_params$big_graph_id[i]),
               node = names(size_vec),
               size = size_vec)
}))
size_truth = get_true_size_mod2(tree_panel_mod$type_graph[[1]])[["11"]][1]
size_all_tb %>% filter(node == "11") %>%
        group_by(cond) %>%
         summarize(rmse = sqrt(mean(log2(size/size_truth))^2))
(size_all_tb %>% filter(node == "11") %>%
                gghistogram(x = "size", fill = "cond", add = "mean") %>%
                facet(facet.by = "cond", ncol = 1) +
                geom_vline(xintercept = get_true_size_mod2(type_graph_mod1)[["11"]][1],
                           color = "red", size = 1)) %>%
        push_pdf("eg_multifurc_size", w = 2.5, h = 3.5, dir = plot_dir)

type_graph = tree_panel_mod$type_graph[[1]]
type_graph_mod1 = tree_panel_mod$type_graph[[2]]
true_split = get_true_size_mod2(type_graph)[["11"]][2:3]
true_split = true_split / sum(true_split)
true_split_mod = get_true_size_mod2(type_graph_mod1)[["11"]][2:4]
true_split_mod = true_split_mod / sum(true_split_mod)
split_all_tb = size_all_tb %>% filter(node != "11") %>%
        group_by(exp_id) %>%
        mutate(frac = size / sum(size))

split_all_tb %>% filter(node == "-1") %>%
        group_by(cond) %>%
        summarize(rmse = sqrt(mean(frac - true_split[1])^2))

split_all_ave = split_all_tb %>% group_by(cond, node) %>%
        summarise(mean_frac = mean(frac))

true_split_tibble = tibble(cond = c(1, 2, 2, 1, 2, 2),
                           mode = c(rep("truth", 3),
                                    rep("estimates", 3)),
                           val = c(true_split[2],
                                   true_split_mod[3],
                                   true_split_mod[3] + true_split_mod[2],
                                   split_all_ave$mean_frac[2],
                                   split_all_ave$mean_frac[5],
                                   split_all_ave$mean_frac[5] + split_all_ave$mean_frac[4]))
(split_all_tb %>% ggbarplot(x = "exp_id",
                            y = "frac", fill = "node",
                            position =position_stack(), color = NA) %>%
                facet(facet.by = "cond", ncol = 1) +
                scale_fill_manual(values = gr_col_vec) +
                geom_hline(aes(yintercept = val, linetype = mode),
                           data = true_split_tibble,
                           color = "red", size = 0.5)) %>%
        push_pdf("eg_multifurc_split", w = 2.5, h = 3.5, dir = plot_dir)

library(qfm)
library(lpSolve)
library(TreeSearch)
library(furrr)
plan(multisession, workers = 6)
source("R/graph_construction.R") # modified
source("analysis/panel/sample_MARC1_hgRNA.R")

# *****
output_dir = "./intermediate_data/panel_mod2_v1/"
plot_dir = "./plots/panel_mod2_v1/"
# dir.create(output_dir)
# dir.create(plot_dir)
# *****
make_type_graph_wrap <- function(num_tip, num_rep = 1000) {
        message(paste0('running num_tip: ', num_tip))
        bal_tbr = future_map(1:num_rep, function(i) {
                adjust_pec_v1(TBR(stree(n = num_tip, "balanced")))
        }, .progress = T, .options = furrr_options(seed = T))
        if (class(make_graph_from_phy_mod2(stree(num_tip, "left"))) == "type_graph") {
                message('can fit pec')
                pec_tbr = c(list(make_graph_from_phy_mod2(stree(num_tip, "left"))),
                            future_map(1:num_rep, function(i) {
                                    adjust_pec_v1(TBR(stree(n = num_tip, "left")))
                                    }, .progress = T, .options = furrr_options(seed = T)))
        } else {
                pec_tbr = future_map(1:num_rep, function(i) {
                        message(paste0('pec ', i))
                        adjust_pec_v1(stree(num_tip, "left"))
                }, .progress = T, .options = furrr_options(seed = T))
        }
        random_trees = future_map(rmtree(num_rep, n = num_tip), adjust_pec_v1, .progress = T)
        all_trees = tibble(num_tip = num_tip,
                           type_graph = c(list(adjust_pec_v1(stree(n = num_tip, "balanced"))),
                                          bal_tbr, random_trees, pec_tbr),
                           category = c("bal",
                                        rep("bal_tbr", length(bal_tbr)),
                                        rep("random", length(random_trees)),
                                        rep("pec_tbr", length(pec_tbr))))
        all_trees = all_trees[map_chr(all_trees$type_graph, class) == "type_graph", ]
        all_trees$bsum = future_map_dbl(all_trees$type_graph, function(x) {
                if (is.null(x)) {
                        return(NA)
                }
                calc_bsum(as.phylo_mod2.type_graph(x))
        }, .progress = T, .options = furrr_options(seed = T))
        all_trees$bsum_bin = cut(all_trees$bsum, breaks = c(-Inf, seq(0, max(all_trees$bsum, na.rm = T), by = 20), Inf))
        all_trees_sel = all_trees %>% nest(data = -c(num_tip, category, bsum_bin)) %>%
                mutate(data_sel = map(data, function(tb) {
                        tb[sample(nrow(tb), size = pmin(nrow(tb), 5), replace = F), ]
                })) %>% select(-data) %>%
                unnest(cols = data_sel)
        all_trees_sel
}
tree_panel = bind_rows(make_type_graph_wrap(16),
                       make_type_graph_wrap(32),
                       make_type_graph_wrap(64))
tree_panel %>% gghistogram(x = "bsum",
                           facet.by = "num_tip",
                           bins = 20,
                           scales = "free_x")
print(tree_panel, n = 500)
plot_type_graph_mod2(tree_panel$type_graph[[257]])
plot_type_graph_mod2(tree_panel$type_graph[[250]])
tree_panel$gr = map(tree_panel$type_graph, as.phylo_mod2.type_graph)
# removing extra pectinate for 16tip
which(tree_panel$num_tip == 16 & tree_panel$bsum == calc_bsum(stree(16, "left")))
tree_panel = tree_panel[-c(36:39), ]
tree_panel$category[which(tree_panel$num_tip == 16 & tree_panel$bsum == calc_bsum(stree(16, "left")))] = "pec"
saveRDS(tree_panel, file = paste0(output_dir, "tree_panel.rds"))

tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))

average_node_length <- function(type_graph) {
        edge_tb = type_graph$edges
        mean(dplyr::filter(edge_tb, out_node %in% type_graph$node_id)$length)
}
tree_panel$ave_node_len = map_dbl(tree_panel$type_graph, average_node_length)
(ggscatter(tree_panel, x = "bsum", y = "ave_node_len", facet.by = "num_tip", scales = "free_x") +
        geom_smooth(method = "lm", se = F) + ylab("Average time between commitment") +
        xlab("Colless index")) %>% push_pdf("node_len_bsum", w = 5.5, h = 2.7, dir = plot_dir)


# plotting distritubtions of tree panel
panel_node_params = bind_rows(map(1:nrow(tree_panel), function(i) {
        x = tree_panel$type_graph[[i]]
        diff_time = get_true_diff_time_mod3(x)
        start_time = get_true_start_time_mod3(x)
        diff_size = get_true_size_mod2(x)
        diff_bias = map_dbl(diff_size, function(x) x[2] / sum(x[2:3]))
        tibble(graph_id = i,
               num_tip = tree_panel$num_tip[i],
               bsum = tree_panel$bsum[i],
               node_rank = rank(diff_time[x$all_id]),
               start_time = start_time[x$all_id],
               end_time = diff_time[x$all_id],
               log2size = log2(map_dbl(diff_size[x$all_id], 1)),
               bias = diff_bias[x$all_id],
               double_time = unlist(x$lifetime[x$all_id]),
               prob_loss = unlist(x$prob_loss[x$all_id]),
               is_node = x$all_id %in% x$node_id)
}))
panel_node_params$node_rank = ceiling(panel_node_params$node_rank)

panel_node_params %>%
        filter(!is_node) %>%
        group_by(graph_id) %>%
        summarize(num_tip = num_tip[1],
                  total_terminal = sum(2^log2size)) %>%
        mutate(sampling_frac = num_tip * 100 / total_terminal) %>%
        summarise(mean(sampling_frac))



basic_theme = theme(text = element_text(size = 11),
              axis.text.x = element_text(angle = 30))
pdf(paste0(plot_dir, "panel_params.pdf"),
    width = 3.5, height = 2.25,
    onefile = T)
gghistogram(tree_panel, x = "bsum", facet.by = "num_tip", bins = 20,
            fill = "#949494",
            xlab = "Colless Index", ylab = "Count", scales = "free_x") +
        basic_theme
panel_node_params %>% gghistogram(x = "double_time",
                                  fill = "#949494",
                                  xlab = "Doubling time") + basic_theme
panel_node_params %>% gghistogram(x = "prob_loss",
                                  fill = "#949494",
                                  xlab = "Cell loss rate") + basic_theme
panel_node_params %>% gghistogram(x = "bias",
                                  fill = "#949494",
                                  xlab = "Commitment bias") + basic_theme
panel_node_params %>% filter(is_node) %>%
        gghistogram(x = "log2size",
                    fill = "#949494",
                    xlab = "log2(Field size)") + basic_theme
panel_node_params %>%
        ggscatter(x = "end_time", y = "double_time",
                  xlab = "Commitment time",
                  ylab = "Doubling time",
                  alpha = 0.025, shape = 16) +
        geom_smooth(se = F, method = "loess")
dev.off()

(panel_node_params %>%
        ggscatter(x = "end_time", y = "double_time",
                  xlab = "Commitment time",
                  ylab = "Doubling time",
                  alpha = 0.025, shape = 16) +
        geom_smooth(se = F, method = "loess")) %>%
        push_png("doubling_time_vs_time", width = 3.5, height = 2.25, dir = plot_dir,
                 res = 500)

(panel_node_params %>%
                filter(is_node) %>%
        ggboxplot(x = "node_rank", y = "end_time", facet.by = "num_tip",
                  size = 0.2,
                  outlier.size = 0.1,
                  xlab = "Rank of commitment event time",
                  ylab = "Time",
                  scales = "free_x", numeric.x.axis = T) +
        theme(text = element_text(size = 11),
              axis.text.x = element_text(size = 8, angle = 90))) %>%
        push_pdf(file_name = "node_rank_time", w = 7.75, h = 2.8, dir = plot_dir)

tr_cat_col = RColorBrewer::brewer.pal(5, "Set2")
names(tr_cat_col) = c("bal", "bal_tbr", "random", "pec_tbr", "pec")
pdf(paste0(plot_dir, "graph_pco.pdf"), w = 6, h = 3.5, onefile = T)
for (num_tip in c(16, 32, 64)) {
        gr_list = tree_panel$gr[tree_panel$num_tip == num_tip]
        class(gr_list) = "multiPhylo"
        tree_sp = treespace::treespace(gr_list, nf = 2)
        pco_mat = as_tibble(tree_sp$pco$tab)
        pco_mat$bsum = tree_panel$bsum[tree_panel$num_tip == num_tip]
        pco_mat$category = tree_panel$category[tree_panel$num_tip == num_tip]
        pco_mat$category = factor(pco_mat$category, levels = rev(names(tr_cat_col)))
        pco_mat = arrange(pco_mat, category)
        print(ggscatter(pco_mat, x = "A1", y = "A2", col = "category", size = 0.5,
                        title = paste0("Number of states: ", num_tip),
                        xlab = "Dimension 1", ylab = "Dimension 2") + scale_color_manual(values = tr_cat_col) +
                      (ggscatter(pco_mat, x = "A1", y = "A2", col = "bsum", size = 0.5,
                                 xlab = "Dimension 1", ylab = "Dimension 2") + scale_color_distiller(palette = "Spectral")) +
                      plot_layout(guide = "collect") & theme(legend.position = "bottom", text = element_text(size = 12)))
}
dev.off()

# generate number of cells over time
generate_time_size <- function(type_graph) {
        gens0 = make_gens_mod2(type_graph)
        count_types = map(gens0, function(x) {
                tibble(start = x$start_time + 0:x$num_gen *x$double_time,
                       end = x$start_time + 1:(x$num_gen+1) *x$double_time,
                       count = rowSums(x$cell_mode_counts[, c(1, 3), drop = F]))
        })
        count_at_time <- function(tt, count_types) {
                map_dbl(count_types,function(x) {
                        out = x$count[x$start <= tt & x$end > tt]
                        if (length(out) == 0) {
                                return(0)
                        } else {
                                return(out)
                        }
                })
        }
        time_size_tb = bind_rows(map(seq(0, type_graph$target_time, by = 0.05), function(tt) {
                tibble(time = tt, count = sum(count_at_time(tt, count_types)))
        })) %>% mutate(log2count = log2(count))
        time_size_tb
}
panel_time_size = bind_rows(map(1:nrow(tree_panel), function(i) {
        x = tree_panel$type_graph[[i]]
        out = generate_time_size(x)
        out$graph_id = i
        out
}))

obs_tb = tibble(time = c(0, 2., 2.5, 3.5, 4.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0),
                log2count = log2(c(c(1, 4, 12, 32, 64), c(120, 250, 660, 4500, 15000, 75000) * 2)))

(panel_time_size %>% ggplot() + geom_line(aes(x = time, y = log2count, group = graph_id), size = 0.25, alpha = 0.15, color = "#bdbdbd") +
        geom_line(data = group_by(panel_time_size, time) %>% summarize(log2_mean_count = log2(mean(count))),
                  aes(x = time, y = log2_mean_count), color = "#377eb8", size = 1.) +
        geom_point(data = obs_tb, aes(x = time, y = log2count), col = "#e41a1c") +
                xlab("Time") + ylab("log2(Cell count)") + theme_pubr()) %>%
        push_pdf("time_size", dir = plot_dir, w= 3.5, h = 2.25)

plot_type_graph_mod2(tree_panel$type_graph[[23]])

map_dbl(tree_panel$type_graph, function(x) {
        if (is.null(x)) {
                return(0)
        }
        length(make_gens_mod2(x))
})














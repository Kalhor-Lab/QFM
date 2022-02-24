# produces phylogenies at smaller number of barcoding elements
output_dir = "./intermediate_data/panel/"
exp_params = readRDS(paste0(output_dir, "exp_data_10rep_wtree.rds"))
all_graphs = readRDS(paste0(output_dir, "all_graphs.rds"))

library(furrr)
plan(multisession, workers = 12)
total_time = 15.0 - 0.6
exp_indices = which(exp_params$num_tip == 16)
tr3_collect = map(exp_indices, function(j) {
        message(j)
        map(c(25, 50), function(num_elements) {
                tr3 = phylotime(exp_params$sc[[j]][, 1:num_elements],
                                t_total = total_time,
                                mut_p = subset_mut_p(exp_params$mut_p[[j]],
                                                     1:num_elements))
        })
})
saveRDS(tr3_collect, "./intermediate_data/panel/phylotime_25_50_hgRNA.rds")

tr5_collect = map(exp_indices, function(j) {
        message(j)
        map(c(25, 50), function(num_elements) {
                x = exp_params$data[[j]]
                if (is.null(x$sc)) {
                        return(NULL)
                }
                sc_onehot = barcode_allele_onehot_new(x$sc[, 1:num_elements])
                sc_onehot = sc_onehot[, colSums(sc_onehot) > 1]
                sc_onehot = sc_onehot[sample(nrow(sc_onehot), replace = F), ]

                dmat = proxyC::dist(sc_onehot, method = "manhattan")
                tr = phangorn::upgma(dmat)

                total_depth = max(node.depth.edgelength(tr))
                tr$edge.length = tr$edge.length / total_depth * total_time
                tr = name_nodes(tr)
                tr
        })
})
saveRDS(tr5_collect, "./intermediate_data/panel/hamming_25_50_hgRNA.rds")

exp_indices = which(exp_params$num_tip == 16)
tr3_collect = readRDS("./intermediate_data/panel/phylotime_25_50_hgRNA.rds")
tr5_collect = readRDS("./intermediate_data/panel/hamming_25_50_hgRNA.rds")

dmat_tb = bind_rows(map(c(25, 50, 100), function(num_elements) {
        dmat = cophenetic(exp_params$tr[[1]])
        if (num_elements == 100) {
                dmat3 = cophenetic(exp_params$tr3[[1]])
        }
        if (num_elements == 50) {
                dmat3 = cophenetic(tr3_collect[[1]][[2]])
        }
        if (num_elements == 25) {
                dmat3 = cophenetic(tr3_collect[[1]][[1]])
        }
        dmat3 = dmat3[rownames(dmat), colnames(dmat)]
        tibble(num_elements = num_elements,
               dist = c(dmat),
               dist3 = c(dmat3))
}))

# Figure 8D
# define dmat and dmat3, then
dmat3 = as.matrix(dmat3)
class(dmat3) = dmat3
Heatmap(dmat3)
dmat_plot_tb = filter(dmat_tb, dist > 0)
dmat_plot_tb = dmat_plot_tb[sample(which(dmat_plot_tb$num_elements == 100), size = 1e5), ]
(ggscatter(dmat_plot_tb,
           x = "dist",
           y = "dist3",
           size = 0.25) + # ,
          #xlim = c(15, 30),
          #ylim = c(15, 30)) +
        xlab("True Cophenetic Dist") +
        ylab("Reconstructed Cophenetic Dist") +
        stat_cor() +
        geom_abline(size = 0.5) + theme(text = element_text(size = 10))) %>%
        push_png("tr_cophenetic_scatter", w = 3, h = 2.5, ps = 10, res = 1200, dir = "plots/panel/phylotime")

# dmat_tb = mutate(dmat_tb, error = dist3 - dist)
# ggboxplot(dmat_tb[sample(nrow(dmat_tb), size = 1e5), ],
#           x = "num_elements", y = "error")

kc_tb3 = bind_rows(map(1:length(exp_indices), function(i) {
        j = exp_indices[i]
        bind_rows(
                tibble(j = j,
                       num_element = 100,
                       kc0 = treespace::treeDist(exp_params$tr3[[j]], exp_params$tr[[j]], lambda = 0),
                       kc1 = treespace::treeDist(exp_params$tr3[[j]], exp_params$tr[[j]], lambda = 1),
                       rf = phangorn::RF.dist(exp_params$tr3[[j]], exp_params$tr[[j]])),
                tibble(j = j,
                       num_element = 50,
                       kc0 = treespace::treeDist(tr3_collect[[i]][[2]], exp_params$tr[[j]], lambda = 0),
                       kc1 = treespace::treeDist(tr3_collect[[i]][[2]], exp_params$tr[[j]], lambda = 1),
                       rf = phangorn::RF.dist(tr3_collect[[i]][[2]], exp_params$tr[[j]])),
                tibble(j = j,
                       num_element = 25,
                       kc0 = treespace::treeDist(tr3_collect[[i]][[1]], exp_params$tr[[j]], lambda = 0),
                       kc1 = treespace::treeDist(tr3_collect[[i]][[1]], exp_params$tr[[j]], lambda = 1),
                       rf = phangorn::RF.dist(tr3_collect[[i]][[1]], exp_params$tr[[j]]))
        )
}))
kc_tb5 = bind_rows(map(1:length(exp_indices), function(i) {
        j = exp_indices[i]
        bind_rows(
                tibble(j = j,
                       num_element = 100,
                       kc0 = treespace::treeDist(exp_params$tr5[[j]], exp_params$tr[[j]], lambda = 0),
                       kc1 = treespace::treeDist(exp_params$tr5[[j]], exp_params$tr[[j]], lambda = 1),
                       rf = phangorn::RF.dist(exp_params$tr5[[j]], exp_params$tr[[j]])),
                tibble(j = j,
                       num_element = 50,
                       kc0 = treespace::treeDist(tr5_collect[[i]][[2]], exp_params$tr[[j]], lambda = 0),
                       kc1 = treespace::treeDist(tr5_collect[[i]][[2]], exp_params$tr[[j]], lambda = 1),
                       rf = phangorn::RF.dist(tr5_collect[[i]][[2]], exp_params$tr[[j]])),
                tibble(j = j,
                       num_element = 25,
                       kc0 = treespace::treeDist(tr5_collect[[i]][[1]], exp_params$tr[[j]], lambda = 0),
                       kc1 = treespace::treeDist(tr5_collect[[i]][[1]], exp_params$tr[[j]], lambda = 1),
                       rf = phangorn::RF.dist(tr5_collect[[i]][[1]], exp_params$tr[[j]]))
        )
}))
exp_cassi = readr::read_csv("./intermediate_data/panel/phylo_cassiopeia_summary.csv")
exp_cassi = exp_cassi %>% dplyr::transmute(j = row_num, num_element = hgrnaNum, rf = rf, kc0 = kc0, kc1 = kc1)
kc_tb = bind_rows(mutate(kc_tb3, method = "Phylidite"),
                  mutate(kc_tb5, method = "Hamming"),
                  mutate(exp_cassi, method = "Cassiopeia"))
# summary table to look at the statistics here.
kc_tb$j = factor(kc_tb$j)
(ggline(kc_tb, x = "num_element", y = "rf", color = "method", add = "mean_se", ylim = c(0, NA), xlab = "Number of Sites") +
                ggline(kc_tb, x = "num_element", y = "kc0", color = "method", add = "mean_se", ylim = c(0, NA), xlab = "Number of Sites") +
                ggline(kc_tb, x = "num_element", y = "kc1", color = "method", add = "mean_se", ylim = c(0, NA), xlab = "Number of Sites") +
                plot_layout(guide = "collect")) %>%
       push_pdf(file_name = "phylo_eval_kc", w = 8, h = 2.5, ps = 10, dir = "./plots/panel/")

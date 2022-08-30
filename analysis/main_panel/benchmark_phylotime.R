library(qfm)
source("./analysis/MARC1/load_MARC1_data.R")

library(furrr)
plan(multisession, workers = 12)

# collecting and plotting Phylotime, Cassiopeia and Hamming comparisons
exp_params = readRDS("./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all.rds")
tree_panel = readRDS("./intermediate_data/panel_mod2_v1/tree_panel.rds")
exp_params$num_tip = tree_panel$num_tip[exp_params$big_graph_id]
exp_id_all = which(exp_params$num_tip == 16)

run_cass_tr <- function(cass_dir, data_dir, exp_id, tr_dir) {
        x = readRDS(paste0(data_dir, stringr::str_pad(exp_id, pad = "0", width = 4), ".rds"))
        recon_file = paste(cass_dir,
                           paste0("exp_", stringr::str_pad(exp_id, pad = "0", width = 4), "0001_recon.txt"), sep = "/")
        if (file.exists(recon_file)) {
                recon <- ape::read.tree(recon_file)
                reconT <- compute.brlen(recon, method = "Grafen")
                reconTest = chronos(reconT, lambda=0)
                reconTest <- ape::multi2di(reconTest, random = T)

                # recon = collapse.singles(recon)
                total_time <- max(node.depth.edgelength(x$tr))
                total_len <- max(node.depth.edgelength(reconTest))

                # scaledCassiopeia1<-rescale.rooted.tree(consensus_ultra, scale.height = old/new)
                scaledreconTest = rlang::duplicate(reconTest)
                scaledreconTest$edge.length = scaledreconTest$edge.length / total_len * total_time
                out_file = paste0(tr_dir, stringr::str_pad(exp_id, width = 4, pad = "0"), ".rds")
                saveRDS(scaledreconTest, file = out_file)
        }
}

cass_dir = "./intermediate_data/panel_mod2_v1/cassi_trees/qfm25_recon/"
tr_dir = "./intermediate_data/panel_mod2_v1/cassi_trees/processed_tr25/"

cass_dir = "./intermediate_data/panel_mod2_v1/cassi_trees/qfm50_recon/"
tr_dir = "./intermediate_data/panel_mod2_v1/cassi_trees/processed_tr50/"

cass_dir = "./intermediate_data/panel_mod2_v1/cassi_trees/qfm100_recon/"
tr_dir = "./intermediate_data/panel_mod2_v1/cassi_trees/processed_tr100/"

# dir.create(tr_dir)
future_walk(exp_id_all, function(j) {
        run_cass_tr(cass_dir,
                    "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/",
                    j,
                    tr_dir)
}, .progress = T, .options = furrr_options(seed = T))

run_hamming_tr <- function(data_dir, j, tr_dir, mut_indices) {
        message(j)
        x = readRDS(paste0(data_dir, stringr::str_pad(j, pad = "0", width = 4), ".rds"))
        sc_onehot = barcode_allele_onehot_new(x$sc[, mut_indices])
        sc_onehot = sc_onehot[, colSums(sc_onehot) > 1]
        sc_onehot = sc_onehot[sample(nrow(sc_onehot), replace = F), ]

        dmat = proxyC::dist(sc_onehot, method = "manhattan")
        tr = phangorn::upgma(dmat)

        total_time = max(node.depth.edgelength(x$tr))
        total_depth = max(node.depth.edgelength(tr))
        tr$edge.length = tr$edge.length / total_depth * total_time
        tr = name_nodes(tr)
        out_file = paste0(tr_dir, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
        saveRDS(tr, file = out_file)
}
mut_indices = readRDS("./intermediate_data/panel_mod2_v1/mut_indices25.rds")
tr_dir = "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/hamming25/"

mut_indices = readRDS("./intermediate_data/panel_mod2_v1/mut_indices50.rds")
tr_dir = "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/hamming50/"

mut_indices = readRDS("./intermediate_data/panel_mod2_v1/mut_indices100.rds")
tr_dir = "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/hamming100/"

dir.create(tr_dir)
future_walk(exp_id_all, function(j) {
        run_hamming_tr("./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/",
                       j,
                       tr_dir,
                       mut_indices)
}, .progress = T, .options = furrr_options(seed = T))

tb_hamming = tibble(num_element = c(25, 50, 100),
       method = "hamming",
       dir = c("./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/hamming25/",
               "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/hamming50/",
               "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/hamming100/")) %>%
        mutate(tr_tb = map(dir, function(x) {
                tibble(exp_id = exp_id_all,
                       tr = map(exp_id_all, function(j) {
                               tr_file = paste0(x, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                               if (file.exists(tr_file)) {
                                       return(readRDS(tr_file))
                               } else {
                                       return(NULL)
                               }
                               }))
        }))
tb_hamming = unnest(tb_hamming, tr_tb)
tb_hamming$tr_eval = map2(tb_hamming$exp_id, tb_hamming$tr, function(exp_id, tr_r) {
        x = readRDS(paste0(data_dir, stringr::str_pad(exp_id, pad = "0", width = 4), ".rds"))
        if (!is.null(tr_r)) {
                tr = x$tr
                evalute_gr(tr, tr_r)
        } else {
                return(NULL)
        }
})
saveRDS(tb_hamming, file = "./intermediate_data/panel_mod2_v1/tb_hamming.rds")

data_dir = "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/"
tb_cass = tibble(num_element = c(25, 50, 100),
                 method = "cass",
                 dir = c("./intermediate_data/panel_mod2_v1/cassi_trees/processed_tr25/",
                         "./intermediate_data/panel_mod2_v1/cassi_trees/processed_tr50/",
                         "./intermediate_data/panel_mod2_v1/cassi_trees/processed_tr100/")) %>%
        mutate(tr_tb = map(dir, function(x) {
                tibble(exp_id = exp_id_all,
                       tr = map(exp_id_all, function(j) {
                               tr_file = paste0(x, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                               if (file.exists(tr_file)) {
                                       return(readRDS(tr_file))
                               } else {
                                       return(NULL)
                               }
                       }))
        }))
tb_cass = unnest(tb_cass, tr_tb)
tb_cass$tr_eval = map2(tb_cass$exp_id, tb_cass$tr, function(exp_id, tr_r) {
        x = readRDS(paste0(data_dir, stringr::str_pad(exp_id, pad = "0", width = 4), ".rds"))
        if (!is.null(tr_r)) {
                tr = x$tr
                evalute_gr(tr, tr_r)
        } else {
                return(NULL)
        }
})
saveRDS(tb_cass, file = "./intermediate_data/panel_mod2_v1/temp_tb_cass.rds")

tb_phylotime = tibble(num_element = c(25, 50, 100),
                      method = "phylotime",
                      dir = c("./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/tr3_25site/",
                              "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/tr3_50site/",
                              "./intermediate_data/panel_mod2_v1/exp_data_5rep_mut_all/tr3_100site/")) %>%
        mutate(tr_tb = map(dir, function(x) {
                tibble(exp_id = exp_id_all,
                       tr = map(exp_id_all, function(j) {
                               tr_file = paste0(x, stringr::str_pad(j, width = 4, pad = "0"), ".rds")
                               if (file.exists(tr_file)) {
                                       return(readRDS(tr_file))
                               } else {
                                       return(NULL)
                               }
                       }))
        }))
tb_phylotime = unnest(tb_phylotime, tr_tb)
tb_phylotime$tr_eval = map2(tb_phylotime$exp_id, tb_phylotime$tr, function(exp_id, tr_r) {
        x = readRDS(paste0(data_dir, stringr::str_pad(exp_id, pad = "0", width = 4), ".rds"))
        if (!is.null(tr_r)) {
                tr = x$tr
                evalute_gr(tr, tr_r)
        } else {
                return(NULL)
        }
})
saveRDS(tb_phylotime, file = "./intermediate_data/panel_mod2_v1/tb_phylotime.rds")

tb_cass = readRDS("./intermediate_data/panel_mod2_v1/temp_tb_cass.rds")
tb_hamming = readRDS("./intermediate_data/panel_mod2_v1/tb_hamming.rds")
tb_phylotime = readRDS("./intermediate_data/panel_mod2_v1/tb_phylotime.rds")

tb_all = bind_rows(tb_hamming, tb_cass, tb_phylotime)
tb_all$kc0 = map_dbl(tb_all$tr_eval, function(x) {
        if (is.null(x)) {
                return(NA)
        } else {
                return(x$kc0)
        }
})
tb_all$kc1 = map_dbl(tb_all$tr_eval, function(x) {
        if (is.null(x)) {
                return(NA)
        } else {
                return(x$kc1)
        }
})

tb_all %>% group_by(method) %>%
        summarize(mean(kc1))


g1 = tb_all %>%
        filter(!is.na(kc0)) %>%
        ggline(x = "num_element", y = "kc0", color = "method", add = "mean_se", ylim = c(0, NA),
               xlab = "Number of sites")
g2 = tb_all %>%
        filter(!is.na(kc0)) %>%
        ggline(x = "num_element", y = "kc1", color = "method", add = "mean_se", ylim = c(0, NA),
               xlab = "Number of sites")
(g1 + g2 + plot_layout(guides = "collect") &
                theme(legend.position = "bottom")) %>%
        push_pdf("phylotime_eval",
                 w = 5.2, h = 3.6,
                 dir = "./plots/panel_mod2_v1/")



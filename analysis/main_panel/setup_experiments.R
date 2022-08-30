# setting up experiment and hgRNA sets
library(qfm)
source("./analysis/MARC1/load_MARC1_data.R")

# large set of barcoding sites
set.seed(73)
mut_p_all = concat_mut_p_list(list(sample_mutp(500, c(0.6, 11.5), "fast"),
                                   sample_mutp(500, c(0.6, 11.5), "mid"),
                                   sample_mutp(500, c(0.6, 11.5), "slow"),
                                   make_som_mut_p(2000, c(0.6, 11.5))))
saveRDS(mut_p_all, file = paste0(output_dir, "mut_p_all.rds"))

# a smaller set of barcoding sites mimicing 50hgRNAs from MARC1
set.seed(73)
mut_p_all = concat_mut_p_list(list(sample_mutp(14, c(0.6, 11.5), "fast"),
                                   sample_mutp(36, c(0.6, 11.5), "mid")))
saveRDS(mut_p_all, file = paste0(output_dir, "mut_p_marc1.rds"))

# set of experiments with 20 replicates of 16-type fate maps
exp_name = "5rep_mut_all"
exp_params = expand.grid(big_graph_id = 1:nrow(tree_panel),
                         sample_size = c(100),
                         sampling = c("fixed", "proportional"),
                         i_sim = 1:5) %>% as_tibble()
saveRDS(exp_params, file = paste0(output_dir, "exp_data_", exp_name, ".rds"))

# set of experiments with 20 replicates of 16-type fate maps
exp_name = "20rep_16type"
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))
exp_params = expand.grid(big_graph_id = which(tree_panel$num_tip == 16),
                         sample_size = c(50, 100, 200),
                         sampling = c("fixed", "proportional"),
                         i_sim = 1:10) %>% as_tibble()
saveRDS(exp_params, file = paste0(output_dir, "exp_data_", exp_name, ".rds"))


# indices of a different set of barcoding sites from "mut_p_all.rds" above
#TODO: slow and fast need to be swapped, manuualy fixed now
all_mut_indices = bind_rows(as_tibble(expand.grid(num_element = c(25, 50, 100),
                                                  class = c("slow", "mid", "fast"))),
                            tibble(num_element = c(250, 500, 1000),
                                   class = "somatic"))
all_mut_indices$mut_indices = c(purrr::reduce(map(c(0, 500, 1000), function(i_start) {
        map(list(1:25,
                 26:75,
                 76:175), function(i_seq) {
                         i_start + i_seq
                 })
}), c),
list(1501:1750,
     1751:2250,
     2251:3250))
saveRDS(all_mut_indices, paste0(output_dir, "all_mut_indices.rds"))

# below are sets of mixed hgRNAs similar to MARC1
mut_indices25 = c(190:196, 712:729)
saveRDS(mut_indices25, file = paste0(output_dir, "mut_indices25.rds"))
mut_indices50 = c(176:189, 676:711)
saveRDS(mut_indices50, file = paste0(output_dir, "mut_indices50.rds"))
mut_indices100 = c(197:224, 730:801)
saveRDS(mut_indices100, file = paste0(output_dir, "mut_indices100.rds"))


output_dir = "./intermediate_data/panel_mod2_v1/"
tree_panel = readRDS(paste0(output_dir, "tree_panel.rds"))
mut_p_all = readRDS(paste0(output_dir, "mut_p_all.rds"))

# testing phylotime on one set of indices
# data_out = readRDS(paste0(output_dir, "exp_data_5rep_mut_all/0001.rds"))
# all_mut_indices$sc = map(all_tr$mut_indices, function(m_indices) {
#         data_out$sc[, m_indices]
# })
# all_mut_indices$mut_p = map(all_tr$mut_indices, function(m_indices) {
#         subset_mut_p(mut_p_all, m_indices)
# })
# all_mut_indices$target_time = 11.5
#
# m_indices = all_tr$mut_indices[[1]]
# tr = phylotime(x$sc[, m_indices],
#                t_total = 11.5 - 0.6,
#                mut_p = subset_mut_p(mut_p_all, m_indices))
# end testing code

# computing average mutation rates
mut_p_all = readRDS("./intermediate_data/panel_mod2_v1/mut_p_all.rds")
mean(mut_p_all$mut_rate[1:500])
mean(mut_p_all$mut_rate[500:1000])
mean(mut_p_all$mut_rate[1000:1500])



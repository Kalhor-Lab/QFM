# construct_type_pair_df <- function(all_celltype) {
#         type_pair_mat = t(combn(all_celltype, 2))
#         colnames(type_pair_mat) = c("a", "b")
#         type_pairs_df = as_tibble(type_pair_mat)
#         type_pairs_df
# }
# find_restriction_nodes <- function(tr, tip_celltype, type_pairs_df) {
#         # this functions add a column to type_paris_df
#         out_list = list_dd_and_tips(tr)
#         tr_dd = out_list$dd
#         tr_tip_list = out_list$tips
#         # tr_dd = map(adephylo::listDD(tr), names)
#         # tr_tip_list = map(adephylo::listTips(tr), names)
#
#         tr_tip_type_list = map(tr_tip_list, function(x) unique(tip_celltype[x]))
#         # tr_tip_type_counts = map(tr_tip_list, function(x) table(tip_celltype[x]))
#         tip_ind = rep(F, length(tr$tip.label))
#         names(tip_ind) = tr$tip.label
#         # building list of restrition events
#         message('processing pairs...')
#         type_pairs_df$restriction_nodes = future_map(1:nrow(type_pairs_df), function(i) {
#                 t1 = type_pairs_df$a[i]
#                 t2 = type_pairs_df$b[i]
#                 node_ind = map_lgl(tr_tip_type_list, function(x) all(c(t1, t2) %in% x))
#                 node_ind = c(node_ind, tip_ind)
#                 out = names(tr_dd)[which(map_lgl(names(tr_dd), function(x) {
#                         d1 = tr_dd[[x]][1]
#                         d2 = tr_dd[[x]][2]
#                         node_ind[x] & !node_ind[d1] & !node_ind[d2]
#                 }))]
#                 out
#         })
#         type_pairs_df
# }

# estimation_pairwise_restrictions <- function(tr1, sc_celltypes = NULL) {
#         tr1 = name_nodes(tr1)
#         tr1$edge.length[tr1$edge.length < 0] = 0
#         if (is.null(sc_celltypes)) {
#                 sc_celltypes = get_type_from_id(tr1$tip.label)
#                 names(sc_celltypes) = tr1$tip.label
#         }
#         type_pairs_df = construct_type_pair_df(unique(sc_celltypes))
#         type_pairs_df = find_restriction_nodes(tr1, sc_celltypes, type_pairs_df)
#         type_pairs_df = type_pairs_df %>% unnest(restriction_nodes) %>% mutate(restriction_time = map_dbl(restriction_nodes, function(r) {
#                 get_node_time(tr1, r)
#         })) %>% nest(restriction_nodes = -c(a, b))
#         type_pairs_df
# }
# restrction_to_graph <- function(type_pairs_df, total_time, stat_func = mean) {
#         dist_df = type_pairs_df %>%
#                 mutate(dist = map_dbl(restriction_nodes, function(x) {
#                         (total_time - stat_func(x$restriction_time))*2
#                 })) %>%
#                 mutate(N1 = a, N2 = b, Var1 = a, Var2 = b)
#         dmat = dist_df2mat(dist_df)
#         gr = phangorn::upgma(dmat)
#         gr
# }
# sample_gaussian <- function(mat, q, cycle) {
#         matrix(map_dbl(c(mat), function(x) {
#                 if (x == 0) {
#                         return(0)
#                 }
#                 mu = x * (1+q)^cycle
#                 phi = (1-q)/(x *(1+q))
#                 sigma_sq = phi * mu^2
#                 rnorm(n = 1, mean = mu, sd = sqrt(sigma_sq))
#         }), nrow = nrow(mat))
# }

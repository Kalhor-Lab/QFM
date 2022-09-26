library(qfm)
library(Biostrings)
# conflicted::conflict_prefer("purrr::reduce", "purrr")
# conflicted::conflict_prefer("filter", "dplyr")
set.seed(121)
get_type <- function(x) {
        out = str_match(x, pattern = "(G\\d)T(\\d)D(.*)")[, 3]
        names(out) = x
        paste0("T", out)
}
get_group <- function(x) {
        out = str_match(x, pattern = "(G\\d)T(\\d)D(.*)")[, 2]
        names(out) = x
        out
}
type_col = c("3" = "#f768a1",
             "-5" = "#dd3497",
             "-6" = "#7a0177",
             "5" = "#c7e9b4",
             "1" = "#80cdc1",
             "-1" = "#35978f",
             "-2" = "#01665e",
             "2" = "#dfc27d",
             "-3" = "#d8b365",
             "-4" = "#8c510a",
             "4" = "#1d91c0")

plot_dir = "./plots/iPSC/"
count_dir_g1 = "./data/iPSC/G1-partial/A-trimmed_raw_data/"
count_dir_g2 = "./data/iPSC/G2-partial/A-trimmed_raw_data/"

# raw count files
pair_count_files = c(list.files(count_dir_g1,
                                pattern = ".*_barcode_gRNA_pairing_counts\\.txt",
                                full = T),
                     list.files(count_dir_g2,
                                pattern = ".*_barcode_gRNA_pairing_counts\\.txt",
                                full = T))
pair_file_match = str_match(pair_count_files, "(G\\d)T(\\d+)[A|D](\\d+)-(.*?)_barcode_gRNA_pairing_counts\\.txt")

all_samples = tibble(file = pair_count_files,
                     batch = pair_file_match[, 5],
                     group = pair_file_match[, 2],
                     cell_type = pair_file_match[, 3],
                     cell_num = pair_file_match[, 4],
                     cell_index = paste0(group, "T", cell_type, "D", cell_num))
all_samples = filter(all_samples, cell_index != "NATNADNA")

#### Merging sequencing data from different runs ####
# all_samples$file_barcode = stringr::str_replace(all_samples$file,
#                                                 "_barcode_gRNA_pairing_counts",
#                                                 "_barcode_counts")
# all_samples$id_spacer_data = map(all_samples$file, read_id_spacer_counts)
# all_samples$id_data = map(all_samples$file_barcode, read_id_counts)

# all_samples_nest = all_samples %>% nest(nested = -cell_index)
# all_samples_nest$id_spacer_data = map(all_samples_nest$nested, function(x) {
#         assertthat::assert_that(nrow(x) <= 3)
#         if (nrow(x) == 3) {
#                 out = merge_id_spacer_counts(x$id_spacer_data[[1]], x$id_spacer_data[[2]])
#                 out = merge_id_spacer_counts(out, x$id_spacer_data[[3]])
#         }
#         if (nrow(x) == 2) {
#                 out = merge_id_spacer_counts(x$id_spacer_data[[1]], x$id_spacer_data[[2]])
#         } else {
#                 out = x$id_spacer_data[[1]]
#         }
#         return(out)
# })
# all_samples_nest$id_data = map(all_samples_nest$nested, function(x) {
#         assertthat::assert_that(nrow(x) <= 3)
#         if (nrow(x) == 3) {
#                 out = merge_id_counts(x$id_data[[1]], x$id_data[[2]])
#                 out = merge_id_counts(out, x$id_data[[3]])
#         }
#         if (nrow(x) == 2) {
#                 out = merge_id_counts(x$id_data[[1]], x$id_data[[2]])
#         } else {
#                 out = x$id_data[[1]]
#         }
#         return(out)
# })
# output_dir = ""./data/iPSC/combined/1-pair_counting/"
# dir.create(output_dir, recursive = T)
# for (i in 1:nrow(all_samples_nest)) {
#         write_tsv(x = all_samples_nest$id_spacer_data[[i]],
#                   file = paste0(output_dir, all_samples_nest$cell_index[i], "_barcode_gRNA_pairing_counts.txt"),
#                   col_names = F)
# }
# for (i in 1:nrow(all_samples_nest)) {
#         write_tsv(x = all_samples_nest$id_data[[i]],
#                   file = paste0(output_dir, all_samples_nest$cell_index[i], "_barcode_counts.txt"),
#                   col_names = F)
# }
##### end merging ####

# read the processed files from combined data
all_samples_nest = all_samples %>% nest(nested = -cell_index)
all_samples_nest$filtered_file = paste0("./data/iPSC/combined/4-pair_filtering/",
                                        all_samples_nest$cell_index,
                                        "_filteredpairs.txt")

# loading additional test plate data
data_files_test_plate = list.files(path = "./data/iPSC/G1_testplate/4-pair_filtering/",
                                   pat = ".*_filteredpairs\\.txt", full = T)
data_match_test_plate = str_match(basename(data_files_test_plate), "(c\\d+)-T1A-KL-CB11-SCtest_filteredpairs\\.txt")
all_samples_test_plate = tibble(filtered_file = data_files_test_plate,
                                cell_index = paste0("G1T1D", data_match_test_plate[, 2]))
all_samples_nest = bind_rows(all_samples_nest, all_samples_test_plate)

# adding the additional test plate manually
all_samples_nest$processed_data = map(all_samples_nest$filtered_file,
                                      function(x) {
                                              out = read_id_spacer_counts(x)
                                              if (is.null(out)) {
                                                      return(NULL)
                                              } else {
                                                      out = out %>% nest(spacer = c(spacer, count))
                                                      return(out)
                                              }
                                      })

all_samples_unest = select(all_samples_nest, cell_index, processed_data) %>% unnest(processed_data)

# loading parent files
parent_file = "./data/iPSC/EC96Lin12-PBGeno_filteredpairs.txt"
parent_data = read_id_spacer_counts(parent_file) %>%
        nest(spacer = c(spacer, count)) %>%
        dplyr::rename(parent_spacer = spacer)


spacer_unmutated_GCCAAAAGCT = "GAAACACCGGTGGTCGCCGTGGAGAGTGGTGGGGTTAGAGCTAGAAATAG"
all_sample_id = all_samples_unest %>% group_by(id) %>% summarize(total = n()) %>% arrange(desc(total))
all_sample_id$parent = all_sample_id$id %in% c(parent_data$id, "GCCAAAAGCT")
print(all_sample_id, n = 35) # checking all the IDs
# parent_miss_id = parent_data$id[!parent_data$id %in% all_sample_id$id[1:32]]

# merging error id sequences
parent_ids = all_sample_id$id[1:32]
mapped_ids = parent_ids[map_dbl(33:nrow(all_sample_id), function(i) {
        message(i)
        id_str = all_sample_id$id[i]
        message(id_str)
        id_str_dist = Biostrings::stringDist(c(id_str, parent_ids), method = "hamming")
        id_str_dist = as.matrix(id_str_dist)[1, -1]
        mapped_id = which(id_str_dist <= 1)
        if (length(mapped_id) == 0) {
                return(NA)
        }
        if (length(mapped_id) > 1) {
                return(NA)
        }
        mapped_id
})]
names(mapped_ids) = all_sample_id$id[33:nrow(all_sample_id)]
names(parent_ids) = parent_ids

id_mapper = c(parent_ids, mapped_ids)
all_samples_unest$id_mapped = id_mapper[all_samples_unest$id]
all_samples_renest = select(all_samples_unest, -id) %>% nest(spacer_nest = -c(cell_index, id_mapped))
all_samples_renest = filter(all_samples_renest, !is.na(id_mapped))

all_samples_renest$spacer = map(all_samples_renest$spacer_nest, function(x) {
        purrr::reduce(x$spacer, merge_spacer_counts)
})

# step1 : correct N errors within each cell+id combinations
# check and repair N identities
all_samples_renest$spacer_merge = map(all_samples_renest$spacer, function(tb) {
        spacer_str = tb$spacer
        spacer_n_ind = grepl("N", spacer_str)
        if (any(spacer_n_ind) & !all(spacer_n_ind)) {
                from_indices = which(spacer_n_ind)
                to_indices = which(!spacer_n_ind)
                remove_indices = numeric()
                for (i in from_indices) {
                        for (j in to_indices) {
                                if (check_n_identity(spacer_str[i], spacer_str[j])) {
                                        message(str_c('adding ', spacer_str[i], ' to ', spacer_str[j]))
                                        tb$count[j] = tb$count[j] + tb$count[i]
                                        remove_indices = c(remove_indices, i)
                                }
                        }
                }
        } else {
                return(tb)
        }
        tb = tb[!(1:nrow(tb)) %in% remove_indices, ]
        tb
})
# step 2: correct N across cells of the same ID
# procedure: unnest correct and renest
all_samples_id_spacer = all_samples_renest %>%
        select(cell_index, id_mapped, spacer_merge) %>%
        unnest(spacer_merge) %>% nest(cell_spacer = -id_mapped)
# cell_spacer_tb = all_samples_id_spacer$cell_spacer[[1]]
all_samples_id_spacer$cell_spacer_proc = map(all_samples_id_spacer$cell_spacer,
                                             # each tb for a single id
                                             function(cell_spacer_tb) {
                                                     spacer_occurance = cell_spacer_tb %>% group_by(spacer) %>%
                                                             summarize(total = n()) %>%
                                                             arrange(desc(total))

                                                     # unique spacers with and without N
                                                     n_spacer_to_map = spacer_occurance$spacer[grepl("N", spacer_occurance$spacer)]
                                                     spacer_reference = spacer_occurance$spacer[!grepl("N", spacer_occurance$spacer)]

                                                     assertthat::assert_that(!anyDuplicated(n_spacer_to_map))
                                                     n_spacer_mapped = map_chr(n_spacer_to_map, function(x) {
                                                             map_ind = map_lgl(spacer_reference, function(a) check_n_identity(x, a))
                                                             if (all(!map_ind)) {
                                                                     return(x)
                                                             }
                                                             if (sum(map_ind) > 1) {
                                                                     warning('multiple maps')
                                                                     print(x)
                                                                     print("reference mapped: ")
                                                                     print(spacer_reference[map_ind])
                                                             }
                                                             # mapping to most abundant by default
                                                             mapped_seq = filter(spacer_occurance, spacer %in% spacer_reference[map_ind])$spacer[1]
                                                             message(str_c("mapping ", x, " to ", mapped_seq))
                                                             return(mapped_seq)
                                                     })
                                                     names(n_spacer_mapped) = n_spacer_to_map
                                                     names(spacer_reference) = spacer_reference
                                                     spacer_mapped_all = c(n_spacer_mapped, spacer_reference)
                                                     cell_spacer_tb$spacer_mapped = spacer_mapped_all[cell_spacer_tb$spacer]

                                                     cell_spacer_dup = filter(cell_spacer_tb, cell_index %in% names(which(table(cell_spacer_tb$cell_index) >1)))
                                                     assertthat::assert_that(all(map_dbl(nest(cell_spacer_dup, data = -cell_index)$data, function(x) anyDuplicated(x$spacer_mapped)) == 0))

                                                     temp_out_tb = transmute(cell_spacer_tb, cell_index = cell_index,
                                                                             spacer = spacer_mapped, count = count) %>% nest(spacer = -cell_index)
                                                     # STEP 3:filter noisy spacers
                                                     temp_out_tb$spacer = map(temp_out_tb$spacer, filter_noise_spacers)
                                                     # STEP 4: filter out counts that are fewer than cutoff: 10
                                                     print(str_c("pre filter total cells: ", nrow(temp_out_tb)))
                                                     temp_out_tb = temp_out_tb[map_lgl(temp_out_tb$spacer, function(x) max(x$count) > 5), ]
                                                     print(str_c("post filter total cells: ", nrow(temp_out_tb)))
                                                     return(temp_out_tb)

                                                     # identify components of spacers sequences within 1 hamming distance of one another
                                                     # TODO: restrict to one hamming away
                                                     # TODO: for this part, need to decide which ones are real errors by looking at if the mistache is
                                                     # close to the cut site
                                                     # spacer_unique = unique(cell_spacer_tb$spacer_mapped)
                                                     # spacer_occurance = cell_spacer_tb %>% group_by(spacer_mapped) %>%
                                                     #         summarize(total = n()) %>%
                                                     #         arrange(desc(total))
                                                     #
                                                     # spacer_dmat = Biostrings::stringDist(spacer_unique, method = "hamming")
                                                     # spacer_dmat = as.matrix(spacer_dmat)
                                                     # spacer_dmat[spacer_dmat > 1] = 0
                                                     # spacer_ig = graph.adjacency(as.matrix(spacer_dmat), mode = "undirected")
                                                     # spacer_membership = components(spacer_ig)$membership
                                                     # # now for each component, create the spacer mappings
                                                     # spacer_components = map(split(names(spacer_membership), spacer_membership),
                                                     #                         function(x) spacer_unique[as.numeric(x)])
                                                     # hamming_mapped = purrr::reduce(map(spacer_components, function(str_vec) {
                                                     #         if (length(str_vec) == 1) {
                                                     #                 names(str_vec) = str_vec
                                                     #                 return(str_vec)
                                                     #         }
                                                     #         spacer_occur_sel = filter(spacer_occurance, spacer_mapped %in% str_vec) %>% arrange(desc(total))
                                                     #         # print(spacer_occur_sel)
                                                     #         out = rep(spacer_occur_sel$spacer_mapped[1], nrow(spacer_occur_sel))
                                                     #         names(out) = spacer_occur_sel$spacer_mapped
                                                     #         out
                                                     # }), c)
                                                     # # table(cell_spacer_tb$spacer_mapped == hamming_mapped[cell_spacer_tb$spacer_mapped])
                                                     # cell_spacer_tb$spacer_mapped_hamming = hamming_mapped[cell_spacer_tb$spacer_mapped]
                                                     #
                                                     # cell_spacer_output = cell_spacer_tb %>% select(cell_index, spacer_mapped_hamming, count) %>%
                                                     #         group_by(cell_index, spacer_mapped_hamming) %>% summarize(count = sum(count))
                                                     # cell_spacer_output = cell_spacer_output %>% rename(spacer = spacer_mapped_hamming) %>%
                                                     #         nest(spacer = c(spacer, count))
                                                     #
                                                     # # next merge by noises threshold
                                                     # cell_spacer_output$spacer = map(cell_spacer_output$spacer, filter_noise_spacers)
                                                     # cell_spacer_output
                                             })
all_samples_id_spacer = all_samples_id_spacer %>% select(-cell_spacer) %>%
        unnest(cell_spacer_proc)
table(map_dbl(all_samples_id_spacer$spacer, nrow))

all_samples_id_spacer = dplyr::rename(all_samples_id_spacer, id = id_mapped)
# all_samples_unest_fil$cell_type = map_chr(strsplit(all_samples_unest_fil$cell_index, "-"), 2)
all_samples_ordered = all_samples_id_spacer %>%
        unnest(spacer) %>%
        nest(cell_spacer_tb = -id) %>%
        left_join(parent_data, by = "id") %>%
        mutate(cell_spacer_tb = map2(cell_spacer_tb, parent_spacer, function(s, p) {
                if (is.null(p)) {
                        s$unmutated = F
                        s$unmutated[which(s$spacer == spacer_unmutated_GCCAAAAGCT)] = T
                        return(s)
                }
                left_join(s, p, by ="spacer", suffix = c("", "_parent")) %>%
                        mutate(unmutated = !is.na(count_parent)) %>%
                        select(-count_parent)
        })) %>% select(-parent_spacer) %>%
        unnest(cell_spacer_tb) %>%
        mutate(spacer_mod = map2_chr(spacer, unmutated, function(a, b) {
                a[b] = "0"
                a
        })) %>%
        mutate(id_spacer = map2_chr(id, spacer_mod, function(a, b) {
                if (b == "0") {
                        return("0")
                } else {
                        return(paste0(a, "_", b))
                }
        })) %>%
        group_by(id, cell_index) %>%
        arrange(desc(count))

all_samples_ordered_nest = all_samples_ordered %>% nest(spacer_tb = -c(cell_index, id))

all_samples_ordered_nest$num_spacer = map_dbl(all_samples_ordered_nest$spacer_tb, nrow)
cell_summary = all_samples_ordered_nest %>% group_by(cell_index) %>% summarise(num_id = n(), mean_spacer = mean(num_spacer))
hist(cell_summary$mean_spacer, main = "Average spacer detected per ID")

# subset all id+cell with more than 1 spacer observed
tb_sel = all_samples_ordered_nest[map_dbl(all_samples_ordered_nest$spacer_tb, nrow) > 1, ]
# printing conflicting sequences
# tb_sel = all_samples_ordered_nest[map_dbl(all_samples_ordered_nest$spacer_tb, nrow) > 1, ]
# tb_sel = arrange(tb_sel, cell_index)
# merge n errors
# spacer_eg = tb_sel[tb_sel$cell_index == "G2T2D74" & tb_sel$id == "CGTTCATCAA", ]$spacer_tb[[1]]$spacer

# table(map_lgl(all_samples_ordered_nest$spacer, function(tb) tb$count[1] > 10),
#       map_dbl(all_samples_ordered_nest$spacer, nrow))
# table(map_dbl(all_samples_ordered_nest$spacer, nrow),
#       grepl("G1", all_samples_ordered_nest$cell_index))
# table(map_dbl(all_samples_ordered_nest$spacer, nrow),
#       grepl("G1", all_samples_ordered_nest$cell_index))

# sort(table(tb_sel$cell_index), decreasing = T)[1:50]
# for (cid in unique(all_samples_ordered_nest$cell_index)) {
#         tb_sel_cell = all_samples_ordered_nest %>% filter(cell_index == cid)
#         if (any(map_dbl(tb_sel_cell$spacer_tb_merge, nrow) > 1)) {
#           tb_sel_cell = tb_sel_cell[map_dbl(tb_sel_cell$spacer_tb_merge, nrow) > 1, ]
#           sink(paste0("temp/", cid, "_conflicting_spacers.txt"))
#           walk2(tb_sel_cell$id,
#                 tb_sel_cell$spacer_tb_merge,
#                 function(x, tb) {
#                   if (!x %in% parent_data$id) {
#                     return(NULL)
#                   }
#                   print(str_c("ID: ", x))
#                   parent_seq = parent_data$parent_spacer[[which(parent_data$id == x)]]$spacer[1]
#                   for (i in 1:nrow(tb)) {
#                     print(str_c("Spacer",  i, ": count ", tb$count[i]))
#                     print(Biostrings::pairwiseAlignment(tb$spacer[i], parent_seq))
#                   }
#                   print(str_c(rep("~", 100), collapse = ""))
#                 })
#           sink()
#         }
# }

tb_sel = all_samples_ordered_nest[map_dbl(all_samples_ordered_nest$spacer_tb, nrow) > 1, ]
tb_sel_summary = tb_sel %>%
        group_by(cell_index) %>%
        summarise(count = n()) %>%
        arrange(desc(count))
hist(tb_sel_summary$count, breaks = 40)
doublet_cells = tb_sel_summary$cell_index[tb_sel_summary$count >= 3]

# compare number of differences across different cells vs same cells
# cell_mat = cell_mat_g1[get_type(rownames(cell_mat_g1)) == "4", ]
# hist(map_dbl(1:1000, function(i) {
#   indices = sample(nrow(cell_mat), size = 2, replace = F)
#   sum(!cell_mat[indices[1], ] == cell_mat[indices[2], ], na.rm = T)
# }), xlab = "differences", main = "", breaks = 20, xlim = c(0, 20))

all_samples_freq = all_samples_ordered_nest %>% mutate(id_spacer_tb = map(spacer_tb, function(tb) {
        if (nrow(tb) == 1) {
                return(tb)
        }
        tb = tb %>% arrange(desc(count))
        tb[1, ]
}))
all_samples_freq = all_samples_freq %>% select(id, cell_index, id_spacer_tb) %>%
        unnest(cols = id_spacer_tb)
all_samples_freq = select(all_samples_freq, id, cell_index, id_spacer)

all_samples_freq = all_samples_freq %>% nest(data = -id) %>% mutate(cell_vec = map(data, function(tb) {
        tb %>%
                select(cell_index, id_spacer) %>% spread(key = cell_index, value = id_spacer)
}))
all_cells = sort(unique(purrr::reduce(map(all_samples_freq$cell_vec, colnames), c)))

# character matrix
all_cell_mat_chr = purrr::reduce(map(all_samples_freq$cell_vec, function(x) {
        x_chr = as.character(x)
        names(x_chr) = colnames(x)
        x_chr[all_cells]
}), cbind)
colnames(all_cell_mat_chr) = all_samples_freq$id
rownames(all_cell_mat_chr) = all_cells

na_cell = rowSums(is.na(all_cell_mat_chr))
cell_tb = tibble(cell_index = rownames(all_cell_mat_chr),
                 n_id = ncol(all_cell_mat_chr) - na_cell,
                 celltype = str_match(rownames(all_cell_mat_chr), pattern = "(G\\d)T(\\d)D(.*)")[, 3],
                 group = paste0("E", str_match(rownames(all_cell_mat_chr), pattern = "G(\\d)T(\\d)D(.*)")[, 2]))
mean(cell_tb$n_id)
(mutate(cell_tb, group_cell_type = paste0(group, "T", celltype)) %>%
        ggboxplot(x = "group_cell_type", y = "n_id", xlab = "") +
        ylab("# of hgRNAs detected") +
        theme(text = element_text(size = 12),
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))) %>%
                push_pdf(file_name = "id_missing_qc",
                         w = 4.5, h = 3.5, dir = plot_dir)

gghistogram(x = "detect", tibble(id = colnames(all_cell_mat_chr),
                                   detect = colMeans(!is.na(all_cell_mat_chr))),
            fill = "grey", ylab = "Count", xlab = "Fraction of cells detected") %>%
        push_pdf(file_name = "cells_missing_qc",
                 w = 2.5, h = 3.5, dir = plot_dir)

# remove doublets
all_cell_mat_chr = all_cell_mat_chr[!rownames(all_cell_mat_chr) %in% doublet_cells, ]
saveRDS(all_cell_mat_chr, file = "./intermediate_data/iPSC/cell_mat_all.rds")

cell_mat_g1 = all_cell_mat_chr[grepl("G1", get_group(rownames(all_cell_mat_chr))), ]
cell_mat_g2 = all_cell_mat_chr[grepl("G2", get_group(rownames(all_cell_mat_chr))), ]

saveRDS(cell_mat_g1, file = "./intermediate_data/iPSC/cell_mat_g1.rds")
saveRDS(cell_mat_g2, file = "./intermediate_data/iPSC/cell_mat_g2.rds")

cell_mat_g1_fil = remove_uniform_id(cell_mat_g1)
cell_mat_g2_fil = remove_uniform_id(cell_mat_g2)
cell_mat_g1_fil = filter_cells(cell_mat_g1_fil)
cell_mat_g2_fil = filter_cells(cell_mat_g2_fil)

saveRDS(cell_mat_g1_fil, file = "./intermediate_data/iPSC/cell_mat_g1_filtered.rds")
saveRDS(cell_mat_g2_fil, file = "./intermediate_data/iPSC/cell_mat_g2_filtered.rds")

cell_mat_g1_fil = readRDS("./intermediate_data/iPSC/cell_mat_g1_filtered.rds")
cell_mat_g2_fil = readRDS("./intermediate_data/iPSC/cell_mat_g2_filtered.rds")

tip_map = as.character((-1):(-6))
names(tip_map) = paste0("T", 1:6)

source("R_mod/plotting.R")
source("R_mod/im_ch.R")
source("R_mod/ice_fase_mod1.R")
cell_mat_g1_im = im_chr(cell_mat_g1_fil)
cell_mat_g2_im = im_chr(cell_mat_g2_fil)

cell_mat_g1_im_tb = as_tibble(cell_mat_g1_im)
cell_mat_g2_im_tb = as_tibble(cell_mat_g2_im)

cell_mat_g1_im_tb$cell_id = rownames(cell_mat_g1_im)
cell_mat_g2_im_tb$cell_id = rownames(cell_mat_g2_im)

cell_mat_g1_im_tb = cell_mat_g1_im_tb[c(32, 1:31)]
cell_mat_g2_im_tb = cell_mat_g2_im_tb[c(30, 1:29)]

write_csv(cell_mat_g1_im_tb, file = "./intermediate_data/iPSC/cell_mat_g1_imputed.csv")
write_csv(cell_mat_g2_im_tb, file = "./intermediate_data/iPSC/cell_mat_g2_imputed.csv")

saveRDS(cell_mat_g1_im, file = "./intermediate_data/iPSC/cell_mat_g1_imputed.rds")
saveRDS(cell_mat_g2_im, file = "./intermediate_data/iPSC/cell_mat_g2_imputed.rds")


library(furrr)
plan(multisession, workers = 12)
tr_g1 = phylotime(cell_mat_g1_im, t_total = 13. - 1., parallel = T)
tr_g2 = phylotime(cell_mat_g2_im, t_total = 13. - 1., parallel = T)
g1_celltypes = get_type(tr_g1$tip.label)
names(g1_celltypes) = tr_g1$tip.label
g2_celltypes = get_type(tr_g2$tip.label)
names(g2_celltypes) = tr_g2$tip.label

dmat_col = make_heatmap_col(max_val = 24., mid_val = 20, min_val = 15., col = rev(RColorBrewer::brewer.pal(3, "RdPu")))
dist_g1 = phylotime(cell_mat_g1_im, t_total = 13. - 1.,
                    return_dist = T,
                    parallel = T)
dmat_g1 = dist_df2mat(dist_g1)
h = Heatmap(dmat_g1[tr_g1$tip.label, tr_g1$tip.label],
            col = dmat_col,
            show_heatmap_legend = F,
            cluster_rows = as.dendrogram(tr_g1),
            cluster_columns = as.dendrogram(tr_g1),
            show_row_names = F,
            show_row_dend = F,
            show_column_names = F,
            show_column_dend = F)
push_png(h, "g1_dist_mat", dir = plot_dir, w = 3, h = 3, ps = 10, res = 1200)

dist_g2 = phylotime(cell_mat_g2_im, t_total = 13. - 1.,
                    return_dist = T,
                    parallel = T)
dmat_g2 = dist_df2mat(dist_g2)
h = Heatmap(dmat_g2[tr_g2$tip.label, tr_g2$tip.label],
            col = dmat_col,
            show_heatmap_legend = F,
            cluster_rows = as.dendrogram(tr_g2),
            cluster_columns = as.dendrogram(tr_g2),
            show_row_names = F,
            show_row_dend = F,
            show_column_names = F,
            show_column_dend = F)
push_png(h, "g2_dist_mat", dir = plot_dir, w = 3, h = 3, ps = 10, res = 1200)


h_legend = recordPlot({
        plot.new()
        draw(Legend(col_fun = dmat_col, title = "Time (days)"))
})
push_pdf(h_legend, "dist_mat_legend", dir = plot_dir)


g1_out = ice_fase_mod(tr_g1,
                      sc_celltypes = g1_celltypes,
                      total_time = 13.,
                      root_time = 1.)
g2_out = ice_fase_mod(tr_g2,
                      sc_celltypes = g2_celltypes,
                      total_time = 13.,
                      root_time = 1.)

plot_gr_data_mod1(g1_out, target_time = 14) +
        plot_gr_data_mod1(g2_out, target_time = 14)

g1_cell_type_mapped = tip_map[get_type(rownames(cell_mat_g1_fil))]
names(g1_cell_type_mapped) = rownames(cell_mat_g1_fil)
plot_barcodes(cell_mat_g1_fil,
              tr = g1_out$tr,
              tip_celltype = g1_cell_type_mapped,
              celltype_col = type_col,
              show_column_names = F) %>%
        push_pdf("g1_barcodes_v1", w = 5, h = 3.5,
                 #res = 1200,
                 dir = plot_dir)

g2_cell_type_mapped = tip_map[get_type(rownames(cell_mat_g2_fil))]
names(g2_cell_type_mapped) = rownames(cell_mat_g2_fil)
plot_barcodes(cell_mat_g2_fil,
              tr = g2_out$tr,
              tip_celltype = g2_cell_type_mapped,
              celltype_col = type_col,
              show_column_names = F) %>%
        push_png("g2_barcodes_v1", w = 5, h = 3.5,
                 res = 1200,
                 dir = plot_dir)
        # push_pdf("g2_barcodes_v1", w = 5, h = 3.5,
        #          # res = 1200,
        #          dir = plot_dir)

g1_node_map = c("Node-1" = "5", "Node-2" = "3", "Node-3" = "4", "Node-4" = "1", "Node-5" = "2")
g2_node_map = c("Node-1" = "5", "Node-2" = "4", "Node-3" = "3", "Node-4" = "1", "Node-5" = "2")
g1_node_map = c(g1_node_map, tip_map)
g2_node_map = c(g2_node_map, tip_map)
g1_col = type_col[g1_node_map]
names(g1_col) = names(g1_node_map)
g2_col = type_col[g2_node_map]
names(g2_col) = names(g2_node_map)
g1_lab = label_map[g1_node_map]
names(g1_lab) = names(g1_node_map)
g2_lab = label_map[g2_node_map]
names(g2_lab) = names(g2_node_map)

(plot_gr_data(g1_out,
              gr_col = g1_col,
              gr_lab = g1_lab, target_time = 14.) +
                plot_gr_data(g2_out,
                             gr_col = g2_col,
                             gr_lab = g2_lab, target_time = 14.)) %>%
        push_pdf("ice_fase_graph", w = 5, h = 3.5, dir = plot_dir)

# resampling cells and analyze again
sample_by_type <- function(cell_type_vec, sample_size, replace = F) {
        purrr::reduce(map(split(1:length(cell_type_vec), cell_type_vec), function(x) {
                sample(x, size = sample_size, replace = replace)
        }), c)
}
rs_tb = tibble(expand.grid(sim_id = 1:50,
                           group = as.character(1:2),
                           sample_size = c(60, 80, 100, 120, 140)))

table(get_type(rownames(cell_mat_g1_fil)))
table(get_type(rownames(cell_mat_g2_fil)))

j <- 0
rs_tb = mutate(rs_tb, data = map2(group, sample_size, function(g, ss) {
        j <<- j + 1
        set.seed(j)
        message(j)
        if (g == "1") {
                sample_indices = sample_by_type(get_type(rownames(cell_mat_g1_fil)), sample_size = ss)
                out = ice_fase(cell_mat_g1_fil[sample_indices, ],
                               sc_celltypes = get_type,
                               total_time = 13.,
                               root_time = 1.)
        }
        if (g == "2") {
                sample_indices = sample_by_type(get_type(rownames(cell_mat_g2_fil)), sample_size = ss)
                out = ice_fase(cell_mat_g2_fil[sample_indices, ],
                               sc_celltypes = get_type,
                               total_time = 13.,
                               root_time = 1.)
        }
        return(out)
}))
# saveRDS(rs_tb, file = "../LTModelData/rs_final.rds")
saveRDS(rs_tb, file = "./intermediate_data/iPSC/ice_fase_resampled.rds")

# reading mutation rate data
mr_tb = readxl::read_excel("./data/iPSC//Lin96_hgRNA_mutation_rates.xlsx")
t_vec = c(4, 8, 11)
mr_tb$mr = map_dbl(1:nrow(mr_tb), function(i) {
        f_vec = as.numeric(mr_tb[i, c("D4", "D8", "D11")]) / 100
        mean(-log(1 - f_vec + 1e-10) / t_vec)
})
write_tsv(mr_tb[c(1, 6)], "./supplementary_data/supplementary_data_2_mutation_rates/iPSC_mutation_rates.tsv")

mr_vec = mr_tb$mr[match(colnames(cell_mat_g1), mr_tb$hgRNA)]
mr_sc_est = estimate_mut_rate(cell_mat_g1, 13.)

library(circlize)
col_func = colorRamp2(breaks = c(-3.2, -2.2, -1, 0, 0.5), colors = c( "#325288", "#325288", "#ffdf6b", "#a20a0a", "#a20a0a"))
id_col = col_func(log10(pmax(mr_tb$mr, 1e-4)))
names(id_col) = mr_tb$hgRNA
mr_tb = filter(mr_tb, !(hgRNA %in% c("CGAAAATAGT", "GCAAGTGAGT", "TGGTCAACAA", "TTAATAGACC")))

# Supplementary Figure
gather(mr_tb[1:5], key = "Day", value = "Percent", -c(hgRNA)) %>%
        ggline(x = "Day", y = "Percent", color = "hgRNA") +
        scale_color_manual(values = id_col) +
        theme(legend.position = "right")

plot(mr_sc_est[1:31], mr_vec[1:31])
abline(0, 1)


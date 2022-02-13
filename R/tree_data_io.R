# dump_sc_data <- function(res_data, target_dir, num_elements = NULL) {
#         if (!dir.exists(target_dir)) {
#                 dir.create(target_dir)
#         }
#         if (is.null(num_elements)) {
#                 message("using all elements..")
#                 num_elements = ncol(res_data[[1]][[1]]$sc)
#         }
#         # write(recur_vec, file = file.path(target_dir, "allele_prob.txt"), ncolumns = 1)
#         # write(mut_rate_vec, file = file.path(target_dir, "mut_rate.txt"), ncolumns = 1)
#         for (j in 1:length(res_data)) {
#                 data_dir = file.path(target_dir, paste0("sim_", j))
#                 data_dir_sc = file.path(data_dir, "sc")
#                 # data_dir_bulk = file.path(data_dir, "bulk")
#                 if (!dir.exists(data_dir_sc)) {
#                         dir.create(data_dir_sc, recursive = T)
#                 }
#                 # if (!dir.exists(data_dir_bulk)) {
#                 #         dir.create(data_dir_bulk, recursive = T)
#                 # }
#                 for (i in 1:length(res_data[[j]])) {
#                         barcode_df = data.frame(res_data[[j]][[i]]$sc[, 1:num_elements])
#                         barcode_df = barcode_df[, !apply(barcode_df, 2, function(z) all(z == "0"))]
#                         barcode_df$cell_label = rownames(res_data[[j]][[i]]$sc)
#                         readr::write_tsv(barcode_df,
#                                          path = file.path(data_dir_sc, paste0(stringr::str_pad(i, width = 3, pad = "0"), "_barcodes.txt")),
#                                          col_names = F)
#                         ape::write.tree(res_data[[j]][[i]]$tr, file.path(data_dir_sc, paste0(stringr::str_pad(i, width = 3, pad = "0"), "_tr.txt")))
#                         # saveRDS(res_data[[j]][[i]]$bulk, file = file.path(data_dir_bulk, paste0(stringr::str_pad(i, width = 3, pad = "0"), ".rds")))
#                 }
#         }
# }
read_cassiopeia_tree <- function(data_dir) {
        recon_files = list.files(data_dir, "*_recon.txt", full = T)
        out = purrr::map(recon_files, function(y) {
                recon = ape::read.tree(y)
                recon = collapse.singles(recon)
                recon
        })
        names(out) = recon_files
        out
}
load_sc_recon <- function(data_dir) {
        dirs = list.dirs(data_dir, recursive = F)
        dirs0 = paste0(data_dir, "/sim_", 1:length(dirs))
        assertthat::assert_that(all(dirs %in% dirs0))
        out = purrr::map(dirs0, function(x) {
                x = paste0(x, "/sc/")
                recon_files = list.files(x, "*_recon.txt", full = T)
                purrr::map(recon_files, function(y) {
                        recon = ape::read.tree(y)
                        recon$edge.length = NULL
                        recon = collapse.singles(recon)
                        recon
                })
        })
        out

}
dump_trees <- function(tree_list, output_dir, suffix = "") {
        for (i in 1:length(tree_list)) {
                if (!is.null(tree_list[[i]])) {
                        ape::write.tree(tree_list[[i]],
                                        file.path(output_dir, paste0(stringr::str_pad(i, width = 4, pad = "0"), suffix, ".txt")))
                }
        }
}
dump_mut_p <- function(mut_p, output_dir) {
        for (i in 1:length(mut_p$recur_vec_list)) {
                write(mut_p$recur_vec_list[[i]], file = file.path(output_dir,
                                                                  paste0("element", stringr::str_pad(i, width = 4, pad = "0"),
                                                                         "_allele_prob.txt")), ncolumns = 1)
        }
        write(mut_p$mut_rate, file = file.path(output_dir, "mut_rate.txt"), ncolumns = 1)
}
dump_mut_p_list <- function(mut_p_list, output_dir) {
        for (j in 1:length(mut_p_list)) {
                mutp_dir = paste0(output_dir, "exp_", stringr::str_pad(j, width = 4, pad = "0"), "/")
                dir.create(mutp_dir)
                dump_mut_p(mut_p_list[[j]], mutp_dir)
        }
}
dump_sc_mat <- function(sc_mat_list, output_dir) {
        for (i in 1:length(sc_mat_list)) {
                sc_mat = sc_mat_list[[i]]
                barcode_df = data.frame(sc_mat)
                barcode_df = barcode_df[, !apply(barcode_df, 2, function(z) all(z == "0"))]
                barcode_df$cell = rownames(sc_mat)
                readr::write_tsv(barcode_df,
                                 path = file.path(output_dir, paste0(stringr::str_pad(i, width = 4, pad = "0"), "_barcodes.txt")),
                                 col_names = F)
        }
}

# original id spacer
table_s1 = read_csv('./metadata/MARC1/TableS1.csv')
table_s3 = read_csv("./metadata/MARC1/SupplementaryTable3 - TableS3-hgRNA Activities.csv")
# begin updated version of processing functions
load_count_file <- function(fn) {
        tb = read_tsv(fn, col_names = F)
        tb = tb[, 1:3] %>% rename(id = X1, spacer = X2, count = X3) %>% nest(spacer = c(spacer, count))
        tb
}
label_unmutated <- function(samples, parent_data) {
        stop("this function is deprecated, use assign_unmutated instead")
}
assign_unmutated <- function(data_tb, parent_data) {
        data_tb %>% left_join(parent_data, by = "id") %>%
                mutate(spacer = map2(spacer, parent_spacer, function(spacer_tb, parent_tb) {
                        if (is.null(parent_tb)) {
                                return(spacer_tb %>% mutate(unmutated = F))
                        }
                        left_join(spacer_tb, parent_tb, by ="spacer", suffix = c("", "_parent")) %>%
                                mutate(unmutated = !is.na(count_parent)) %>%
                                select(-count_parent)
                })) %>% select(-parent_spacer)
}
# end updated version of processing functions

table_s1 = table_s1[1:112, ]
table_s1_active = table_s1[table_s1$`Identifier (ID)` %in% table_s3$identifier[table_s3$Class != "inactive"], ]
active_id = table_s1_active$`Identifier (ID)`
inactive_id = table_s3$identifier[table_s3$Class == "inactive"]
id_pb3 = table_s1$`Identifier (ID)`[table_s1$Founder == "PB3"]
id_pb7 = table_s1$`Identifier (ID)`[table_s1$Founder == "PB7"]
id_spacer_df_pb3 = data.frame(id = id_pb3[!is.na(id_pb3)])
id_spacer_df_pb3$spacer = table_s1$`Spacer sequence (TSS to PAM)`[match(id_spacer_df_pb3$id, table_s1$`Identifier (ID)`)]
id_spacer_df_pb7 = data.frame(id = id_pb7[!is.na(id_pb7)])
id_spacer_df_pb7$spacer = table_s1$`Spacer sequence (TSS to PAM)`[match(id_spacer_df_pb7$id, table_s1$`Identifier (ID)`)]

#' read a list of files, order by id then sample
process_files <- function(dir0, f0) {
        sample_df_list = lapply(paste0(dir0, f0), function(fn) {
                read_tsv(fn, col_names = F)
        })
        all_id_seq = unique(unlist(lapply(sample_df_list, function(x) x$X1)))
        id_sample_list = lapply(all_id_seq, function(id_seq) {
                out = lapply(sample_df_list, function(x) {
                        x[x$X1 == id_seq, ]
                })
                names(out) = f0
                out
        })
        names(id_sample_list) = all_id_seq
        id_sample_list
}
#' get all unique spacer sequences for each id from all samples
collect_all_spacers <- function(id_spacer_count_list, parent_list) {
        all_id_seq = names(id_spacer_count_list)
        spacer_df_all_id = lapply(1:length(all_id_seq), function(i) {
                id_seq = all_id_seq[i]
                parent_seq = parent_list[[id_seq]][[1]]$X2
                # spacer_seq = id_spacer_df$spacer[id_spacer_df$id == id_seq]
                spacer_sample_list = id_spacer_count_list[[i]]
                all_spacer_seq = unique(unlist(lapply(spacer_sample_list, function(x) x$X2)))
                # if (sum(all_spacer_seq %in% parent_seq) != 1) {
                #         print(sum(all_spacer_seq %in% parent_seq))
                # }
                spacer_df = data.frame(spacer = all_spacer_seq,
                                       unmutated = all_spacer_seq %in% parent_seq)
                spacer_df$spacer = as.character(spacer_df$spacer)
                spacer_df
        })
        names(spacer_df_all_id) = all_id_seq
        spacer_df_all_id
}
#' get spacer count matrix for each id
make_spacer_mat <- function(id_spacer_count_list, all_spacer_list) {
        all_id_seq = names(id_spacer_count_list)
        spacer_df_all_id = lapply(1:length(all_id_seq), function(i) {
                spacer_sample_list = id_spacer_count_list[[i]]
                spacer_mat = do.call(rbind, lapply(spacer_sample_list, function(x) {
                        all_spacer_seq = all_spacer_list[[i]]$spacer
                        sapply(all_spacer_seq, function(spacer_seq) {
                                if (spacer_seq %in% x$X2) {
                                        # should be only one element, but sum for exceptions for now
                                        return(sum(x$X3[x$X2 == spacer_seq]))
                                } else {
                                        return(0)
                                }
                        })
                }))
                spacer_mat
        })
}



file_list<-list.files("./../6state_exp_test_data_v1/trees/", pattern= "cassiopeia_v\\d+__resolved10000.txt*")
tree_nums <-data_frame(as.numeric(str_replace(file_list, "_cassiopeia_v\\d+__resolved10000.txt", "")))
colnames(tree_nums) <- "tree"
tree_nums$chain <- as.numeric(str_replace(str_replace(file_list, "\\d+_cassiopeia_v", ""), "__resolved10000.txt", ""))
tree_nums$mode <- "cassiopeia"
tree_nums$tr<- paste0("./../6state_exp_test_data_v1/trees/", file_list) %>% map(ape::read.tree)
unique(tree_nums$tree)
write_rds(tree_nums, "datasets/cassiopeia_compiled_with_bl.rds")
file_list<-list.files("./../6state_exp_test_data_v1/trees/", pattern= "__phylip.txt*")

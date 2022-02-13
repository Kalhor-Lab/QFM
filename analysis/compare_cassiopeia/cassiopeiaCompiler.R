runtime_dir <- "../Python/outputs/cassiopeia_longer/"
mappings <- read.csv("mapping.csv")
df2<- read.csv("../Python/outputs/cassiopeia_longer/6state_exp_test_cassi_run_time_2.txt",header = F)
colnames(df2) <- c("index", "rnaCount", "runtime")
df2 <- merge(df2,mappings, by.x = "index", by.y = "tree_num")
df2$nnode <- 1
df2$hgrnanum <- 2 
for (i in 1:nrow(df2)){
  tree_num <- df2[i,]$qreenum
  df2[i,]$nnode <- exp_params_test$sample_size[[i]] * 6
  df2[i,]$hgrnanum <- exp_params$num_element[[i]]
}
values 
treechain <- df2$index 
cassiopeia0121_recontructed.txt
fileList<-list.files(path = "../Python/outputs", pattern = "*_100_recontructed.txt")
index[1]
files <- paste0("6state_exp_test_cassi_run_time_bak",2,".txt")
fiftyfiles<-list.files(path=runtime_dir, pattern = "6state_exp_test_cassi_run_time_bak_50_")
index <- as.numeric(str_remove(str_remove(fiftyfiles, pattern = "6state_exp_test_cassi_run_time_bak_50_"), pattern=".txt" ))
df2$hgrnaNum <- 100
df2$runtime <- 0
df<-as.data.frame(index)
for (i in 1:nrow(df2)){
  tree_num <- df2$index[[i]]
  df2[i,]$runtime <- read.csv(paste0(runtime_dir,"6state_exp_test_cassi_run_time_bak",tree_num,".txt"), header=F)
}
read.csv("../Python/outputs/cassiopeia_longer/6state_exp_test_cassi_run_time_bak_50_1.txt")
"../Python/outputs"
list.files(runtime_dir)
index <- as.numeric(str_remove(str_remove(fileList, pattern = "cassiopeia_longer"), pattern="_100_recontructed.txt" ))
tree_chain_files <- paste0("../Python/outputs/cassiopeia_longer" ,str_pad(index, 4, pad = "0"), "_", 100 ,"_recontructed.txt")
df2 <- as.data.frame(index)
df2$tree <- paste0(tree_chain_files) %>% map(ape::read.tree)
df2<- df2 %>% select(!qreenum)
df2$mode <- "Cassiopeia"
write_rds(df, "Cass50SupRuns.Rds")
write_rds(df2, "Cass100SupRuns.Rds")


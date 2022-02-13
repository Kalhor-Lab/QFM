## Breaks the trees that have 1600 tips and 100 barcodes into a 100 barcode, 50 barcode, and 25 barcode tree
# then saves that tree
library(tidyverse)
exp_params<- readRDS("exp_data_rep1_2_proc.rds")
reduced <- exp_params %>% filter(num_tip == 16)
mapping <- reduced %>% select(j)
mapping$tree_num<- 1:nrow(mapping)
write.csv(mapping, "mapping.csv")
dirpath<- "../6state_exp_test_data_v2/data/"
i<- 1 
for (j in 1:nrow(reduced)){
  write.table(as.data.frame(matrix(reduced$sc[[j]][,1:25], ncol = 25)) %>% 
                mutate(names = rownames(reduced$sc[[j]])),
              paste0(dirpath, str_pad(j, 4, pad = "0"),"_25_barcodes.txt"), col.names = F, row.names = F, sep = "\t", quote = F)
  write.table(as.data.frame(matrix(reduced$sc[[j]][,1:50], ncol = 50)) %>% 
                mutate(names = rownames(reduced$sc[[j]])),
              paste0(dirpath, str_pad(j, 4, pad = "0"),"_50_barcodes.txt"), col.names = F, row.names = F, sep = "\t", quote = F)
  write.table(as.data.frame(matrix(reduced$sc[[j]][,1:100], ncol = 100)) %>% 
                mutate(names = rownames(reduced$sc[[j]])),
              paste0(dirpath, str_pad(j, 4, pad = "0"),"_100_barcodes.txt"), col.names = F, row.names = F, sep = "\t", quote = F)
                                      i<- i+1
}

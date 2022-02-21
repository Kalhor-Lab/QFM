## ----setup, message=F, include=FALSE------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(qfm)

## ---- message=FALSE-----------------------------------------------------------
# read mutant allele matrix
cell_allele_table = readr::read_csv("../data/example/mutant_allele_matrix.csv")
cell_allele_table[1:5]

## -----------------------------------------------------------------------------
chr_mat = as.matrix(cell_allele_table[-c(1:2)])
rownames(chr_mat) = cell_allele_table$cell

library(furrr)
plan(multisession, workers = 8)
tr = phylotime(chr_mat[, 1:50], t_total = 15. - 0.6)

## ---- fig.width=8, fig.height=5-----------------------------------------------
plot_barcodes(chr_mat[, 1:50], tr, show_column_names = F)

## ---- fig.width=10, fig.height=5----------------------------------------------
plot_tr(tr,
        end_alpha_terminal = 0.05,
        edge_width_terminal = 0.02)

## ---- fig.width=8, fig.height=5-----------------------------------------------
sc_celltypes = cell_allele_table$type
names(sc_celltypes) = cell_allele_table$cell
print(sc_celltypes[1:10])
plot_barcodes(chr_mat[, 1:50], tr, tip_celltype = sc_celltypes, show_column_names = F)

## -----------------------------------------------------------------------------
res = ice_fase(tr, sc_celltypes, total_time = 15 - 0.6, root_time = 0.6)

## -----------------------------------------------------------------------------
res = set_color_palette(res)

## ---- fig.width=10, fig.height=5----------------------------------------------
plot_gr(gr = res$gr,
        total_time = 15,
        gr_node_time = res$gr_trans_time,
        type_col = res$col_pal)

## ---- fig.width=10, fig.height=5----------------------------------------------
plot_tr(res$tr,
        node_types = res$tr_node_assign,
        type_col = res$col_pal,
        end_alpha_terminal = 0.05,
        edge_width_terminal = 0.02)

## ---- fig.width=8, fig.height=5-----------------------------------------------
plot_ice_times(res)

## ---- fig.width=8, fig.height=5-----------------------------------------------
plot_node_sizes(res)

## -----------------------------------------------------------------------------
output_estimates(res)


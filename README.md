# Quantitative Fate Mapping (QFM) via retrospective lineage barcodes

Quick start:
To install:
```
devtools::install_github("Kalhor-Lab/QFM")
```
To run ICE_FASE:
```
result = ice_fase(cell_mat,
                  sc_celltypes,
                  total_time,
                  root_time)
```
To run Phylotime:
```
tr = phylotime(cell_mat,
               total_time)
```

To reproduce results in manuscript:
  1. MARC1 estimations ('analysis/MARC1/')
  2. Evaluation of ICE_FASE and Phylotime on simulation of panel of fate maps ('analysis/panel/')
  3. Process and analyze split well experiment ('analysis/splitwell/')

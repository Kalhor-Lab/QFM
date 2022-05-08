# Quantitative Fate Mapping (QFM) via retrospective lineage barcodes

To install:
```
if (!require("devtools"))
  install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

# install two dependencies from Bioconductor
BiocManager::install("Biostrings")
BiocManager::install("ComplexHeatmap")

devtools::install_github("Kalhor-Lab/QFM")
```

For a basic tutorial, check out the vignette:
https://kalhor-lab.github.io/QFM/

To reproduce results in manuscript:
  1. MARC1 estimations ('analysis/MARC1/')
  2. Evaluation of ICE_FASE and Phylotime on simulation of panel of fate maps ('analysis/panel/')
  3. Process and analyze split well experiment ('analysis/splitwell/')

Data from the manuscript can be found at:
https://github.com/Kalhor-Lab/QFM-Data

Report bugs and provide suggestions by sending email to:

Author and maintainer: Weixiang Fang (wfang9@jh.edu)

Or open a new issue on this Github page

Corresponding authors:
Reza Kalhor, Hongkai Ji

Preprint:
https://biorxiv.org/cgi/content/short/2022.02.13.480215v1

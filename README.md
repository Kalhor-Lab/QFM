# Quantitative Fate Mapping (QFM) via retrospective lineage barcodes

To install:
```
if (!require("devtools"))
  install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

# install dependencies from Bioconductor
BiocManager::install("ComplexHeatmap")

devtools::install_github("Kalhor-Lab/QFM")
```

For a basic tutorial, check out the vignette:
https://kalhor-lab.github.io/QFM/articles/basic_tutorial.html

To reproduce results in manuscript: Check 'analysis/' folder.

Data from the manuscript can be found at:
https://github.com/Kalhor-Lab/QFM-Data

Report bugs and provide suggestions by sending email to:

Author and maintainer: Weixiang Fang (wfang9@jh.edu)

Or open a new issue on this Github page

Corresponding authors:
Reza Kalhor, Hongkai Ji

Preprint:
https://biorxiv.org/cgi/content/short/2022.02.13.480215v1

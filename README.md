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
<a href="https://doi.org/10.5281/zenodo.7112097"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7112097.svg" alt="DOI"></a>

Report bugs and provide suggestions by sending email to:

Author and maintainer: Weixiang Fang (wfang9@jh.edu)

Or open a new issue on this Github page

Corresponding authors:
Reza Kalhor, Hongkai Ji

## Citation
Fang, Weixiang, et al. "Quantitative fate mapping: A general framework for analyzing progenitor state dynamics via retrospective lineage barcoding." Cell 185.24 (2022): 4604-4620.

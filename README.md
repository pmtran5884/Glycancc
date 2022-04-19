# Glycancc
Package to analyze anti-carbohydrate antibody population data

# To Install
The packages "[PCAtools](https://bioconductor.org/packages/release/bioc/html/PCAtools.html)" and "[ComplexHeatmap](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)" have to be installed seperately from Bioconductor first before package installation.
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

BiocManager::install("PCAtools")
```
Now install the "Glycancc" package
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("pmtran5884/Glycancc")
```

# Usage
See "code for figures" folder for example usage

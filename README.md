# Glycancc
Package to analyze anti-carbohydrate antibody population data

# To Install (total time < 5min)
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

# Usage (time <5 min)
See "code for figures" folder for detailed code for generating manuscript figures

### Load library and set seed
```
library(Glycancc)

####### set seed ##########
set.seed(123)
```

### Generate K-nearest neighbors network and louvain cluster ACAs
```
#umap
daisy_umap<-make_umap(t(daisy$igg),n_neighbors = 5,metric = "cosine")

#network
daisy_knn_graph<-knn_network_from_umap(daisy_umap,metric="cosine")

#plot
plot_cluster_graph(daisy_knn_graph, daisy_umap, cluster_method = "louvain")
```

### Linear models of ACA clusters against type 1 diabetes related phenotypes
```
#calculate the first principal component for each ACA cluster
daisy$eigenmatrix<-make_eigen_matrix(daisy$igg,co_cluster)

#Combine eigenmatrix with daisy phenotype data and add to daisy list and reorder factor levels
daisy$combined<-cbind.data.frame(daisy$eigenmatrix,daisy$pheno)
daisy$combined$Group<-factor(daisy$combined$Group,levels = c("Control","Non-progressor","Progressor"))

#Calculate linear models for each ACA cluster X ~ Group + Sex + FDR + HLA_risk + Draw_Age
daisy_eigen_models<-univariate_models(daisy$combined,
                                      c("Group", 
                                        "Sex", "FDR", 
                                        "HLA_risk","Draw_Age"))
                                        
#Convert linear models results to ggplot2 friendly format                                        
daisy_model<-convert_to_ggplotformat(daisy_eigen_models$fulltable)

# forest plot for all ACA clusters
forest_plot_lm_model(daisy_model)
```

# End
For more information, please read the preprint on Research Square "[The Anti-carbohydrate antibody repertoire in type 1 diabetes](https://assets.researchsquare.com/files/rs-1490184/v1_covered.pdf?c=1649958275)"

#Citation
If you use the package and its code, please cite the release from Zenodo (DOI: 10.5281/zenodo.7143430).

#Glycan class enrichment in clusters
#Date: March 9, 2022
#Author: Paul Tran
#Updated: March 9, 2022
rm(list = ls())


########## HELPER FUNCTIONS #####################
source("C:/Users/spurohit/Box/Sharad_DataShare/Papers/Students/T1D glycan array data/aca_rep_t1d/scripts/helper_functions_v2.R")

########## Libraries #############
# BiocManager::install("PCAtools")
# devtools::install_github("jokergoo/ComplexHeatmap")

library(pacman)
library(PCAtools)
library(ComplexHeatmap)
p_load("xlsx")
p_load("umap")
p_load("igraph")
p_load("ggplot2")
p_load("alluvial")
p_load("Hmisc")
p_load("fastDummies")
p_load("dplyr")
p_load("gmodels")
p_load("gtsummary")
p_load("dendextend")
p_load("circlize")
p_load("Ryacas")
p_load("PCAtools")
p_load("BiocManager")
p_load("pROC")
p_load("fmsb")
p_load("glmnet")

####### set seed ##########
set.seed(123)

library(Glycancc)
#cluster enrichment analysis
cluster_class_enrichment<-glycan_class_enrichment(co_cluster,
                                                  glycan_classes)

plot_significant_cluster_classes(cluster_class_enrichment)

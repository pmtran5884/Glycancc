#Graph clustering 
#Date: March 9, 2022
#Author: Paul Tran
#Updated: March 9, 2022
rm(list = ls())

########## HELPER FUNCTIONS #####################
library(Glycancc)

####### set seed ##########
set.seed(123)

######## LOAD AND cLEAN DATA ################
daisy<-load_and_clean_DAISY_data(dataloc="c:/Users/spurohit/Box/Sharad_DataShare/Papers/Students/T1D glycan array data/aca_rep_t1d/raw data/1b_NGCBgMFI_Pheno data for_IgG_DAISY.csv",
                                 phenoloc="c:/Users/spurohit/Box/Sharad_DataShare/Papers/Students/T1D glycan array data/aca_rep_t1d/raw data/1_DAISY112_phenotype data.xls",
                                 removeoutliers = F, 
                                 backgroundThesholdval = 1,
                                 log2=T)

fileloc="c:/Users/spurohit/Box/Sharad_DataShare/Papers/Students/T1D glycan array data/aca_rep_t1d/raw data/1b_NGCBgMFI_Pheno data for_IgG_T1D samples_20220222.csv"
pagoda<-load_and_clean_PAGODA_data(fileloc,removeoutliers = F, 
                                   backgroundThesholdval = 1,
                                   log2=T)



#umap
daisy_umap<-make_umap(t(daisy$igg),n_neighbors = 5,metric = "cosine")
pagoda_umap<-make_umap(t(pagoda$igg),n_neighbors = 5,metric = "cosine")

#network
daisy_knn_graph<-knn_network_from_umap(daisy_umap,metric="cosine")
pagoda_knn_graph<-knn_network_from_umap(pagoda_umap,metric="cosine")

#plot
plot_cluster_graph(daisy_knn_graph, daisy_umap, cluster_method = "louvain")
plot_cluster_graph(pagoda_knn_graph, pagoda_umap, cluster_method = "louvain")

#merge networks
merged_graph<-join_graphs(daisy_knn_graph,pagoda_knn_graph)

#cluster
merged_clusters<-cluster_graph(merged_graph,cluster_method = "louvain")

#plot merged cluster
plot(merged_clusters,merged_graph,layout=igraph::layout_with_drl(merged_graph),vertex.size=10, vertex.label=NA)


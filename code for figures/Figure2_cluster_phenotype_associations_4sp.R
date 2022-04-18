#Linear associations of glycan clusters with phenotypic variables
#Date: March 9 2022
#Author: Paul Tran
#Updated: March 9, 2022
rm(list = ls())

library(Glycancc)

####### set seed ##########
set.seed(123)

######## LOAD AND cLEAN DATA ################
daisy<-load_and_clean_DAISY_data(dataloc="c:/Users/spurohit/Box/Sharad_DataShare/Papers/Students/T1D glycan array data/aca_rep_t1d/raw data/1b_NGCBgMFI_Pheno data for_IgG_DAISY.csv",
                                 phenoloc="c:/Users/spurohit/Box/Sharad_DataShare/Papers/Students/T1D glycan array data/aca_rep_t1d/raw data/1_DAISY112_phenotype data.xls",
                                 removeoutliers = F, 
                                 backgroundThesholdval = 1,
                                 log2=T)

co_cluster<-readRDS("c:/Users/spurohit/Box/Sharad_DataShare/Papers/Students/T1D glycan array data/aca_rep_t1d/raw data/glycan_cluster_membership.rds")

#eigenvalue
daisy$eigenmatrix<-make_eigen_matrix(daisy$igg,co_cluster)

#lm
daisy$combined<-cbind.data.frame(daisy$eigenmatrix,daisy$pheno)
daisy$combined$Group<-factor(daisy$combined$Group,levels = c("Control","Non-progressor","Progressor"))
daisy_eigen_models<-univariate_models(daisy$combined,
                                      c("Group", 
                                        "Sex", "FDR", 
                                        "HLA_risk","Draw_Age"))
#check lm assumptions
# par(mfrow=c(2,2))
# plot(daisy_eigen_models$models$Cluster7)
daisy_model<-convert_to_ggplotformat(daisy_eigen_models$fulltable)

# forest plot for all ACA clusters
forest_plot_lm_model(daisy_model)

#ACA predictive models
#Date: March 9, 2022
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


#predictive mod
daisy<-readRDS("c:/Users/spurohit/Box/Sharad_DataShare/Papers/Students/T1D glycan array data/aca_rep_t1d/raw data/DAISY_with_eigenvals.rds")
fullmod <- glm(Group ~ 
                 Sex + 
                 FDR +
                 HLA_risk +
                 Draw_Age +
                 Cluster1 +
                 Cluster2 +
                 Cluster3 +
                 Cluster4 +
                 Cluster5 +
                 Cluster6 +
                 Cluster7 +
                 Cluster8 +
                 Cluster9 +
                 Cluster10 +
                 Cluster11, 
               data = daisy$combined,
               subset = which(daisy$combined$Group!="Non-progressor"),
               na.action = na.omit,
               family = "binomial")

cl8mod <- glm(Group ~ 
                Sex + 
                FDR +
                HLA_risk +
                Draw_Age +
                Cluster1+
                Cluster2+
                Cluster3+
                Cluster5+
                Cluster7+
                Cluster8, 
              data = daisy$combined,
              subset = which(daisy$combined$Group!="Non-progressor"),
              na.action = na.omit,
              family = "binomial")

nullmod <- glm(Group ~ 
                 Sex + 
                 FDR +
                 HLA_risk +
                 Draw_Age, 
               data = daisy$combined,
               subset = which(daisy$combined$Group!="Non-progressor"),
               na.action = na.omit,
               family = "binomial")

#LRT to compare models

lmtest::lrtest(fullmod, nullmod)
lmtest::lrtest(cl8mod, nullmod)
lmtest::lrtest(cl8mod, fullmod)

#ROC plot
test_roc<-plot(smooth(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                          predict(nullmod),ci=T)),print.auc=T,col="black")
# plot(ci.se(test_roc),type="shape",col="black")
test_roc<-plot(smooth(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                          predict(fullmod),ci=T)),print.auc=T,add=T,col="red",print.auc.y=0.4)
test_roc<-plot(smooth(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                          predict(cl8mod),ci=T)),print.auc=T,add=T,col="blue",print.auc.y=0.3)

#partial r-squared
myrsq<-data.frame(rsq::rsq.partial(fullmod,adj = T))
myrsq2<-rbind.data.frame(myrsq[1:4,],
                         c("TRUE","GLYCANS",sum(myrsq$partial.rsq[grep("Cluster",myrsq$variable)])))
myrsq2$partial.rsq<-as.numeric(myrsq2$partial.rsq)

ggplot2::ggplot(myrsq2,ggplot2::aes(x=partial.rsq,y=variable))+ggplot2::geom_col()+ggplot2::theme_classic()


#AUC saturation analysis
daisy_prog<-list()
daisy_prog$igg<-daisy$igg[daisy$pheno$Group!="Non-progressor",]
daisy_prog$pheno<-daisy$pheno[daisy$pheno$Group!="Non-progressor",]

#random gly ridge model
auc_df<-data.frame("no_gly"=0,"mean"=0,"lowerCI"=0,"upperCI"=0)[-1,]

for (i in seq(1,40,10)){
  x <- replicate(10, {glmnet_auc(daisy_prog,i)})
  auc_df<-rbind.data.frame(auc_df,
                           c(i,mean(x),mean(x)-1.96*sd(x),mean(x)+1.96*sd(x)))
}
colnames(auc_df)<-c("no_gly","AUC","lowerCI","upperCI")

write.csv(auc_df,"/Users/paultran/Downloads/Gly_saturation_aucs.csv")

ggplot2::ggplot(data=auc_df,ggplot2::aes(x=no_gly,y=AUC))+
  ggplot2::geom_point()+
  ggplot2::geom_errorbar(ggplot2::aes(ymin=lowerCI,ymax=upperCI))+
  ggplot2::theme_classic()#+
# geom_function(fun=function(x)log(x,3))
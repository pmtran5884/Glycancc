#ACA predictive models
#Date: March 9, 2022
#Author: Paul Tran
#Updated: October 4, 2022

rm(list = ls())

library(Glycancc)

#calculate the first principal component for each ACA cluster
daisy$eigenmatrix<-make_eigen_matrix(daisy$igg,co_cluster)

#Combine eigenmatrix with daisy phenotype data and add to daisy list and reorder factor levels
daisy$combined<-cbind.data.frame(daisy$eigenmatrix,daisy$pheno)
daisy$combined$Group<-factor(daisy$combined$Group,levels = c("Control","Non-progressor","Progressor"))

#predictive mod
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

#randmod function
randmodAUC<-function(mydata=daisy$combined){
  mydata$Grouprand1<-sample(mydata$Group)
  randmod1<-glm(Grouprand1 ~
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
                data = mydata,
                subset = which(mydata$Grouprand1!="Non-progressor"),
                na.action = na.omit,
                family = "binomial")

  #ROC plot
  test_roc<-pROC::roc(fastDummies::dummy_cols(mydata$Grouprand1[mydata$Grouprand1!="Non-progressor"])[,4],
                      predict(randmod1))
  as.numeric(test_roc$auc)
}


#Calculate ROC for models
null_roc<-plot(pROC::smooth(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                         predict(nullmod),ci=T)),print.auc=T,add=F,col="black")
full_roc<-plot(pROC::smooth(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                                predict(fullmod),ci=T)),print.auc=T,add=T,col="red",print.auc.y=0.4)
cl8_roc<-plot(pROC::smooth(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                                predict(cl8mod),ci=T)),print.auc=T,add=T,col="blue",print.auc.y=0.3)

#Calculate q val for each model
AUCs<-replicate(10000,randmodAUC(mydata=daisy$combined))

qval<-function(randAUCs=AUCs,testAUC=full_roc$auc){
qval<-1-(sum(randAUCs<as.numeric(testAUC))/10000)
qval
}

qval(testAUC = full_roc$auc)
qval(testAUC = cl8_roc$auc)
qval(testAUC = null_roc$auc)

#LRT to compare models
lmtest::lrtest(fullmod, nullmod)
lmtest::lrtest(cl8mod, nullmod)
lmtest::lrtest(cl8mod, fullmod)

# test_roc<-plot(pROC::smooth(pROC::roc(fastDummies::dummy_cols(daisy$combined$Grouprand1[daisy$combined$Grouprand1!="Non-progressor"])[,4],
#                          predict(randmod1),ci=T)),print.auc=T,col="grey")


############### Figure 5b #############
#partial r-squared
myrsq<-data.frame(rsq::rsq.partial(fullmod,adj = T))
myrsq2<-rbind.data.frame(myrsq[1:4,],
                         c("TRUE","GLYCANS",sum(myrsq$partial.rsq[grep("Cluster",myrsq$variable)])))
myrsq2$partial.rsq<-as.numeric(myrsq2$partial.rsq)

ggplot2::ggplot(myrsq2,ggplot2::aes(x=partial.rsq,y=variable))+ggplot2::geom_col()+ggplot2::theme_classic()


############## Figure 5c #############

#AUC saturation analysis
daisy_prog<-list()
daisy_prog$igg<-daisy$igg[daisy$pheno$Group!="Non-progressor",]
daisy_prog$pheno<-daisy$pheno[daisy$pheno$Group!="Non-progressor",]

#random gly ridge model
auc_df<-data.frame("no_gly"=0,"mean"=0,"lowerCI"=0,"upperCI"=0)[-1,]

for (i in seq(1,200,10)){
  x <- replicate(1000, {glmnet_auc(daisy_prog,i)})
  auc_df<-rbind.data.frame(auc_df,
                           c(i,mean(x),mean(x)-1.96*sd(x),mean(x)+1.96*sd(x)))
}
colnames(auc_df)<-c("no_gly","AUC","lowerCI","upperCI")

#plot
ggplot2::ggplot(data=auc_df,ggplot2::aes(x=no_gly,y=AUC))+
  ggplot2::geom_point()+
  ggplot2::geom_errorbar(ggplot2::aes(ymin=lowerCI,ymax=upperCI))+
  ggplot2::theme_classic()#+
# geom_function(fun=function(x)log(x,3))



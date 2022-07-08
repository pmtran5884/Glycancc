#ACA predictive models
#Date: March 9, 2022
#Author: Paul Tran
#Updated: July 7, 2022

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

### rand mod 1
daisy$combined$Grouprand1<-sample(daisy$combined$Group)
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
              data = daisy$combined,
              subset = which(daisy$combined$Grouprand1!="Non-progressor"),
              na.action = na.omit,
              family = "binomial")

#2
daisy$combined$Grouprand2<-sample(daisy$combined$Group)
randmod2<-glm(Grouprand2 ~
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
              subset = which(daisy$combined$Grouprand2!="Non-progressor"),
              na.action = na.omit,
              family = "binomial")
#3
daisy$combined$Grouprand3<-sample(daisy$combined$Group)
randmod3<-glm(Grouprand3 ~
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
              subset = which(daisy$combined$Grouprand3!="Non-progressor"),
              na.action = na.omit,
              family = "binomial")

#4
daisy$combined$Grouprand4<-sample(daisy$combined$Group)
randmod4<-glm(Grouprand4 ~
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
              subset = which(daisy$combined$Grouprand4!="Non-progressor"),
              na.action = na.omit,
              family = "binomial")

#5
daisy$combined$Grouprand5<-sample(daisy$combined$Group)
randmod5<-glm(Grouprand5 ~
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
              subset = which(daisy$combined$Grouprand5!="Non-progressor"),
              na.action = na.omit,
              family = "binomial")

#ROC plot
test_roc<-plot(pROC::roc(sample(fastDummies::dummy_cols(daisy$combined$Grouprand1[daisy$combined$Grouprand1!="Non-progressor"])[,4]),
                         predict(randmod1),ci=T),print.auc=T,col="grey")
test_roc<-plot(pROC::roc(sample(fastDummies::dummy_cols(daisy$combined$Grouprand2[daisy$combined$Grouprand2!="Non-progressor"])[,4]),
                         predict(randmod2),ci=T),print.auc=T,add=T,col="grey")
test_roc<-plot(pROC::roc(sample(fastDummies::dummy_cols(daisy$combined$Grouprand3[daisy$combined$Grouprand3!="Non-progressor"])[,4]),
                         predict(randmod3),ci=T),print.auc=T,add=T,col="grey")
test_roc<-plot(pROC::roc(sample(fastDummies::dummy_cols(daisy$combined$Grouprand4[daisy$combined$Grouprand4!="Non-progressor"])[,4]),
                         predict(randmod4),ci=T),print.auc=T,add=T,col="grey")
test_roc<-plot(pROC::roc(sample(fastDummies::dummy_cols(daisy$combined$Grouprand5[daisy$combined$Grouprand5!="Non-progressor"])[,4]),
                         predict(randmod5),ci=T),print.auc=T,add=T,col="grey")



test_roc<-plot(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                         predict(nullmod),ci=T),print.auc=T,add=T,col="black")
test_roc<-plot(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                                predict(fullmod),ci=T),print.auc=T,add=T,col="red",print.auc.y=0.4)
test_roc<-plot(pROC::roc(fastDummies::dummy_cols(daisy$combined$Group[daisy$combined$Group!="Non-progressor"])[,4],
                                predict(cl8mod),ci=T),print.auc=T,add=T,col="blue",print.auc.y=0.3)

#LRT to compare models

lmtest::lrtest(fullmod, nullmod)
lmtest::lrtest(cl8mod, nullmod)
lmtest::lrtest(cl8mod, fullmod)

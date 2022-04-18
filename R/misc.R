
#' entropy/shannon diversity index function
#' 
#' This function calculates the entropy of shannon diversity index for a given vector
#' 
#' @param pop vector
#' 
#' @return entropy value
#' @export
entropy<-function(pop){
  n = pop
  N = sum(pop, na.rm = T)
  p = n/N
  H = -sum(p*log(p),na.rm = T)
  H
}


#' outliers to na
#' 
#' This function takes a data frame and converts all column-wise outliers to na
#' An outlier is a value less than the 1st quartile minus the interquartile range
#' or a value greater than the 3rd quartile plus the interquartile range
#' 
#' @param iggdata numeric matrix
#' 
#' @return numeric matrix with outliers converted to "NA"
#' @export
outliers_to_na<-function(iggdata){
  sum_igg<-apply(iggdata,2,quantile,na.rm=T)
  igg_iqr<-sum_igg[4,]-sum_igg[2,]
  igg_lower<-sum_igg[2]-igg_iqr
  igg_higher<-sum_igg[4]+igg_iqr
  
  for (i in 1:dim(iggdata)[2]){
    iggdata[which(iggdata[,i]<igg_lower[i]),i]<-NA
    iggdata[which(iggdata[,i]>igg_higher[i]),i]<-NA
  }
  iggdata
}


#' rader/spider chart per glycan class
#' 
#' This function plots the significance of association for each anti-carbohydrate
#' antibody in a class against islet autoimmunity and progression to type 1 diabetes
#' accounting for covariates sex, first-degree relative, HLA risk, and age. This uses the 
#' ??univariate_models function
#' 
#' @param daisy igg data
#' @param gly_class glycan class to assess. Default = "Aminoglycoside antibiotics"
#' 
#' @return spider plot
#' @export
univariate_radar_chart<-function(daisy=daisy,gly_class="Aminoglycoside antibiotics",label=T){
  
  aminoglycosides<-glycan_classes$Glycan.ID[glycan_classes$Glycan.Class==gly_class]
  aminoglycosides<-aminoglycosides[!is.na(aminoglycosides)]
  aminoglycosides<-intersect(aminoglycosides,colnames(daisy$igg))
  daisy$combined<-cbind.data.frame(apply(daisy$igg[,aminoglycosides],2,scale),daisy$pheno)
  daisy$combined$Group<-factor(daisy$combined$Group,levels = c("Control","Non-progressor","Progressor"))
  daisy_aminoglyc_models<-univariate_models(daisy$combined,
                                            c("Group", 
                                              "Sex", "FDR", 
                                              "HLA_risk","Draw_Age"))
  
  radar_df<-rbind.data.frame(1.5,0,
                             -log10(daisy_aminoglyc_models$small_table$prog_p),
                             -log10(0.05),
                             -log10(daisy_aminoglyc_models$small_table$nonprog_p)
                             
                             # daisy_aminoglyc_models$small_table$nonprog_est)
  )
  colnames(radar_df)<-rownames(daisy_aminoglyc_models$small_table)
  areas <- c(rgb(0, 0, 1, 0.15),rgb(0, 0, 0, 0.1),rgb(1, 0, 0, 0.15))
  if (label==F){
    fmsb::radarchart(radar_df,
               pfcol = areas,
               pty = c(16,32,16),
               pcol = c("blue","black","red"),
               axistype = 0,
               # caxislabels = c("-log10(p-value)",NA,NA,"0.05"),
               calcex = 1,
               vlabels = NA,
               axislabcol = "black")
  }
  else {
    fmsb::radarchart(radar_df,
               pfcol = areas,
               pty = c(16,32,16),
               pcol = c("blue","black","red"),
               axistype = 0,
               # caxislabels = c("-log10(p-value)",NA,NA,"0.05"),
               calcex = 1,
               axislabcol = "black")
  }
}


#' Calculate the AUC of binomial glmnet model from random number of glycans
#' 
#' @param daisy_prog igg data subseted two only include two factors for comparison. In this case, control subjects and progressors.
#' @param no_glycans Number of glycans to randomly sample
#' @return AUC value
#' @export
glmnet_auc<-function(daisy_prog=daisy_prog,no_glycans){
  subigg<-daisy_prog$igg[,sample(1:dim(daisy_prog$igg)[2],no_glycans)]
  x<-data.matrix(cbind(daisy_prog$pheno[,-1],subigg))
  cv_fit <- glmnet::cv.glmnet(y=as.vector(daisy_prog$pheno$Group),x=x , alpha = 0,family="binomial")
  mypred<-predict(cv_fit$glmnet.fit,newx=as.matrix(x),type = "response",s=cv_fit$lambda.min)
  
  suppressMessages(pROC::roc(fastDummies::dummy_cols(daisy_prog$pheno$Group)[,2],as.numeric(mypred))$auc)
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{glmnet_auc}
\alias{glmnet_auc}
\title{Calculate the AUC of binomial glmnet model from random number of glycans}
\usage{
glmnet_auc(daisy_prog = daisy_prog, no_glycans)
}
\arguments{
\item{daisy_prog}{igg data subseted two only include two factors for comparison. In this case, control subjects and progressors.}

\item{no_glycans}{Number of glycans to randomly sample}
}
\value{
AUC value
}
\description{
Calculate the AUC of binomial glmnet model from random number of glycans
}

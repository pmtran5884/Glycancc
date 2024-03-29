% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_clustering.R
\name{plot_cluster_graph}
\alias{plot_cluster_graph}
\title{run umap algorithm on igg data}
\usage{
plot_cluster_graph(knn_graph, myumap, layout = "umap", cluster_method = "prop")
}
\arguments{
\item{knn_graph}{igraph object}

\item{cluster_method}{graph cluster method. Default = "louvain"}
}
\value{
igraph cluster object
}
\description{
run umap algorithm on igg data
}

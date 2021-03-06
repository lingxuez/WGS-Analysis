% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_size.R
\name{cluster_size}
\alias{cluster_size}
\title{Cluster size}
\usage{
cluster_size(vec, clustering, flag_vec)
}
\arguments{
\item{vec}{vector of test statistics}

\item{clustering}{vector of clustering assignments}

\item{flag_vec}{boolean vector of same length as \code{vec}, where indices with
\code{TRUE} will be omitted from calculation}
}
\value{
vector
}
\description{
Cluster size
}

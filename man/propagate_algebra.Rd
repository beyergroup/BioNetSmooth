% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetSmooth.R
\name{propagate_algebra}
\alias{propagate_algebra}
\title{propagates node scores on a graph}
\usage{
propagate_algebra(W, z, alpha, iter)
}
\arguments{
\item{W}{igraph object respresenting the laplacian matrix (degree normalized adjacency matrix) of graph generated from protein interaction list. It is a symmetric matrix with equal number of rows and columns which is equal to the number of nodes in graph.}

\item{z}{vector of mapped initial scores (gene expression values of each condition) corresponding to nodes in graph "G".}

\item{alpha}{ranging between 0 and 1, the fraction of intial scores that has to be diffused to the adjacent nodes. This is a tunable parameter. If 0.5 then 50 percent of intial score of the nodes will be spread to its neighbours, if it is 1 then 100 percent of initial scores will be spread to the adjacent nodes (complete loss of initial information).}

\item{iter}{an integer indicating number of runs of propagation of initial scores on network. At iter "t" every node propagates the scores received in the t-1 iteration to its neighbors. The node scores will get converged over the progression of iterations. The converged smooth matrix is found by taking the norm of smooth mat at t subtracted with smooth mat at t-1, the norm is exprected to be a very small value (in the ranges of 10^-6).}
}
\value{
vector of propagated node scores at the specified alpha and iterations
}
\description{
This is the actual function for diffusing the scores on the network. The sharing of initial scores to adjacent nodes in network is constrained by two parameters alpha and iter.
}
\examples{

}


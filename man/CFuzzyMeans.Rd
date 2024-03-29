% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo-c-mean.R
\name{CFuzzyMeans}
\alias{CFuzzyMeans}
\title{calssical c-mean algorithm}
\usage{
CFuzzyMeans(data, k, m, maxiter = 500, tol = 0.01, standardize = TRUE)
}
\arguments{
\item{data}{A dataframe with only numerical variable}

\item{k}{An integer describing the number of cluster to find}

\item{m}{An integer for the fuzzyness degree}

\item{maxiter}{An integer for the maximum number of iteration}

\item{tol}{The tolerance criterion used in the evaluateMatrices function for convergence assessment}

\item{standardize}{A boolean to specify if the variable must be centered and reduce (default = True)}
}
\description{
calssical c-mean algorithm
}
\examples{
mat1 <- rbind(c(0.45,0.55),c(0.35,0.65),c(0.8,0.2),c(0.5,0.5),c(0.9,0.1),c(0.7,0.3))
mat2 <- rbind(c(0.45,0.55),c(0.4,0.60),c(0.8,0.2),c(0.45,0.55),c(0.95,0.05),c(0.7,0.3))
evaluateMatrices(mat1,mat2,0.01)
}

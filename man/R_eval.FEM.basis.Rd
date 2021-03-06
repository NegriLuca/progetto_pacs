% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdaPDE.smoothing_R.R
\name{R_eval.FEM.basis}
\alias{R_eval.FEM.basis}
\title{Evaluate Finite Element bases and their Derivatives at a set of locations}
\usage{
R_eval.FEM.basis(FEMbasis, locations, nderivs = matrix(0,1,2))
}
\arguments{
\item{FEMbasis}{An \code{FEMbasis} object representing the Finite Element basis; See \link{create.FEM.basis}.}

\item{locations}{A 2-columns matrix with the spatial locations where the bases should be evaluated.}

\item{nderivs}{A vector of lenght 2 specifying the order of the partial derivatives of the bases to be evaluated. The vectors' entries can
be 0,1 or 2, where 0 indicates that only the basis functions, and not their derivatives, should be evaluated.}
}
\value{
A matrix of basis function values. Each row indicates the location where the evaluation has been taken, the column indicates the 
basis function evaluated
}
\description{
Only executed when the function \code{smooth.FEM.basis} is run with the option \code{CPP_CODE} = \code{FALSE}. It evaluates the Finite Element basis functions and their derivatives up to order 2 at the specified set of locations. 
This version of the function is implemented using only R code. It is called by \link{R_smooth.FEM.basis}.
}
\seealso{
\code{\link{R_eval.FEM}}
}

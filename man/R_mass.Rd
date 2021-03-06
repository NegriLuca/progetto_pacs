% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdaPDE.smoothing_R.R
\name{R_mass}
\alias{R_mass}
\title{Compute the mass matrix}
\usage{
R_mass(FEMbasis)
}
\arguments{
\item{FEMbasis}{A \code{FEM} object representing the Finite Element basis. See \code{\link{create.FEM.basis}}.}
}
\value{
A square matrix with the integrals of all the basis' functions pairwise products.
The dimension of the matrix is equal to the number of the nodes of the mesh.
}
\description{
Only executed when \code{smooth.FEM.basis} is run with the option  \code{CPP_CODE} = \code{FALSE}. It computes the mass matrix. The element (i,j) of this matrix contains the integral over the domain of the product between the ith and kth element 
of the Finite Element basis. As common practise in Finite Element Analysis, this quantities are computed iterating over all the mesh triangles.
}
\seealso{
\code{\link{R_stiff}}
}

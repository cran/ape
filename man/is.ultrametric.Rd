\name{is.ultrametric}
\alias{is.ultrametric}
\title{Test if a Tree is Ultrametric}
\usage{
is.ultrametric(phy, tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{tol}{a numeric >= 0, variation below this value are considered
    non-significant (see details).}
}
\description{
  This function computes the distances from each tip to the root: if the
  variance of these distances is null, the tree is considered as
  ultrametric.
}
\value{
  a logical: \code{TRUE} if the tree is ultrametric, \code{FALSE}
  otherwise.
}
\details{
  The default value for \code{tol} is based on the numerical
  characteristics of the machine R is running on.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{is.binary.tree}}, \code{\link[base]{.Machine}}
}
\keyword{utilities}

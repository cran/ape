\name{skyline.plot}
\alias{skyline.plot}
\title{Skyline Plot of Estimated Effective Population Size}
\usage{
skyline.plot(phy, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{...}{further arguments to be passed to \code{plot()}.}
}
\description{
  This plots, on the current graphical device, the effective population
  sizes through time as estimated from a phylogenetic tree (see Pybus et
  al., 2000, for the rationale of the method).
}
\references{
  Pybus, O. G., Rambaut, A. and Harvey, P. H. (2000) An integrated
  framework for the inference of viral population history from
  reconstructed genealogies. \emph{Genetics}, \bold{155}, 1429--1437.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link[base]{plot}} for the basic plotting function in R,
  \code{\link{ltt.plot}}, \code{\link{branching.times}}
}
\keyword{hplot}

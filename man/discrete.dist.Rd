\name{discrete.dist}
\alias{discrete.dist}
\title{Taxonomic Level from a Pairwise Distances Matrix}
\usage{
  discrete.dist(dist, inf, sup)
}
\arguments{
  \item{dist}{The pairwise distances matrix}
  \item{inf}{The lower bound of the class to use}
  \item{sup}{The upper bound of the class to use}
}
\description{
  Takes all distances in \code{dist} and set TRUE if the distance is < \code{sup} and >= \code{inf}.
  Sends the resulting boolean matrix.
}
\value{
  a boolean matrix of same size than \code{dist}.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
	\code{\link{weight.taxo}} and \code{\link{cophenetic.phylo}} for other distances.
  This function is used in \code{\link{correlogram.phylo}}.
}
\keyword{regression}
\keyword{manip}

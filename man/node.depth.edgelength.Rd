\name{node.depth.edgelength}
\alias{node.depth.edgelength}
\title{Distance from Root to Nodes and Tips}
\description{
  This function returns the distance from the root to each node and tip
  with respect to the edge lengths.
}
\usage{
node.depth.edgelength(x, el)
}
\arguments{
  \item{x}{the matrix `edge' of an object of class "phylo".}
  \item{el}{the vector `edge.length' of an object of class "phylo".}
}
\value{
  A numeric vector with `names' being the ``numbers'' of the matrix
  `edge'.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\note{
  This function is actually called by plot.phylo(), and is unlikely to
  be called directly. However, it is documented separately since it may
  be useful for developpers.
}
\seealso{
  \code{\link{node.depth}}, \code{\link{plot.phylo}}
}
\keyword{manip}

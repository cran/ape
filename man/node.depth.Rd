\name{node.depth}
\alias{node.depth}
\title{Depth of Nodes and Tips}
\description{
  This function returns the depth of nodes and tips given by the number
  of descendants (1 is returned for tips).
}
\usage{
node.depth(x)
}
\arguments{
  \item{x}{the matrix `edge' of an object of class "phylo".}
}
\details{
  The depth of a node is computed as the number of tips which are its
  descendants. The value of 1 is given to the tips.
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
  \code{\link{plot.phylo}}
}
\keyword{manip}

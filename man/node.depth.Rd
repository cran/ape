\name{node.depth}
\alias{node.depth}
\title{Depth of Nodes and Tips}
\description{
  This function returns the depth of nodes and tips given by the number
  of descendants (1 is returned for tips).
}
\usage{
node.depth(phy)
}
\arguments{
  \item{phy}{an object of class "phylo".}
}
\details{
  The depth of a node is computed as the number of tips which are its
  descendants. The value of 1 is given to the tips.
}
\value{
  A numeric vector indexed with the node numbers of the matrix `edge' of
  \code{phy}.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{plot.phylo}}
}
\keyword{manip}

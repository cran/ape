\name{node.leafnumber}
\alias{node.leafnumber}
\title{Get the Number of Leaves for a Node}
\usage{
	node.leafnumber(phy, node)
}
\arguments{
	\item{phy}{An object of class \code{phylo} representing the tree to parse.}
	\item{node}{A string giving the node index.}
}
\description{
	Compute the number of leaves of a subtree defined by a particular node.
}
\value{
	The number of leaves of the subtree defined by \code{node}.
	If \code{node} is a leaf or node a valid node name, 0 is returned.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
  \code{\link{node.depth}}, \code{\link{node.sons}}
}
\keyword{manip}

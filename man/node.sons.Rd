\name{node.sons}
\alias{node.sons}
\title{Get the list of Sons for a Node}
\usage{
	node.sons(phy, node)
}
\arguments{
	\item{phy}{An object of class \code{phylo} representing the tree to parse.}
	\item{node}{A string giving the node index.}
}
\description{
	Get all tags of son nodes from a given parent node.
}
\value{
	A vector of strings with the names of all son nodes.
	If \code{node} is a leaf or node a valid node name, a void vector is returned.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
  \code{\link{node.depth}}, \code{\link{node.depth.edgelength}}, \code{\link{node.leafnumber}}.
}
\keyword{manip}


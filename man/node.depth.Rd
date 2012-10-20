\name{node.depth}
\alias{node.depth}
\alias{node.depth.edgelength}
\alias{node.height}
\alias{node.height.clado}
\title{Depth and Heights of Nodes and Tips}
\description{
  These functions return the depth or height of nodes and tips.
}
\usage{
node.depth(phy)
node.depth.edgelength(phy)
node.height(phy)
node.height.clado(phy)
}
\arguments{
  \item{phy}{an object of class "phylo".}
}
\details{
  \code{node.depth} computes the depth of a node as the number of tips
  which are its descendants. The value of 1 is given to the tips.

  \code{node.depth.edgelength} does the same but using branch lengths.

  \code{node.height} computes the heights of nodes and tips as plotted
  by a phylogram.

  \code{node.height.clado} does the same but for a cladogram.
}
\value{
  A numeric vector indexed with the node numbers of the matrix `edge' of
  \code{phy}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}
}
\keyword{manip}

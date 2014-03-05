\name{is.binary.tree}
\alias{is.binary.tree}
\title{Test for Binary Tree}
\description{
  This function tests whether a phylogenetic tree is binary.
}
\usage{
is.binary.tree(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
}
\details{
  The test differs slightly whether the tree is rooted or not. An
  urooted tree is considered binary if all its nodes are of degree three
  (i.e., three edges connect to each node). A rooted tree is considered
  binary if all nodes (including the root node) have exactly two
  descendant nodes, so that they are of degree three expect the root
  which is of degree 2.
}
\value{
  \code{is.binary.tree} returns \code{TRUE} if \code{tree} is a fully
  binary phylogenetic tree, otherwise it returns \code{FALSE}.
}
\seealso{
\code{\link{is.rooted}}, \code{\link{is.ultrametric}}, \code{\link{multi2di}}
}
\author{Korbinian Strimmer (\url{http://www.stat.uni-muenchen.de/~strimmer/})}
\examples{
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick)
is.binary.tree(tree.hiv)
plot(tree.hiv)
}
\keyword{logic}

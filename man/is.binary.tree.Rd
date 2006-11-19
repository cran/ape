\name{is.binary.tree}
\alias{is.binary.tree}
\title{Test for Binary Tree}
\usage{
is.binary.tree(phy)
}
\arguments{
  \item{phy}{phylogenetic tree (i.e. an object of class \code{"phylo"}).}
}
\description{
 This function tests whether a phylogenetic tree is binary, i.e. whether every node
 (including the root node) has exactly two descendant nodes.  In this case the number
 of edges in the tree equals 2 n - 2 where n is the number of taxa (tips) in the tree.
}
\value{
  \code{is.binary.tree} returns \code{TRUE} if \code{tree}
   is a fully binary phylogenetic tree, otherwise it returns \code{FALSE}.
}
\seealso{
\code{\link{read.tree}}, \code{\link{is.ultrametric}}
}

\author{Korbinian Strimmer (\url{http://www.stat.uni-muenchen.de/~strimmer/})}

\examples{
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick) # load tree

is.binary.tree(tree.hiv)

plot(tree.hiv)

unlink("hivtree.phy") # delete the file "hivtree.phy"
}
\keyword{logic}

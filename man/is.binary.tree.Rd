\name{is.binary.tree}
\alias{is.binary.tree}
\title{Test for Binary Tree}
\usage{
is.binary.tree(tree)
}
\arguments{
  \item{tree}{phylogenetic tree (i.e. an object of class \code{"phylo"}).}
}
\description{
 This function tests whether a phylogenetic tree is binary.
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
library(ape)
data("hivtree.newick") # example tree
cat(hivtree.newick, file = "hivtree.phy", sep = "\n") # generate file "hivtree.phy"

tree.hiv <- read.tree("hivtree.phy")
is.binary.tree(tree.hiv)

plot(tree.hiv)

unlink("hivtree.phy") # delete the file "hivtree.phy"
}
\keyword{logic}

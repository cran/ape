\name{collapse.singles}
\alias{collapse.singles}
\title{Collapse Single Nodes}
\description{
  This function deletes the single nodes (i.e., with a single
  descendant) in a tree.
}
\usage{
collapse.singles(tree, root.edge = FALSE)
}
\arguments{
  \item{tree}{an object of class \code{"phylo"}.}
  \item{root.edge}{whether to get the singleton edges from the root
    until the first bifurcating node and put them as \code{root.edge} of
    the returned tree. By default, this is ignored or if the tree has no
    edge lengths (see examples).}
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{read.tree}}
}
\examples{
## a tree with 3 tips and 3 nodes:
e <- c(4L, 6L, 6L, 5L, 5L, 6L, 1L, 5L, 3L, 2L)
dim(e) <- c(5, 2)
tr <- structure(list(edge = e, tip.label = LETTERS[1:3], Nnode = 3L),
                class = "phylo")
tr
## the following shows that node #4 (ie, the root) is a singleton
## and node #6 is the first bifurcating node
tr$edge
## Normally, a tree has less nodes than it has tips:
## > plot(tr)
## Error in plot.phylo(tr) :
##   there are single (non-splitting) nodes in your tree;
##   you may need to use collapse.singles()
collapse.singles(tr) # only 2 nodes
## give branch lengths to use the 'root.edge' option:
tr$edge.length <- runif(5)
str(collapse.singles(tr, TRUE)) # has a 'root.edge'
}
\keyword{manip}

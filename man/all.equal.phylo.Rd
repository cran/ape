\name{all.equal.phylo}
\alias{all.equal.phylo}
\title{Global Comparison of two Phylogenies}
\usage{
\method{all.equal}{phylo}(target, current, \dots)
}
\arguments{
  \item{target}{an object of class \code{"phylo"}.}
  \item{current}{an object of class \code{"phylo"}.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  This function makes a global comparison of two phylogenetic trees.
}
\details{
  This function is meant to be an adaptation of the generic function
  \code{all.equal} for the comparison of phylogenetic trees.

  A single phylogenetic tree may have several representations in the Newick
  format and in the \code{"phylo"} class of objects used in `ape'. One
  aim of the present function is to be able to identify whether two
  objects of class \code{"phylo"} represent the same phylogeny.
}
\value{
  If both trees are strictly similar, the logical \code{TRUE} is
  returned. Otherwise, the different elements of the trees (numbers of
  tips, numbers of nodes, tip labels, tree topologies) are compared, and
  the results of these comparisons are printed.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link[base]{all.equal}} for the generic R function
}
\examples{
### maybe the simplest example of two representations
### for the same rooted tree
cat("((a:1,b:1):1,c:2);", file = "t1.tre")
cat("(c:2,(a:1,b:1):1);", file = "t2.tre")
t1 <- read.tree("t1.tre")
t2 <- read.tree("t2.tre")
all.equal(t1, t2)
### clean-up
unlink("t1.tre")
unlink("t2.tre")
}
\keyword{manip}

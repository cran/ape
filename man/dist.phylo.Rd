\name{dist.phylo}
\alias{dist.phylo}
\title{Pairwise Distances from a Phylogenetic Tree}
\usage{
dist.phylo(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
}
\description{
  This function computes the pairwise distances between the pairs of
  tips from a phylogenetic tree using the branch lengths of the latter.
}
\value{
  a numeric matrix with the names of the tips (as given by the element
  \code{tip.label} of the argument \code{phy}) as colnames and rownames.
}
\author{Julien Claude \email{claude@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}} to read tree files in Newick format
}
\keyword{manip}

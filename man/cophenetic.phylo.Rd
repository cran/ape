\name{cophenetic.phylo}
\alias{cophenetic.phylo}
\alias{dist.phylo}
\title{Pairwise Distances from a Phylogenetic Tree}
\usage{
\method{cophenetic}{phylo}(x, full = FALSE)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{full}{a logical indicating whether to return the distances among
    all tips and nodes (if \code{TRUE}); the default is to return only
    the distances among tips.}
}
\description{
  This function computes the pairwise distances between the pairs of
  tips from a phylogenetic tree using its branch lengths.
}
\value{
  a numeric matrix with colnames and rownames set, by default, to the
  names of the tips (as given by the element \code{tip.label} of the
  argument \code{phy}), or, if \code{full = TRUE}, the numbers of the
  tips and the nodes (as given by the element \code{edge}).
}
\author{
  Julien Claude \email{claude@isem.univ-montp2.fr} and
  Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}
}
\seealso{
  \code{\link{read.tree}} to read tree files in Newick format,
  \code{\link[stats]{cophenetic}} for the generic function
}
\keyword{manip}

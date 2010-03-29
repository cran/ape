\name{reorder.phylo}
\alias{reorder.phylo}
\title{Internal Reordering of Trees}
\description{
  This function changes the internal structure of a phylogeny stored as
  an object of class \code{"phylo"}. The tree returned is the same than
  the one input, but the ordering of the edges could be different.
}
\usage{
\method{reorder}{phylo}(x, order = "cladewise", ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{order}{a character string: either \code{"cladewise"} (the
    default), or \code{"pruningwise"}, or any unambiguous abbreviation
    of these.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  Because in a tree coded as an object of class \code{"phylo"} each
  branch is represented by a row in the element `edge', there is an
  arbitrary choice for the ordering of these rows. \code{reorder} allows
  to reorder these rows according to two rules: in the
  \code{"cladewise"} order each clade is formed by a series of
  contiguous rows; this is the order returned by
  \code{\link{read.tree}}. In the \code{"pruningwise"} order, rows are
  arranged so that ``pruning'' the tree (or post-order tree traversal)
  can be done by descending along the rows of `edge'. The possible
  multichotomies and branch lengths are preserved.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}} to read tree files in Newick format,
  \code{\link[stats]{reorder}} for the generic function
}
\examples{
data(bird.families)
tr <- reorder(bird.families, "p")
all.equal(bird.families, tr) # uses all.equal.phylo actually
all.equal.list(bird.families, tr) # bypasses the generic
}

\keyword{manip}

\name{reorder.phylo}
\alias{reorder.phylo}
\title{Internal Reordering of Trees}
\description{
  This function changes the internal structure of a phylogeny stored as
  an object of class \code{"phylo"}. The tree returned is the same than
  the one input, but the ordering of the edges could be different.
}
\usage{
\method{reorder}{phylo}(x, order = "cladewise", index.only = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{order}{a character string: either \code{"cladewise"} (the
    default), \code{"postorder"}, \code{"pruningwise"}, or any
    unambiguous abbreviation of these.}
  \item{index.only}{should the function return only the ordered indices
    of the rows of the edge matrix?}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  Because in a tree coded as an object of class \code{"phylo"} each
  branch is represented by a row in the element `edge', there is an
  arbitrary choice for the ordering of these rows. \code{reorder} allows
  to reorder these rows according to three rules: in the
  \code{"cladewise"} order each clade is formed by a series of
  contiguous rows. In the \code{"postorder"} order, the rows are
  arranged so that computations following pruning-like algorithm the
  tree (or postorder tree traversal) can be done by descending along
  these rows (conversely, a preorder tree traversal can be performed by
  moving from the last to the first row). The \code{"pruningwise"} order
  is an alternative ``pruning'' order which is actually a bottom-up
  traversal order (Valiente 2002). (This third choice might be removed
  in the future as it merely duplicates the second one which is more
  efficient.) The possible multichotomies and branch lengths are preserved.

  Note that for a given order, there are several possible orderings of
  the rows of `edge'.
}
\value{
  an object of class \code{"phylo"} (with the attribute \code{"order"}
  set accordingly), or a numeric vector if \code{index.only = TRUE}.
}
\references{
  Valiente, G. (2002) \emph{Algorithms on Trees and Graphs.} New York:
  Springer.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}} to read tree files in Newick format,
  \code{\link[stats]{reorder}} for the generic function
}
\examples{
data(bird.families)
tr <- reorder(bird.families, "postorder")
all.equal(bird.families, tr) # uses all.equal.phylo actually
all.equal.list(bird.families, tr) # bypasses the generic
}
\keyword{manip}

\name{c.phylo}
\alias{c.phylo}
\alias{c.multiPhylo}
\alias{.compressTipLabel}
\alias{.uncompressTipLabel}
\title{Building Lists of Trees}
\description{
  These functions help to build lists of trees of class \code{"multiPhylo"}.
}
\usage{
\method{c}{phylo}(..., recursive = FALSE)
\method{c}{multiPhylo}(..., recursive = FALSE)
.compressTipLabel(x)
.uncompressTipLabel(x)
}
\arguments{
  \item{\dots}{one or several objects of class \code{"phylo"} or
    \code{"multiPhylo"}.}
  \item{recursive}{for compatibily with the generic (do not change).}
  \item{x}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
}
\details{
  These \code{c} methods do not check all the arguments, so it is the
  user's responsibility to make sure that only objects of the same class
  (either \code{"phylo"} or \code{"multiPhylo"}) are used.

  \code{.compressTipLabel} transforms an object of class
  \code{"multiPhylo"} by checking that all trees have the same tip
  labels and renumbering the tips in the \code{edge} matrix so that the
  tip numbers are also the same taking the first tree as the reference
  (duplicated labels are not allowed). The returned object has a unique
  vector of tip labels (\code{attr(x, "TipLabel")}).

  \code{.uncompressTipLabel} does the reverse operation.
}
\value{
  An object of class \code{"multiPhylo"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{summary.phylo}}, \code{\link{multiphylo}}
}
\examples{
x <- c(rtree(4), rtree(2))
x
y <- c(rtree(4), rtree(4))
z <- c(x, y)
z
print(z, TRUE)
try(.compressTipLabel(x)) # error
a <- .compressTipLabel(y)
.uncompressTipLabel(a) # back to y
## eventually compare str(a) and str(y)
}
\keyword{manip}

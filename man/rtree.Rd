\name{rtree}
\alias{rtree}
\alias{rcoal}
\alias{rmtree}
\title{Generates Random Trees}
\description{
  These functions generate trees by splitting randomly the edges
  (\code{rtree}) or randomly clustering the tips (\code{rcoal}).
  \code{rtree} generates general trees, and \code{rcoal} generates
  coalescent trees. The algorithms are described in Paradis (2012).
}
\usage{
rtree(n, rooted = TRUE, tip.label = NULL, br = runif, ...)
rcoal(n, tip.label = NULL, br = "coalescent", ...)
rmtree(N, n, rooted = TRUE, tip.label = NULL, br = runif, ...)
}
\arguments{
  \item{n}{an integer giving the number of tips in the tree.}
  \item{rooted}{a logical indicating whether the tree should be rooted
    (the default).}
  \item{tip.label}{a character vector giving the tip labels; if not
    specified, the tips "t1", "t2", ..., are given.}
  \item{br}{one of the following: (i) an \R function used to generate the
    branch lengths (\code{rtree}; use \code{NULL} to simulate only a
    topology), or the coalescence times (\code{rcoal}); (ii) a character
    to simulate a genuine coalescent tree for \code{rcoal} (the
    default); or (iii) a numeric vector for the branch lengths or the
    coalescence times.}
  \item{\dots}{further argument(s) to be passed to \code{br}.}
  \item{N}{an integer giving the number of trees to generate.}
}
\details{
  The trees generated are bifurcating. If \code{rooted = FALSE} in
  (\code{rtree}), the tree is trifurcating at its root.

  The default function to generate branch lengths in \code{rtree} is
  \code{runif}. If further arguments are passed to \code{br}, they need
  to be tagged (e.g., \code{min = 0, max = 10}).

  \code{rmtree} calls successively \code{rtree} and set the class of
  the returned object appropriately.
}
\value{
  An object of class \code{"phylo"} or of class \code{"multiPhylo"} in
  the case of \code{rmtree}.
}
\references{
  Paradis, E. (2012) \emph{Analysis of Phylogenetics and Evolution with
    R (Second Edition).} New York: Springer.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{stree}}, \code{\link{rlineage}}
}
\examples{
layout(matrix(1:9, 3, 3))
### Nine random trees:
for (i in 1:9) plot(rtree(20))
### Nine random cladograms:
for (i in 1:9) plot(rtree(20, FALSE), type = "c")
### generate 4 random trees of bird orders:
data(bird.orders)
layout(matrix(1:4, 2, 2))
for (i in 1:4)
  plot(rcoal(23, tip.label = bird.orders$tip.label), no.margin = TRUE)
layout(1)
par(mar = c(5, 4, 4, 2))
}
\keyword{datagen}

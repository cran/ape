\name{rtree}
\alias{rtree}
\alias{rcoal}
\title{Generates Random Trees}
\usage{
rtree(n, rooted = TRUE, tip.label = NULL, br = runif, ...)
rcoal(n, tip.label = NULL, br = rexp, ...)
}
\arguments{
  \item{n}{an integer giving the number of tips in the tree.}
  \item{rooted}{a logical indicating whether the tree should be rooted
    (the default).}
  \item{tip.label}{a character vector giving the tip labels; if not
    specified, the tips "t1", "t2", ..., are given.}
  \item{br}{either an R function used to generate the branch lengths
    (\code{rtree}) or the coalescence times (\code{rcoal}), or
    \code{NULL} to gave no branch lengths in the tree.}
  \item{...}{further argument(s) to be passed to \code{br}.}
}
\description{
  These functions generate trees by splitting randomly the edges
  (\code{rtree}) or randomly clustering the tips (\code{rcoal}).
  \code{rtree} generates general (non-ultrametric) trees, and
  \code{rcoal} generates coalescent (ultrametric) trees.
}
\details{
  The trees generated are bifurcating. If \code{rooted = FALSE} in
  (\code{rtree}), the tree is trifurcating at its `root'.

  The default function to generate branch lengths in \code{rtree} is
  \code{runif}. In \code{rcoal} \code{rexp} is used to generate the
  inter-node distances. If further arguments are passed to \code{br},
  they need to be tagged (e.g., \code{min = 0, max = 10}).
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
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
layout(matrix(1))
}
\keyword{datagen}

\name{rotate}
\alias{rotate}
\title{Rotate an Internal Branch of a Tree}
\usage{
rotate(phy, group)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{group}{a vector of mode numeric or character specifying the tips
    defining the clade to rotate, or the string \code{"all"} meaning to
    rotate the whole tree.}
}
\description{
  This function rotates a clade within a tree.
}
\details{
  The argument \code{outgroup} can be one of the three followings: (1) a
  vector of mode character (= the labels of the tips); (2) a vector of
  mode numeric (= the numbers of these labels in the vector
  \code{phy$tip.label}), or (3) the string \code{"all"}. In cases (1)
  and (2), the clade is rotated around the most recent common ancestor
  of \code{group}. If the latter is not monophyletic, the operation
  fails and an error message is issued.

  Note that the tree returned is the same than the original one, just
  its representation is changed.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{bind.tree}}, \code{\link{drop.tip}}, \code{\link{root}}
}
\examples{
data(bird.orders)
layout(matrix(1:4, 2, 2))
op <- par(cex = 0.5)
plot(bird.orders)
plot(rotate(bird.orders, 1:5))
plot(rotate(bird.orders, 6:23))
plot(rotate(bird.orders, "all"))
layout(matrix(1))
par(op)
### The tree hasn't changed:
all.equal(bird.orders, rotate(bird.orders, "all"))
}
\keyword{manip}

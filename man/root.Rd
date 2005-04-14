\name{root}
\alias{root}
\title{Roots Phylogenetic Trees}
\usage{
root(phy, outgroup)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{outgroup}{a vector of mode numeric or character specifying the
    new outgroup.}
}
\description{
  This function reroots a phylogenetic tree with respect to the
  specified outgroup.
}
\details{
  The argument \code{outgroup} can be either character or numeric. In
  the first case, it gives the labels of the tips of the new outgroup;
  in the second case the numbers of these labels in the vector
  \code{phy$tip.label} are given.

  If \code{outgroup} is of length one (i.e. a single value), then the
  tree is rerooted using the node below this tip as the new root.

  If \code{outgroup} is of length two or more, the most recent common
  ancestor (MRCA) is used as the new root. Note that the tree is really
  unrooted before being rerooted, so that if \code{outgroup} is already
  the outgroup, then the returned tree is not the same than the original
  one (see examples). If \code{outgroup} is not monophyletic, the
  operation fails and an error message is issued.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{bind.tree}}, \code{\link{drop.tip}}, \code{\link{rotate}}
}
\examples{
data(bird.orders)
plot(root(bird.orders, 1))
plot(root(bird.orders, 1:5))
}
\keyword{manip}

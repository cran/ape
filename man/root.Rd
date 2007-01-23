\name{root}
\alias{root}
\alias{unroot}
\alias{is.rooted}
\title{Roots Phylogenetic Trees}
\usage{
root(phy, outgroup)
unroot(phy)
is.rooted(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{outgroup}{a vector of mode numeric or character specifying the
    new outgroup.}
}
\description{
  \code{root} reroots a phylogenetic tree with respect to the specified
  outgroup.

  \code{unroot} unroots a phylogenetic tree, or returns it unchanged if
  it is already unrooted.

  \code{is.rooted} tests whether a tree is rooted.
}
\details{
  The argument \code{outgroup} can be either character or numeric. In
  the first case, it gives the labels of the tips of the new outgroup;
  in the second case the numbers of these labels in the vector
  \code{phy$tip.label} are given.

  If \code{outgroup} is of length one (i.e., a single value), then the
  tree is rerooted using the node below this tip as the new root.

  If \code{outgroup} is of length two or more, the most recent common
  ancestor (MRCA) is used as the new root. Note that the tree is really
  unrooted before being rerooted, so that if \code{outgroup} is already
  the outgroup, then the returned tree is not the same than the original
  one (see examples). If \code{outgroup} is not monophyletic, the
  operation fails and an error message is issued.

  A tree is considered rooted if either only two branches connect to the
  root, or if there is a \code{root.edge} element. In all other cases,
  \code{is.rooted} returns \code{FALSE}.
}
\value{
  an object of class \code{"phylo"} for \code{root} and \code{unroot}; a
  single logical value for \code{is.rooted}.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{bind.tree}}, \code{\link{drop.tip}}
}
\examples{
data(bird.orders)
plot(root(bird.orders, 1))
plot(root(bird.orders, 1:5))

tr <- root(bird.orders, 1)
is.rooted(bird.orders) # yes!
is.rooted(tr)          # no!
### This is because the tree has been unrooted first before rerooting.
### You can delete the outgroup...
is.rooted(drop.tip(tr, "Struthioniformes"))
### ... or resolve the basal trichotomy:
is.rooted(multi2di(tr))
### To keep the basal trichotomy but forcing the tree as rooted:
tr$root.edge <- 0
is.rooted(tr)
}
\keyword{manip}

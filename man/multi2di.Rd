\name{multi2di}
\alias{multi2di}
\alias{di2multi}
\title{Collapse and Resolve Multichotomies}
\description{
  These two functions collapse or resolve multichotomies in phylogenetic
  trees (objects of class \code{"phylo"}).
}
\usage{
multi2di(phy)
di2multi(phy, tol = 1e-08)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{tol}{a numeric value giving the tolerance to consider a branch
    length significantly greater than zero.}
}
\details{
  \code{multi2di} transforms all multichotomies into a series of
  dichotomies with one (or several) branch(es) of length zero.

  \code{di2multi} deletes all branches smaller than \code{tol} and
  collapses the corresponding dichotomies into a multichotomy.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\value{
  Both functions return an object of class \code{"phylo"}.
}
\examples{
data(bird.families)
is.binary.tree(bird.families)
is.binary.tree(multi2di(bird.families))
all.equal(di2multi(multi2di(bird.families)), bird.families)
}
\keyword{manip}

\name{multi2di}
\alias{multi2di}
\alias{di2multi}
\title{Collapse and Resolve Multichotomies}
\description{
  These two functions collapse or resolve multichotomies in phylogenetic
  trees.
}
\usage{
multi2di(phy, random = TRUE)
di2multi(phy, tol = 1e-08)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{random}{a logical value specifying whether to resolve the
    multichotomies randomly (the default) or in the order they appear in
    the tree (if \code{random = FALSE}).}
  \item{tol}{a numeric value giving the tolerance to consider a branch
    length significantly greater than zero.}
}
\details{
  \code{multi2di} transforms all multichotomies into a series of
  dichotomies with one (or several) branch(es) of length zero.

  \code{di2multi} deletes all branches smaller than \code{tol} and
  collapses the corresponding dichotomies into a multichotomy.
}
\seealso{
\code{\link{is.binary.tree}}
}
\author{Emmanuel Paradis}
\value{
  Both functions return an object of class \code{"phylo"}.
}
\examples{
data(bird.families)
is.binary.tree(bird.families)
is.binary.tree(multi2di(bird.families))
all.equal(di2multi(multi2di(bird.families)), bird.families)
### To see the results of randomly resolving a trichotomy:
tr <- read.tree(text = "(a:1,b:1,c:1);")
layout(matrix(1:4, 2, 2))
for (i in 1:4)
  plot(multi2di(tr), use.edge.length = FALSE, cex = 1.5)
layout(matrix(1))
}
\keyword{manip}

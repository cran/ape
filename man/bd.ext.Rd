\name{bd.ext}
\alias{bd.ext}
\title{Extended Version of the Birth-Death Models to Estimate Speciation
  and Extinction Rates}
\usage{
bd.ext(phy, S)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{S}{a numeric vector giving the number of species for each tip.}
}
\description{
  This function fits by maximum likelihood a birth-death model to the
  combined phylogenetic and taxonomic data of a given clade. The
  phylogenetic data are given by a tree, and the taxonomic data by the
  number of species for the its tips.
}
\details{
  A re-parametrization of the birth-death model studied by Kendall
  (1948) so that the likelihood has to be maximized over \emph{d/b} and
  \emph{b - d}, where \emph{b} is the birth rate, and \emph{d} the death
  rate.

  The standard-errors of the estimated parameters are computed using a
  normal approximation of the maximum likelihood estimates.

  If the argument \code{S} has names, then they are matched to the tip
  labels of \code{phy}. The user must be careful here since the function
  requires that both series of names perfectly match, so this operation
  may fail if there is a typing or syntax error. If both series of names
  do not match, the values \code{S} are taken to be in the same order
  than the tip labels of \code{phy}, and a warning message is issued.

  Note that the function does not check that the tree is effectively
  ultrametric, so if it is not, the returned result may not be meaningful.
}
\references{
  Paradis, E. (2003) Analysis of diversification: combining phylogenetic
  and taxonomic data. \emph{Proceedings of the Royal Society of
    London. Series B. Biological Sciences}, \bold{270}, 2499--2505.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{birthdeath}}, \code{\link{branching.times}},
  \code{\link{diversi.gof}}, \code{\link{diversi.time}},
  \code{\link{ltt.plot}}, \code{\link{yule}}, \code{\link{yule.cov}}
}
\examples{
### An example from Paradis (2003) using the avian orders:
data(bird.orders)
### Number of species in each order from Sibley and Monroe (1990):
S <- c(10, 47, 69, 214, 161, 17, 355, 51, 56, 10, 39, 152,
       6, 143, 358, 103, 319, 23, 291, 313, 196, 1027, 5712)
bd.ext(bird.orders, S)
}
\keyword{models}

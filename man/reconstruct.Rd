\name{reconstruct}
\alias{reconstruct}
\title{Continuous Ancestral Character Estimation}
\description{
  This function estimates ancestral character states, and the associated
  uncertainty, for continuous characters. It mainly works as the
  \code{\link{ace}} function, from which it differs in the fact that
  optimisations are not performed by numerical algorithms but through
  matrix computations.
}
\usage{
reconstruct(x, phyInit, method = "ML", CI = TRUE)
}
\arguments{
  \item{x}{a vector or a factor.}
  \item{phyInit}{an object of class \code{"phylo"}.}
  \item{method}{a character specifying the method used for
    estimation. Three choices are possible: \code{"ML"}, \code{"REML"}
    or \code{"GLS"}.}
  \item{CI}{a logical specifying whether to return the 95\% confidence
    intervals of the ancestral state estimates.}
}
\details{
  The default model is Brownian motion where characters evolve randomly
  following a random walk. This model can be fitted by maximum
  likelihood (Felsenstein 1973, Schluter et al. 1997 - the default),
  residual maximum likelihood, or generalized least squares
  (\code{method = "GLS"}, Martins and Hansen 1997, Cunningham et
  al. 1998).
}
\value{
  a list with the following elements:

  \item{ace}{the estimates of the ancestral character values.}
  \item{CI95}{the estimated 95\% confidence intervals.}
  \item{sigma2}{if \code{method = "ML"}, the maximum likelihood estimate
    of the Brownian parameter.}
  \item{loglik}{if \code{method = "ML"}, the maximum log-likelihood.}
}
\references{
  Cunningham, C. W., Omland, K. E. and Oakley, T. H. (1998)
  Reconstructing ancestral character states: a critical
  reappraisal. \emph{Trends in Ecology & Evolution}, \bold{13},
  361--366.

  Felsenstein, J. (1973) Maximum likelihood estimation of evolutionary
  trees from continuous characters. \emph{American Journal of Human
  Genetics}, \bold{25}, 471--492.

  Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the
  comparative method: a general approach to incorporating phylogenetic
  information into the analysis of interspecific data. \emph{American
  Naturalist}, \bold{149}, 646--667.

  Schluter, D., Price, T., Mooers, A. O. and Ludwig, D. (1997)
  Likelihood of ancestor states in adaptive radiation. \emph{Evolution},
  \bold{51}, 1699--1711.

  Yang, Z. (2006) \emph{Computational Molecular Evolution}. Oxford:
  Oxford University Press.

  Royer-Carenzi, M. and Didier, G. (2014) Comparison of ancestral state
  reconstruction methods for continuous characters under directional
  evolution. Submitted.
}
\author{Manuela Royer-Carenzi, Gilles Didier}
\seealso{
  \code{\link{MPR}}, \code{\link{corBrownian}}, \code{\link{compar.ou}}

  Reconstruction of ancestral sequences can be done with the package
  \pkg{phangorn} (see function \code{?ancestral.pml}).
}
\examples{
### Some random data...
data(bird.orders)
x <- rnorm(23)
### Compare the three methods for continuous characters:
reconstruct(x, bird.orders)
reconstruct(x, bird.orders, method = "REML")
}
\keyword{models}

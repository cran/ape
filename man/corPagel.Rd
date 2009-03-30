\name{corPagel}
\alias{corPagel}
\alias{coef.corPagel}
\alias{corMatrix.corPagel}
\title{Pagel's ``lambda'' Correlation Structure}
\usage{
corPagel(value, phy, form = ~1, fixed = FALSE)
\method{corMatrix}{corPagel}(object, covariate = getCovariate(object),
                   corr = TRUE, ...)
\method{coef}{corPagel}(object, unconstrained = TRUE, \dots)
}
\arguments{
  \item{value}{the (initial) value of the parameter
    \eqn{\lambda}{lambda}.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{form}{(ignored).}
  \item{fixed}{a logical specifying whether \code{gls} should
    estimate \eqn{\lambda}{lambda} (the default) or keep it fixed.}
  \item{object}{an (initialized) object of class \code{"corPagel"}.}
  \item{covariate}{(ignored).}
  \item{corr}{a logical value specifying whether to return the
    correlation matrix (the default) or the variance-covariance matrix.}
  \item{unconstrained}{a logical value. If \code{TRUE} (the default),
    the coefficients are returned in unconstrained form (the same used
    in the optimization algorithm). If \code{FALSE} the coefficients are
    returned in ``natural'', possibly constrained, form.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  The correlation structure from the present model is derived from the
  Brownian motion model by multiplying the off-diagonal elements (i.e.,
  the covariances) by \eqn{\lambda}{lambda}. The variances are thus the
  same than for a Brownian motion model.
}
\value{
  an object of class \code{"corPagel"}, the coefficients from an object
  of this class, or the correlation matrix of an initialized object of
  this class. In most situations, only \code{corPagel} will be called
  by the user.
}
\author{Emmanuel Paradis}
\references{
  Freckleton, R. P., Harvey, P. H. and M. Pagel, M. (2002) Phylogenetic
  analysis and comparative data: a test and review of evidence.
  \emph{American Naturalist}, \bold{160}, 712--726.

  Pagel, M. (1999) Inferring the historical patterns of biological
  evolution. \emph{Nature}, \bold{401},877--884.
}
\keyword{models}

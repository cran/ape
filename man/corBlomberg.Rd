\name{corBlomberg}
\alias{corBlomberg}
\alias{coef.corBlomberg}
\alias{corMatrix.corBlomberg}
\title{Blomberg et al.'s Correlation Structure}
\usage{
corBlomberg(value, phy, form = ~1, fixed = FALSE)
\method{corMatrix}{corBlomberg}(object, covariate = getCovariate(object),
                   corr = TRUE, ...)
\method{coef}{corBlomberg}(object, unconstrained = TRUE, \dots)
}
\arguments{
  \item{value}{the (initial) value of the parameter \eqn{g}{g}.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{form}{(ignored).}
  \item{fixed}{a logical specifying whether \code{gls} should
    estimate \eqn{\gamma}{gamma} (the default) or keep it fixed.}
  \item{object}{an (initialized) object of class \code{"corBlomberg"}.}
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
  The ``ACDC'' (accelerated/decelerated) model assumes that continuous
  traits evolve under a Brownian motion model which rates accelerates
  (if \eqn{g}{g} < 1) or decelerates (if \eqn{g}{g} > 1) through
  time. If \eqn{g}{g} = 1, then the model reduces to a Brownian motion
  model.
}
\value{
  an object of class \code{"corBlomberg"}, the coefficients from an
  object of this class, or the correlation matrix of an initialized
  object of this class. In most situations, only \code{corBlomberg} will
  be called by the user.
}
\author{Emmanuel Paradis}
\references{
  Blomberg, S. P., Garland, Jr, T., and Ives, A. R. (2003) Testing for
  phylogenetic signal in comparative data: behavioral traits are more
  labile. \emph{Evolution}, \bold{57}, 717--745.
}
\keyword{models}

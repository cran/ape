\name{as.bitsplits}
\alias{as.bitsplits}
\alias{as.bitsplits.prop.part}
\alias{print.bitsplits}
\title{Conversion Among Split Classes}
\description{
  \code{as.bitsplits} is a generic function with a method for objects of
  class \code{"prop.part"}.
}
\usage{
as.bitsplits(x)
\method{as.bitsplits}{prop.part}(x)
\method{print}{bitsplits}(x, ...)
}
\arguments{
  \item{x}{an object of the appropriate class.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  an object of class \code{"bitsplits"} or NULL for \code{print}.
}
\author{Emmanuel Paradis}
\examples{
tr <- rtree(20)
pp <- prop.part(tr)
as.bitsplits(pp)
}
\keyword{manip}

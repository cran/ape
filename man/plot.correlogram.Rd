\name{plot.correlogram}
\alias{plot.correlogram}
\title{Plot a Correlogram}
\usage{
  plot.correlogram(x, test.level=0.05, ...)
}
\arguments{
  \item{x}{an object of class 'correlogram'}
  \item{test.level}{The level above which points will be considered as non-significant}
  \item{...}{other parameters to be passed to the \code{xyplot} function}
}
\description{
  Draw a correlogram previously computed.
}
\details{
  The methods rely on the \code{\link[lattice]{xyplot}} function of the
  \pkg{lattice} package. Each significant point is plotted by a filled
  circle. Level of significance may be passed as an argument to the
  function.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
  \code{\link[lattice]{xyplot}}, \code{\link{Correlogram}} for examples.
}
\keyword{hplot}

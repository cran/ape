\name{panel.superpose.correlogram}
\alias{panel.superpose.correlogram}
\title{Lattice panel used in 'plot.correlogramList'.}
\description{
  This a slightly modified version of the \code{\link[lattice]{panel.superpose}} class
  to allow the use of different characters in a group line.
  Only two characters are used: filled circle for significant autocorrelation,
  and open circle for non significant ones.
  The level of significance and all p-values are passed as arguments.
}
\usage{
	panel.superpose.correlogram(x, y=NUL, subscripts, groups,
             panel.groups = "panel.xyplot",
             col,
             col.line = superpose.line$col,
             col.symbol = superpose.symbol$col,
             pch = superpose.symbol$pch,
             p-values = NULL,
             test.level = 0.05,
             cex = superpose.symbol$cex,
             font = superpose.symbol$font,
             fontface = superpose.symbol$fontface,
             fontfamily = superpose.symbol$fontfamily,
             lty = superpose.line$lty,
             lwd = superpose.line$lwd,
             ...)
}
\arguments{
  \item{p-values}{ The p-values for each point. If NULL, then act as the \code{panel.superpose} class and use the \code{pch} argument.}
  \item{test.level}{ The level of significance to consider.}
  Other parameters have the same meaning as in \code{\link[lattice]{panel.superpose}}.
}
\details{
  This class is used in \code{\link{plot.correlogramList}}.
}
\author{Julien Dutheil \email{julien.dutheil@uni-montp2.fr}}
\seealso{
  \code{\link[lattice]{panel.superpose}}, \code{\link{plot.correlogramList}}.
}
\keyword{dplot}

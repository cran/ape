\name{panel.superpose.correlogram}
\alias{panel.superpose.correlogram}
\title{Lattice panel used in 'plot.correlogramList'.}
\description{
  This a slightly modified version of the
  \code{\link[lattice]{panel.superpose}} class to allow the use of
  different characters in a group line. Only two characters are used:
  filled circle for significant autocorrelation, and open circle for non
  significant ones. The level of significance and all p-values are
  passed as arguments.
}
\usage{
panel.superpose.correlogram(x, y = NULL, subscripts, groups,
     panel.groups = "panel.xyplot",
     col, col.line = superpose.line$col,
     col.symbol = superpose.symbol$col,
     pch = superpose.symbol$pch,
     p.values = NULL, test.level = 0.05,
     cex = superpose.symbol$cex,
     font = superpose.symbol$font,
     fontface = superpose.symbol$fontface,
     fontfamily = superpose.symbol$fontfamily,
     lty = superpose.line$lty,
     lwd = superpose.line$lwd, ...)
}
\arguments{
  \item{x}{coordinates of the points to be displayed.}
  \item{y}{(optional) another series of coordinates of the points to be
    displayed.}
  \item{subscripts}{subscripts giving indices in original data frame.}
  \item{groups}{a grouping variable.}
  \item{panel.groups}{panel function to be used for each group.}
  \item{col}{colors of the plotting symbols or the lines.}
  \item{col.line}{colors of the lines.}
  \item{col.symbol}{colors of the plotting symbols.}
  \item{pch}{the type of plotting symbols.}
  \item{p.values}{the p-values for each point. If NULL, then act as the
    \code{panel.superpose} class and use the \code{pch} argument.}
  \item{test.level}{the level of significance to consider.}
  \item{cex}{character expansion coefficient.}
  \item{font}{font for the text.}
  \item{fontface}{font face.}
  \item{fontfamily}{font family.}
  \item{lty}{type of lines.}
  \item{lwd}{line width.}
  \item{...}{other arguments passed to
    \code{\link[lattice]{panel.superpose}}.}
}
\details{
  This class is used in \code{\link{plot.correlogramList}}.
}
\author{Julien Dutheil \email{julien.dutheil@uni-montp2.fr}}
\seealso{
  \code{\link[lattice]{panel.superpose}},
  \code{\link{plot.correlogramList}}
}
\keyword{dplot}

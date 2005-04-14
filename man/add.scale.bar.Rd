\name{add.scale.bar}
\alias{add.scale.bar}
\title{Add a Scale Bar to a Phylogeny Plot}
\usage{
add.scale.bar(x = 0, y = 1, length = NULL, ...)
}
\arguments{
  \item{x}{x location of the bar.}
  \item{y}{y location of the bar.}
  \item{length}{a numeric value giving the length of the scale bar. If
    none is supplied, a value is calculated from the data.}
  \item{...}{further arguments to be passed to \code{text}.}
}
\description{
  This function adds a horizontal bar giving the scale of the branch
  lengths to a plot of a phylogenetic tree on the current graphical
  device.
}
\details{
  As from version 1.4 of ape, the options of this function have been
  redefined, and have now default values. By default, the scale bar is
  plotted on the left bottom corner of the plot.

  The further arguments (\code{...}) are used to format the text. They
  may be \code{font}, \code{cex}, \code{col}, and so on (see examples
  below, and the help page on \code{\link[graphics]{text}}).

  The function \code{\link[graphics]{locator}}  may be used to
  determine the \code{x} and \code{y} arguments.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{axisPhylo}},
  \code{\link[graphics]{locator}}
}
\examples{
tr <- rtree(10)
layout(matrix(1:2, 2, 1))
plot(tr)
add.scale.bar()
plot(tr)
add.scale.bar(cex = 0.7, font = 2, col = "red")
layout(matrix(1, 2, 1))
}
\keyword{aplot}

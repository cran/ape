\name{add.scale.bar}
\alias{add.scale.bar}
\title{Add a Scale Bar to a Phylogeny Plot}
\usage{
add.scale.bar(x, y, phy, bar = NULL)
}
\arguments{
  \item{x}{the x location of the legend.}
  \item{y}{the y location of the legend.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{bar}{a numeric value giving the length of the scale bar. If none
    is supplied (the default), a value is calculated from the data.}
}
\description{
  This function adds a horizontal bar giving the scale of the branch
  lengths to a plot of a phylogenetic tree on the current graphical
  device.
}
\details{
  This function can be used in two ways: either an object of class
  \code{"phylo"} is given and its branch lengths are used to determine
  the length of the scale bar, or the argument \code{bar} is specified
  and this is used as the length of the bar. In the latter case, the
  argument \code{phy} may be omitted.
}
\note{
  Since the coordinate system is not displayed when plotting a tree, it
  may be useful to use the function \code{\link[base]{locator}} to
  determine the \code{x} and \code{y} arguments.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{plot.phylo}}, \code{\link[base]{locator}}
}
\keyword{hplot}

\name{evonet}
\alias{evonet}
\alias{print.evonet}
\alias{plot.evonet}
\alias{as.phylo.evonet}
\alias{as.networx.evonet}
\alias{as.network.evonet}
\alias{as.igraph.evonet}
\title{Evolutionary Networks}
\description{
  \code{evonet} builds a network from a tree of class
  \code{"phylo"}. There are \code{print} and \code{plot} methods as well
  as a few conversion functions.
}
\usage{
evonet(phy, from, to = NULL)
\method{print}{evonet}(x, ...)
\method{plot}{evonet}(x, col = "blue", lty = 1, lwd = 1, alpha = 0.5,
              arrows = 0, arrow.type = "classical", ...)
\method{as.phylo}{evonet}(x, ...)
\method{as.networx}{evonet}(x, weight = NA, ...)
\method{as.network}{evonet}(x, directed = TRUE, ...)
\method{as.igraph}{evonet}(x, directed = TRUE, use.labels = TRUE, ...)
}
\arguments{
  \item{phy, x}{an object of class \code{"phylo"} or \code{"evonet"}.}
  \item{from}{a vector (or a matrix if \code{to = NULL}) giving the node
    or tip numbers involved in the reticulations.}
  \item{to}{a vector of the same length than \code{from}.}
  \item{col, lty, lwd}{colors, line type and width of the reticulations
    (recycled if necessary).}
  \item{alpha}{a value between 0 and 1 specifying the transparency of
    the reticulations.}
  \item{arrows}{see \code{\link{fancyarrows}}.}
  \item{arrow.type}{idem.}
  \item{weight}{a numeric vector giving the weights for the
    reticulations when converting to the class \code{"networx"}
    (recycled or shortened if needed).}
  \item{directed}{a logical: should the network be considered as
    directed? \code{TRUE} by default.}
  \item{use.labels}{a logical specifying whether to use the tip and node
    labels when building the network of class \code{"igraph"}.}
  \item{\dots}{arguments passed to other methods.}
}
\details{
  \code{evonet} is a constructor function that checks the arguments.

  The classes \code{"networx"}, \code{"network"}, and \code{"igraph"}
  are defined in the packages \pkg{phangorn}, \pkg{network}, and
  \pkg{igraph}, respectively.
}
\value{
  an object of class \code{c("evonet", "phylo")} which is made of an
  object of class \code{"\link{phylo}"} plus an element
  \code{reticulation} coding additional edges among nodes and uses the
  same coding rules than the \code{edge} matrix.

  The conversion functions return an object of the appropriate class.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link[phangorn]{as.networx}} in package \pkg{phangorn}
}
\examples{
tr <- rcoal(5)
(x <- evonet(tr, 6:7, 8:9))
plot(x)
}
\keyword{manip}
\keyword{hplot}

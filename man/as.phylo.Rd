\name{as.phylo}
\alias{as.phylo}
\alias{as.phylo.hclust}
\alias{as.phylo.phylog}
\alias{as.hclust.phylo}
\title{Conversion Among Tree Objects}
\usage{
as.phylo(x, ...)
\method{as.phylo}{hclust}(x, ...)
\method{as.phylo}{phylog}(x, ...)
\method{as.hclust}{phylo}(x, ...)
}
\arguments{
  \item{x}{an object of the correct to be converted into another class.}
  \item{...}{further arguments to be passed to or from other methods.}
}
\description{
  \code{as.phylo} is a generic function which converts an object into a
  tree of class \code{"phylo"}. There are currently two methods for this
  generic for objects of class \code{"hclust"} and of class
  \code{"phylog"} (implemented in the package ade4).

  \code{as.hclust.phylo} is a method of the generic
  \code{\link[stats]{as.hclust}} which converts an object of class
  \code{"phylo"} into one of class \code{"hclust"}. This can used to
  convert an object of class \code{"phylo"} into one of class
  \code{"dendrogram"} (see examples).
}
\value{
  An object of class \code{"hclust"} or \code{"phylo"}.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link[stats]{hclust}}, \code{\link[stats]{as.hclust}},
  \code{\link[stats]{dendrogram}}, \code{\link[ade4]{phylog}}
}
\examples{
data(bird.orders)
hc <- as.hclust(bird.orders)
tr <- as.phylo(hc)
identical(bird.orders, tr) # FALSE, but...
all.equal(bird.orders, tr) # ... TRUE

### shows the three plots for tree objects:
dend <- as.dendrogram(hc)
layout(matrix(c(1:3, 3), 2, 2))
plot(bird.orders, font = 1)
plot(hc)
par(mar = c(8, 0, 0, 0)) # leave space for the labels
plot(dend)

### how to get (nearly) identical plots with
### plot.phylo and plot.dendrogram:
layout(matrix(1:2, 2, 1))
plot(bird.orders, font = 1, no.margin = TRUE)
par(mar = c(0, 0, 0, 8))
plot((dend), horiz = TRUE)
layout(matrix(1, 1, 1))
}
\keyword{manip}

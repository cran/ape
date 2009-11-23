\name{plot.ancestral}
\alias{plot.ancestral}
\title{Plot Ancestral Character Values on a Tree}
\description{
  Plot a phylogenetic tree with edge colors picked according to the
  corresponding node ancestral character value.
}
\usage{
\method{plot}{ancestral}(x, which = names(x$node.character), n.col = 10,
col.fun = function(n) rainbow(n, start = 0.4, end = 0),
plot.node.values = FALSE,
ask = prod(par("mfcol")) < length(which) && dev.interactive(),
...)
}
\arguments{
  \item{x}{An object of class 'ancestral'.}
  \item{which}{Which characters to plot. Can be a vecotr of names, or a
    vector of indices.}
  \item{n.col}{The number of colors to use in the gradient.}
  \item{col.fun}{the color function to use.}
  \item{plot.node.values}{Should character values used as node labels?}
  \item{ask}{Ask before each plot?}
  \item{\dots}{Further parameters to pass to the plot.phylo function.}
}
\details{
  This function produces one plot by selected ancestral character. It
  uses the plot.phylo function with particular arguments to display edge
  colors according to ancestral character values.
}
\author{Julien Dutheil \email{Julien.Dutheil@univ-montp2.fr}}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{evolve.phylo}}
}
\examples{
data(bird.orders)
x <- rep(0, 4)
names(x) <- c("A", "B", "C", "D")
anc <- evolve.phylo(bird.orders, x, 1)
plot(anc, edge.width = 3, plot.node.values = TRUE)
par(mfrow = c(2, 2), mar = c(5.5, 0, 0, 0))
plot(anc, edge.width = 3, type = "r")
}
\keyword{hplot}

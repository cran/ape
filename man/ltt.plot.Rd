\name{ltt.plot}
\alias{ltt.plot}
\title{Lineages Through Time Plot}
\usage{
ltt.plot(phy, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{...}{further arguments to be passed to \code{plot()} (see
    \code{Details:} on how to transform the axes).}
}
\description{
  This function plots, on the current graphical device, the minimum
  numbers of lineages through time from a phylogenetic tree.
}
\details{
  The \eqn{y}-axis can be log-transformed by adding the following
  option: \code{log = "y"}.
}
\references{
  Harvey, P. H., May, R. M. and Nee, S. (1994) Phylogenies without
  fossils. \emph{Evolution}, \bold{48}, 523--529.

  Nee, S., Holmes, E. C., Rambaut, A. and Harvey, P. H. (1995) Inferring
  population history from molecular phylogenies. \emph{Philosophical
    Transactions of the Royal Society of London. Series B. Biological
    Sciences}, \bold{349}, 25--31.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{skyline}}, \code{\link{branching.times}},
  \code{\link{birthdeath}}, \code{\link[base]{plot}} for the basic
  plotting function in R
}
\examples{
data(bird.families)
data(bird.orders)
opar <- par(mfrow = c(2, 1))
ltt.plot(bird.families)
title("Lineages Through Time Plot of the Bird Families")
ltt.plot(bird.families, log = "y")
title(main = "Lineages Through Time Plot of the Bird Families",
      sub = "(with logarithmic transformation of the y-axis)")
par(opar)

### to plot the tree and the LTT plot together
tmp <- bird.families
tmp$tip.label <- rep("", length(bird.families$tip.label))
tmp2 <- bird.orders
tmp2$tip.label <- rep("", length(bird.orders$tip.label))
layout(matrix(1:4, 2, 2))
plot(tmp)
ltt.plot(bird.families, main = "Bird families")
plot(tmp2)
ltt.plot(bird.orders, main = "Bird orders")
}
\keyword{hplot}

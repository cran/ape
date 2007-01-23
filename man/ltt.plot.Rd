\name{ltt.plot}
\alias{ltt.plot}
\alias{ltt.lines}
\alias{mltt.plot}
\title{Lineages Through Time Plot}
\usage{
ltt.plot(phy, xlab = "Time", ylab = "N", ...)
ltt.lines(phy, ...)
mltt.plot(phy, ..., dcol = TRUE, dlty = FALSE,
          legend = TRUE, xlab = "Time", ylab = "N")
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}; this could be an object
    of class \code{c("phylo", "multi.tree")} in the case of
    \code{mltt.plot}.}
  \item{xlab}{a character string (or a variable of mode character)
    giving the label for the x-axis (default is "Time").}
  \item{ylab}{idem for the y-axis (default is "N").}
  \item{...}{in the cases of \code{ltt.plot()} and \code{ltt.lines()}
    these are further (graphical) arguments to be passed to
    \code{plot()} or \code{lines()}, respectively (see \code{Details:}
    on how to transform the axes); in the case \code{mltt.plot()} these
    are additional trees to be plotted (see \code{Details:}).}
  \item{dcol}{a logical specifying whether the different curves should
    be differentiated with colors (default is \code{TRUE}).}
  \item{dlty}{a logical specifying whether the different curves should
    be differentiated with patterns of dots and dashes (default is
    \code{FALSE}).}
  \item{legend}{a logical specifying whether a legend should be
    plotted.}
}
\description{
  These functions plot, on the current graphical device, the minimum
  numbers of lineages through time from phylogenetic trees.
}
\details{
  \code{ltt.plot} does a simple lineages through time (LTT)
  plot. Additional arguments (\code{...}) may be used to change, for
  instance, the limits on the axes (with \code{xlim} and/or
  \code{ylim}) or other graphical settings (\code{col} for the color,
  \code{lwd} for the line thickness, \code{lty} for the line type may be
  useful; see \code{\link[graphics]{par}} for an exhaustive listing of
  graphical parameters). The \eqn{y}-axis can be log-transformed by
  adding the following option: \code{log = "y"}.

  \code{ltt.lines} adds a LTT curve to an existing plot. Additional
  arguments (\code{...}) may be used to change the settings of the added
  line. Of course, the settings of the already existing LTT plot cannot
  be altered this way.

  \code{mltt.plot} does a multiple LTT plot taking as arguments one or
  several trees. These trees may be given as objects of class
  \code{"phylo"} (single trees) or \code{c("phylo", "multi.tree")}
  (multiple trees). Any number of objects may be given. This function is
  mainly for exploratory analyses with the advantages that the axes are
  set properly to view all lines, and the legend is plotted by
  default. For more flexible settings of line drawings, it is probably
  better to combine \code{ltt.plot()} with successive calls of
  \code{ltt.lines()} (see \code{Examples:}).
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
  \code{\link{birthdeath}}, \code{\link{bd.ext}}, \code{\link{yule.cov}}
  \code{\link[graphics]{plot}} for the basic plotting function in R
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
layout(matrix(1:4, 2, 2))
plot(bird.families, show.tip.label = FALSE)
ltt.plot(bird.families, main = "Bird families")
plot(bird.orders, show.tip.label = FALSE)
ltt.plot(bird.orders, main = "Bird orders")
layout(matrix(1))
mltt.plot(bird.families, bird.orders)
### Generates 10 random trees with 23 tips:
TR <- replicate(10, rcoal(23), FALSE)
### Give names to each tree:
names(TR) <- paste("random tree", 1:10)
### And specify the class of the list so that mltt.plot()
### does not trash it!
class(TR) <- c("phylo", "multi.tree")
mltt.plot(TR)
### And now for something (not so) completely different:
ltt.plot(bird.orders, lwd = 2)
lapply(TR, ltt.lines, lty = 2) # same as:
### for (i in 1:10) ltt.lines(TR[[i]], lty = 2)
legend(-10, 5, lwd = c(2, 1), lty = c(1, 2), bty = "n",
       legend = c("Bird orders", "Random trees"))
}
\keyword{hplot}
\keyword{aplot}

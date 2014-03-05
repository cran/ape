\name{phydataplot}
\alias{phydataplot}
\alias{ring}
\title{Tree Annotation}
\description{
  \code{phydataplot} plots data on a tree in a way that adapts to the
  type of tree. \code{ring} does the same for circular trees.

  Both functions match the data with the labels of the tree.
}
\usage{
phydataplot(x, phy, style = "bars", offset = 1, ...)
ring(x, phy, style = "ring", offset = 1, ...)
}
\arguments{
  \item{x}{a vector, a factor, a matrix, or a data frame.}
  \item{phy}{the tree (which must be already plotted).}
  \item{style}{a character string specifying the type of graphics; can
    be abbreviated (see details).}
  \item{offset}{the space between the tips of the tree and the plot.}
  \item{\dots}{further arguments passed to the graphical functions.}
}
\details{
  The possible values for \code{style} are ``bars'', ``segments'',
  ``image'', or ``arrows'' for \code{phydataplot}, and ``ring'',
  ``segments'', or ``arrows'' for \code{ring}.
}
\note{
  For the moment, only rightwards trees are supported (does not apply to
  circular trees).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{nodelabels}},
  \code{\link{fancyarrows}}
}
\examples{
## demonstrates matching with names:
tr <- rcoal(n <- 10)
x <- 1:n
names(x) <- tr$tip.label
plot(tr, x.lim = 11)
phydataplot(x, tr)
## shuffle x but matching names with tip labels reorders them:
phydataplot(sample(x), tr, "s", lwd = 3, lty = 3)

## adapts to the tree:
plot(tr, "f", x.l = c(-11, 11), y.l = c(-11, 11))
phydataplot(x, tr, "s")

## leave more space with x.lim to show 2 barplots:
plot(tr, x.lim = 22)
phydataplot(x, tr, "b")
phydataplot(x, tr, "b", offset = 12, col = "yellow")

ts <- rcoal(N <- 100)
X <- rTraitCont(ts) # names are set
dd <- dist(X)
op <- par(mar = rep(0, 4))
plot(ts, x.lim = 10, cex = 0.4, font = 1)
phydataplot(as.matrix(dd), ts, "i", offset = 0.2)

par(xpd = TRUE, mar = op$mar)
co <- c("blue", "red"); l <- c(-2, 2)
X <- X + abs(min(X)) # move scale so X >= 0
plot(ts, "f", show.tip.label = FALSE, x.lim = l, y.lim = l, open.angle = 30)
phydataplot(X, ts, "s", col = co, offset = 0.05)
ring(X, ts, "ring", col = co, offset = max(X) + 0.1) # the same info as a ring

## as many rings as you want...
co <- c("blue", "yellow")
plot(ts, "r", show.tip.label = FALSE, x.l = c(-1, 1), y.l = c(-1, 1))
for (o in seq(0, 0.4, 0.2)) {
    co <- rev(co)
    ring(0.2, ts, "r", col = rep(co, each = 5), offset = o)
}

lim <- c(-5, 5)
co <- rgb(0, 0.4, 1, alpha = 0.1)
y <- seq(0.01, 1, 0.01)
plot(ts, "f", x.lim = lim, y.lim = lim, show.tip.label = FALSE)
ring(y, ts, offset = 0, col = co, lwd = 0.1)
for (i in 1:3) {
    y <- y + 1
    ring(y, ts, offset = 0, col = co, lwd = 0.1)
}

## rings can be in the background
plot(ts, "r", plot = FALSE)
ring(1, ts, "r", col = rainbow(100), offset = -1)
par(new = TRUE)
plot(ts, "r", font = 1, edge.color = "white")

## might be more useful:
co <- c("lightblue", "yellow")
plot(ts, "r", plot = FALSE)
ring(0.1, ts, "r", col = sample(co, size = N, rep = TRUE), offset = -.1)
par(new = TRUE)
plot(ts, "r", font = 1)

## if x is matrix:
tx <- rcoal(m <- 20)
X <- runif(m, 0, 0.5); Y <- runif(m, 0, 0.5)
X <- cbind(X, Y, 1 - X - Y)
rownames(X) <- tx$tip.label
plot(tx, x.lim = 6)
co <- rgb(diag(3))
phydataplot(X, tx, col = co)
## a variation:
plot(tx, show.tip.label = FALSE, x.lim = 5)
phydataplot(X, tx, col = co, offset = 0.05, border = NA)

plot(tx, "f", show.tip.label = FALSE, open.angle = 180)
ring(X, tx, col = co, offset = 0.05)
}
\keyword{aplot}

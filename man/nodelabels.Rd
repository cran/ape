\name{nodelabels}
\alias{nodelabels}
\alias{tiplabels}
\title{Labelling the Nodes and the Tips of a Tree}
\description{
  These functions add labels to or near the nodes or the tips of a tree
  using text or plotting symbols. The text can be framed.
}
\usage{
nodelabels(text, node, adj = c(0.5, 0.5), frame = "rect",
           pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
           col = "black", bg = "lightblue", ...)
tiplabels(text, tip, adj = c(0.5, 0.5), frame = "rect",
          pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
          col = "black", bg = "yellow", ...)
}
\arguments{
  \item{text}{a vector of mode character giving the text to be
    printed. Can be left empty.}
  \item{node}{a vector of mode numeric giving the numbers of the nodes
    where the text or the symbols are to be printed. Can be left empty.}
  \item{tip}{a vector of mode numeric giving the numbers of the tips
    where the text or the symbols are to be printed. Can be left empty.}
  \item{adj}{one or two numeric values specifying the horizontal and
    vertical, respectively, justification of the text. By default, the
    text is centered horizontally and vertically. If a single value is
    given, this alters only the horizontal position of the text.}
  \item{frame}{a character string specifying the kind of frame to be
    printed around the text. This must be one of "rect" (the default),
    "circle", "none", or any unambiguous abbreviation of these.}
  \item{pch}{a numeric giving the type of plotting symbol to be used;
    this is eventually recycled. See \code{\link[graphics]{par}} for R's
    plotting symbols. If \code{pch} is used, then \code{text} is
    ignored.}
  \item{thermo}{a numeric vector giving some proportions (values between
    0 and 1) for each node, or a numeric matrix giving some proportions
    (the rows must sum to one).}
  \item{pie}{same than \code{thermo}.}
  \item{piecol}{a list of colours (given as a character vector) to be
    used by \code{thermo} or \code{pie}; if left \code{NULL}, a series
    of colours given by the function \code{rainbow} is used.}
  \item{col}{a character string giving the color to be used for the
    text or the plotting symbols; this is eventually recycled.}
  \item{bg}{a character string giving the color to be used for the
    background of the text frames or of the plotting symbols if it
    applies; this is eventually recycled.}
  \item{\dots}{further arguments passed to the \code{text} or
    \code{points} functions (e.g. \code{cex} to alter the size of the
    text or the symbols, or \code{font} for the text; see the examples
    below).}
}
\details{
  These two functions have the same optional arguments and the same
  functioning.

  If the arguments \code{text} is missing and \code{pch} and
  \code{thermo} are left as \code{NULL}, then the numbers of the nodes
  (or of the tips) are printed.

  If \code{node} (or \code{tip}) is missing, then the text or the
  symbols are printed on all nodes (or tips).

  The option \code{cex} can be used to change the size of all types of
  labels.

  A simple call of these functions with no arguments (e.g.,
  \code{nodelabels()}) prints the numbers of all nodes (or tips).

  In the case of \code{tiplabels}, it will be useful to play with the
  options \code{x.lim} and \code{label.offset} (and possibly
  \code{show.tip.label}) of \code{plot.phylo} in most cases (see the
  examples).
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}, Ben Bolker
  \email{bolker@zoo.ufl.edu}, and Jim Lemon}
\seealso{
  \code{\link{plot.phylo}}
}
\examples{
tr <- read.tree(text = "((Homo,Pan),Gorilla);")
plot(tr)
nodelabels("7.3 Ma", 1, frame = "r", bg = "yellow", adj = 0)
nodelabels("5.4 Ma", 2, frame = "c", bg = "tomato", font = 3)

data(bird.orders)
plot(bird.orders, use.edge.length = FALSE, font = 1)
bs <- round(runif(22, 90, 100), 0) # some imaginary bootstrap values
bs2 <- round(runif(22, 90, 100), 0)
bs3 <- round(runif(22, 90, 100), 0)
nodelabels(bs, adj = 1.2)
nodelabels(bs2, adj = -0.2, bg = "yellow")

### something more classical
plot(bird.orders, use.edge.length = FALSE, font = 1)
nodelabels(bs, adj = -0.2, frame = "n", cex = 0.8)
nodelabels(bs2, adj = c(1.2, 1), frame = "n", cex = 0.8)
nodelabels(bs3, adj = c(1.2, -0.2), frame = "n", cex = 0.8)

### the same but we play with the font
plot(bird.orders, use.edge.length = FALSE, font = 1)
nodelabels(bs, adj = -0.2, frame = "n", cex = 0.8, font = 2)
nodelabels(bs2, adj = c(1.2, 1), frame = "n", cex = 0.8, font = 3)
nodelabels(bs3, adj = c(1.2, -0.2), frame = "n", cex = 0.8)

plot(bird.orders, "c", use.edge.length = FALSE, font = 1)
nodelabels(thermo = runif(22), cex = .8)

plot(bird.orders, "u", FALSE, font = 1, lab4ut = "a")
nodelabels(cex = .75, bg = "yellow")

### representing two characters at the tips (you could have as many
### as you want)
plot(bird.orders, "c", FALSE, font = 1, label.offset = 3,
     x.lim = 31, no.margin = TRUE)
tiplabels(pch = 21, bg = gray(1:23/23), cex = 2, adj = 1.4)
tiplabels(pch = 19, col = c("yellow", "red", "blue"), adj = 2.5, cex = 2)
### Some random data to compare piecharts and thermometres:
tr <- rcoal(15)
x <- runif(14, 0, 0.33)
y <- runif(14, 0, 0.33)
z <- runif(14, 0, 0.33)
x <- cbind(x, y, z, 1 - x - y - z)
layout(matrix(1:2, 1, 2))
plot(tr, "c", FALSE, no.margin = TRUE)
nodelabels(pie = x, cex = 1.3)
text(4.5, 15, "Are you \"pie\"...", font = 4, cex = 1.5)
plot(tr, "c", FALSE, no.margin = TRUE)
nodelabels(thermo = x, col = rainbow(4), cex = 1.3)
text(4.5, 15, "... or \"thermo\"?", font = 4, cex = 1.5)
layout(matrix(1))
}
\keyword{aplot}

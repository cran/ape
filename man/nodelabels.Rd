\name{nodelabels}
\alias{nodelabels}
\title{Labelling the Nodes of a Tree}
\description{
  This function adds labels to or near the nodes of a tree using text or
  plotting symbols. The text can be framed.
}
\usage{
nodelabels(text, node, adj = c(0.5, 0.5), frame = "rect",
           pch = NULL, col = "black", bg = "lightblue", ...)
}
\arguments{
  \item{text}{a vector of mode character giving the text to be
    printed. Can be left empty.}
  \item{node}{a vector of mode numeric giving the numbers of the nodes
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
  \item{col}{a character string giving the color to be used for the
    plotting symbols; this is eventually recycled.}
  \item{bg}{a character string giving the color to be used for the
    background of the text frames or of the plotting symbols if it
    applies; this is eventually recycled.}
  \item{\dots}{further arguments passed to the \code{text} or
    \code{points} functions (e.g. \code{cex} to alter the size of the
    text or the symbols, or \code{font} for the text; see the examples
    below).}
}
\details{
  If the arguments \code{text} is missing and \code{pch} is left as
  \code{NULL}, then the numbers of the nodes are printed.

  If \code{node} is missing, then the text or the symbols are printed on
  all nodes.

  A simple call of the function with no arguments
  (i.e. \code{nodelabels()}) prints the numbers of all nodes.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
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

plot(bird.orders, "u", FALSE, font = 1, lab4ut = "a")
nodelabels(cex = .75, bg = "yellow")
}
\keyword{aplot}

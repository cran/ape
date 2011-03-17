\name{image.DNAbin}
\alias{image.DNAbin}
\title{Plot of DNA Sequence Alignement}
\description{
  This function plots an image of an alignment of nucleotide sequences.
}
\usage{
\method{image}{DNAbin}(x, what, col, bg = "white", xlab = "", ylab = "",
       show.labels = TRUE, cex.lab = 1, legend = TRUE, ...)
}
\arguments{
  \item{x}{a matrix of DNA sequences (class \code{"DNAbin"}).}
  \item{what}{a vector of characters specifying the bases to
    visualize. If missing, this is set to ``a'', ``g'', ``c'', ``t'',
    ``n'', and ``-'' (in this order).}
  \item{col}{a vector of colours. If missing, this is set to ``red'',
    ``yellow'', ``green'', ``blue'', ``grey'', and ``black''. If it is
    shorter (or longer) than \code{what}, it is recycled (or shortened).}
  \item{bg}{the colour used for nucleotides whose base is not among
    \code{what}.}
  \item{xlab}{the label for the \emph{x}-axis; none by default.}
  \item{ylab}{Idem for the \emph{y}-axis. Note that by default, the
    labels of the sequences are printed on the \emph{y}-axis (see next option).}
  \item{show.labels}{a logical controlling whether the sequence labels
    are printed (\code{TRUE} by default).}
  \item{cex.lab}{a single numeric controlling the size of the sequence labels.
    Use \code{cex.axis} to control the size of the annotations on the \emph{x}-axis.}
  \item{legend}{a logical controlling whether the legend is plotted
    (\code{TRUE} by default).}
  \item{\dots}{further arguments passed to
    \code{\link[graphics]{image.default}} (e.g., \code{cex.axis}).}
}
\details{
  The idea of this function is to allow fleixble plotting and colouring
  of a nucleotide alignment. By default, the most common bases (a, g, c,
  t, and n) and alignment gap are plotted using a standard colour
  scheme.

  It is possible to plot only one base specified as \code{what} with a
  chosen colour: this might be useful to check, for instance, the
  distribution of alignment gaps (\code{image(x, "-")}) or missing data
  (see examples).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{DNAbin}}, \code{\link{del.gaps}}, \code{\link{clustal}},
  \code{\link[graphics]{grid}}
}
\examples{
data(woodmouse)
image(woodmouse)
image(woodmouse, "n", "blue") # show missing data
image(woodmouse, c("g", "c"), "green") # G+C
par(mfcol = c(2, 2))
### barcoding style:
for (x in c("a", "g", "c", "t"))
    image(woodmouse, x, "black", cex.lab = 0.5, cex.axis = 0.7)
par(mfcol = c(1, 1))
### zoom on a portion of the data:
image(woodmouse[11:15, 1:50], c("a", "n"), c("blue", "grey"))
grid(50, 5, col = "black")
### see the guanines on a black background:
image(woodmouse, "g", "yellow", "black")
}
\keyword{hplot}


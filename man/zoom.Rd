\name{zoom}
\alias{zoom}
\title{Zoom on a Portion of a Phylogeny}
\description{
  This function plots simultaneously a whole phylogenetic tree
  (supposedly large) and a portion of it.
}
\usage{
zoom(phy, focus, type = "phylogram", use.edge.length = TRUE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{focus}{a vector, either numeric or character, specifying the
    tips to be focused on.}
  \item{type}{a character string specifying the type of phylogeny to be
    drawn; it must be one of "phylogram" (the default), "cladogram",
    "unrooted", or any unambiguous abbreviation of these.}
  \item{use.edge.length}{a logical indicating whether to use the edge
    lengths of the phylogeny to draw the branches (the default) or not
    (if \code{FALSE}). This option has no effect if the object of class
    \code{"phylo"} has no `edge.length' element.}
}
\details{
  This is an experimental function aimed at exploring very large
  trees. The main argument is a phylogenetic tree, and the second one is
  a vector specifying the tips to be focused on. The latter can be
  either numeric and is thus taken as the index of the tip labels, or
  character in which case it is taken as the corresponding tip labels.

  The plot consists of two trees: the whole tree is plotted on the
  left-hand side in a narrower sub-window (about a quarter of the
  device) without any tip labels; the sub-tree consisting of the tips
  in `focus' is extracted and plotted on the right-hand side. The tips
  in `focus' are marked on the left-hand side tree with a blue ``*''.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{drop.tip}},
  \code{\link[graphics]{layout}}
}
\examples{
\dontrun{
data(chiroptera)
zoom(chiroptera, 1:20)
zoom(chiroptera, grep("Plecotus", chiroptera$tip.label))
}
}
\keyword{hplot}

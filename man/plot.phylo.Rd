\name{plot.phylo}
\alias{plot.phylo}
\title{Plot Phylogenies}
\usage{
\method{plot}{phylo}(x, show.node.label = FALSE, edge.color = NULL,
                     edge.width = NULL, font = 3, adj = 0, srt = 0,
                     no.margin = FALSE, label.offset = NULL,
                     underscore = FALSE, x.lim = NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{show.node.label}{a logical indicating whether to show the node
    labels on the phylogeny (defaults to \code{"FALSE"}, i.e. the labels
    are not shown).}
  \item{edge.color}{a vector of mode character giving the colours used
    to draw the branches of the plotted phylogeny. These are taken to be
    in the same order than the component \code{edge} of \code{phy}. If
    fewer colours are given than the length of \code{edge}, then the
    colours are recycled.}
  \item{edge.width}{a numeric vector giving the width of the branches of
    the plotted phylogeny. These are taken to be in the same order than
    the component \code{edge} of \code{phy}. If fewer widths are given
    than the length of \code{edge}, then these are recycled.}
  \item{font}{an integer specifying the type of font for the labels: 1
    (plain text), 2 (bold), 3 (italic, the default), or 4 (bold italic).}
  \item{adj}{a numeric specifying the justification of the text strings
    of the labels: 0 (left-justification, the default), 0.5 (centering),
    or 1 (right-justification).}
  \item{srt}{a numeric giving how much the labels are rotated in degrees
    (negative values are allowed resulting in clock-like rotation).}
  \item{no.margin}{a logical. If \code{TRUE}, the margins are set to
    zero and the plot uses all the space of the device (note that this
    was the behaviour of \code{plot.phylo} up to version 0.2-1 of `ape'
    with no way to modify it by the user, at least easily).}
  \item{label.offset}{a numeric giving the space between the nodes and
    the tips of the phylogeny and their corresponding labels. If
    \code{NULL}, this is computed with respect to various parameters
    such as the string lengths of the labels.}
  \item{underscore}{a logical specifying whether the underscores in tip
    labels should be written as spaces (the default) or left as are (if
    \code{TRUE}.}
  \item{x.lim}{a numeric giving the limit of the x-axis. If \code{NULL},
    this is computed with respect to various parameters such as the
    string lengths of the labels and the branch lengths.}
  \item{...}{further arguments to be passed to \code{plot()}.}
}
\description{
  This function plots a phylogenetic tree on the current graphical device.
}
\details{
  The formatting of the labels of both the nodes and the tips is the
  same.

  The new default behaviour is to use the current margins rather
  than to set them to zero. If \code{no.margin = TRUE}, the margins are
  set to zero and are not restored after plotting the tree, so that the
  user can access the coordinates system of the plot.

  If the arguments \code{label.offset} and \code{x.lim} are not
  specified by the user, they are determined roughly by the
  function. This may not always give a nice result: the user may check
  these values with the (invisibly) returned list (see ``Value:'').
}
\value{
  The function returns invisibly a list with the following components
  which values are those used for the current plot:

  \item{show.node.label}
  \item{edge.color}
  \item{edge.width}
  \item{font}
  \item{adj}
  \item{srt}
  \item{no.margin}
  \item{label.offset}
  \item{x.lim}  
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link{add.scale.bar}},
  \code{\link[base]{plot}} for the basic plotting function in R
}
\examples{
### An extract from Sibley and Ahlquist (1990)
cat("(((Strix_aluco:4.2,Asio_otus:4.2):3.1,",
   "Athene_noctua:7.3):6.3,Tyto_alba:13.5);",
   file = "ex.tre", sep = "\n")
tree.owls <- read.tree("ex.tre")
plot(tree.owls)
unlink("ex.tre") # delete the file "ex.tre"

data(xenarthra)
op <- par()
\dontrun{par(ask = TRUE)}
plot(xenarthra)
### remove the margins and put the tip labels closer
### to the tips
plot(xenarthra, no.margin = TRUE, label.offset = .005)
### use a smaller font size
par(cex = 0.8)
plot(xenarthra, no.margin = TRUE, label.offset = .005)
par(op)

data(bird.orders)
### using random colours and thichness
plot(bird.orders,
     edge.color = sample(colors(), length(bird.orders$edge)/2),
     edge.width = sample(1:10, length(bird.orders$edge)/2, replace = TRUE))
### rainbow colouring...
plot(bird.orders,
     edge.color = sample(c("red", "yellow", "blue", "green", "orange", "purple"), length(bird.orders$edge)/2, replace = TRUE),
     edge.width = sample(1:10, length(bird.orders$edge)/2, replace = TRUE))
}
\keyword{hplot}

\name{plot.phylo}
\alias{plot.phylo}
\title{Plot Phylogenies}
\usage{
\method{plot}{phylo}(x, type = "phylogram", use.edge.length = TRUE,
             node.pos = NULL, show.node.label = FALSE,
             edge.color = NULL, edge.width = NULL, font = 3,
             adj = 0, srt = 0, no.margin = FALSE, root.edge = FALSE,
             label.offset = NULL, underscore = FALSE, x.lim = NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{type}{a character string specifying the type of phylogeny to be
    drawn; it must be one of "phylogram" (the default), "cladogram",
    "unrooted", or any unambiguous abbreviation of these.}
  \item{use.edge.length}{a logical indicating whether to use the edge
    lengths of the phylogeny to draw the branches (the default) or not
    (if \code{FALSE}). This option has no effect if the object of class
    \code{"phylo"} has no `edge.length' element.}
  \item{node.pos}{a numeric taking the value 1 or 2 which specifies the
    vertical position of the nodes with respect to their descendants. If
    \code{NULL} (the default), then the value is determined in relation
    to `type' and `use.edge.length' (see details).}  
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
  \item{root.edge}{a logical indicating whether to draw the root edge
    (defaults to FALSE); this has no effect if `use.edge.length = FALSE'
    or if `type = "unrooted"'.}
  \item{label.offset}{a numeric giving the space between the nodes and
    the tips of the phylogeny and their corresponding labels. If
    \code{NULL}, this is computed with respect to various parameters
    such as the string lengths of the labels.}
  \item{underscore}{a logical specifying whether the underscores in tip
    labels should be written as spaces (the default) or left as are (if
    \code{TRUE}).}
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

  The option `node.pos' allows the user to alter the vertical position
  (i.e. ordinates) of the nodes. If \code{node.pos = 1}, then the
  ordinate of a node is the mean of the ordinates of its direct
  descendants (nodes and/or tips). If \code{node.pos = 2}, then the
  ordinate of a node is the mean of the ordinates of all the tips of
  which it is the ancestor. If \code{node.pos = NULL} (the default),
  then its value is determined with respect to other options: if
  \code{type = "phylogram"} then `node.pos = 1'; if \code{type =
    "cladogram"} and \code{use.edge.length = FALSE} then `node.pos = 2';
  if \code{type = "cladogram"} and \code{use.edge.length = TRUE} then
  `node.pos = 1'.

  If the arguments \code{label.offset} and \code{x.lim} are not
  specified by the user, they are determined roughly by the
  function. This may not always give a nice result: the user may check
  these values with the (invisibly) returned list (see ``Value:'').
}
\note{
  The placement of the tip labels in unrooted trees are not very
  nice. This will be sorted out in future versions.
}
\value{
  The function returns invisibly a list with the following components
  which values are those used for the current plot:

  \item{type}
  \item{use.edge.length}
  \item{node.pos}
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

### Show the types of trees.
layout(matrix(1:6, 3, 2))
plot(tree.owls, main = "With branch lengths")
plot(tree.owls, type = "c")
plot(tree.owls, type = "u")
plot(tree.owls, use.edge.length = FALSE, main = "Without branch lengths")
plot(tree.owls, type = "c", use.edge.length = FALSE)
plot(tree.owls, type = "u", use.edge.length = FALSE)
layout(matrix(1))

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
plot(xenarthra, type = "c", no.margin = TRUE, label.offset = .005, use.edge.length = FALSE)
par(op)

data(bird.orders)
### using random colours and thichness
plot(bird.orders,
     edge.color = sample(colors(), length(bird.orders$edge)/2),
     edge.width = sample(1:10, length(bird.orders$edge)/2, replace = TRUE))
title("Random colours and branch thickness")
### rainbow colouring...
X <- c("red", "orange", "yellow", "green", "blue", "purple")
plot(bird.orders,
     edge.color = sample(X, length(bird.orders$edge)/2, replace = TRUE),
     edge.width = sample(1:10, length(bird.orders$edge)/2, replace = TRUE))
title("Rainbow colouring")
plot(bird.orders, type = "c", use.edge.length = FALSE,
     edge.color = sample(X, length(bird.orders$edge)/2, replace = TRUE),
     edge.width = rep(5, length(bird.orders$edge)/2))
segments(rep(0, 6), 6.5:1.5, rep(2, 6), 6.5:1.5, lwd = 5, col = X)
text(rep(2.5, 6), 6.5:1.5, paste(X, "..."), adj = 0)
title("Character mapping is possible...")
}
\keyword{hplot}

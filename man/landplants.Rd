\name{landplants}
\alias{landplants}
\alias{landplants.newick}

\title{Gene Tree of 36 Landplant rbcL Sequences}

\description{
  This data set describes a gene tree estimated from 36 landplant
  \emph{rbc}L sequences. 
}

\usage{
data(landplants.newick)
}

\format{
  \code{landplants.newick} is a string with the tree in Newick format.
}

\source{
  This tree is described in Sanderson (1997) and is also  a
  data example in the software package r8s
  (\url{http://ginger.ucdavis.edu/r8s/}).
}

\seealso{
\code{\link{chronogram}}, \code{\link{ratogram}}, \code{\link{NPRS.criterion}}.
}

\references{
  Sanderson, M. J. (1997) A nonparametric approach to estimating
    divergence times in the absence of rate constancy. \emph{Molecular
    Biology and Evolution}, \bold{14}, 1218--1231.
}

\examples{
library(ape)

# example tree in NH format (a string)
data("landplants.newick") 
landplants.newick

# get corresponding phylo object
tree.landplants <- read.tree(text = landplants.newick)

# plot tree
plot(tree.landplants, label.offset = 0.001)
}
\keyword{datasets}


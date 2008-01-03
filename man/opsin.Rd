\name{opsin}
\alias{opsin}
\alias{opsin.newick}
\title{Gene Tree of 32 opsin Sequences}
\description{
  This data set describes a gene tree estimated from 32 opsin
  sequences.
}
\usage{
data(opsin.newick)
}
\format{
  \code{opsin.newick} is a string with the tree in Newick format.
}
\source{
  This tree is described in Misawa and Tajima (2000) as an example
  for application of the Klastorin (1982) classification method.
}
\seealso{
\code{\link{klastorin}}.
}
\references{
  Misawa, K. (2000) A simple method for classifying genes and a bootstrap
  test for classifications. \emph{Molecular
    Biology and Evolution}, \bold{17}, 1879--1884.
}
\examples{
# example tree in NH format (a string)
data("opsin.newick")
opsin.newick

# get corresponding phylo object
tree.opsin <- read.tree(text = opsin.newick)

# plot tree
plot(tree.opsin, label.offset = 0.01)
}
\keyword{datasets}

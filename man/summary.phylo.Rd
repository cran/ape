\name{summary.phylo}
\alias{summary.phylo}
\title{Print Summary of a Phylogeny}
\usage{
\method{summary}{phylo}(object, \dots)
}
\arguments{
  \item{object}{an object of class \code{"phylo"}.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  This function prints a compact summary of a phylogenetic tree (an
  object of class \code{"phylo"}).
}
\details{
  The summary includes the numbers of tips and of nodes, summary
  statistics of the branch lengths (if they are available) with mean,
  variance, minimum, first quartile, median, third quartile, and
  maximum, listing of the first ten tip labels, and (if available) of
  the first ten node labels. It is also printed whether some of these
  optional elements (branch lengths, node labels, and root edge) are not
  found in the tree.

  This function simply prints its results on the standard output and is
  not meant for programming.
}
\value{
  A NULL value is returned, the results are simply printed.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link[base]{summary}} for the generic R
  function
}
\examples{
data(bird.families)
summary(bird.families)
}
\keyword{manip}

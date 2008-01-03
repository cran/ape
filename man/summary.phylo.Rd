\name{summary.phylo}
\alias{summary.phylo}
\alias{Ntip}
\alias{Nnode}
\alias{Nedge}
\title{Print Summary of a Phylogeny}
\usage{
\method{summary}{phylo}(object, \dots)
Ntip(phy)
Nnode(phy, internal.only = TRUE)
Nedge(phy)
}
\arguments{
  \item{object, phy}{an object of class \code{"phylo"}.}
  \item{...}{further arguments passed to or from other methods.}
  \item{internal.only}{a logical indicating whether to return the number
    of internal nodes only (the default), or of internal and terminal
    (tips) nodes (if \code{FALSE}).}
}
\description{
  The first function prints a compact summary of a phylogenetic tree (an
  object of class \code{"phylo"}). The three other functions return the
  number of tips, nodes, or edges, respectively.
}
\details{
  The summary includes the numbers of tips and of nodes, summary
  statistics of the branch lengths (if they are available) with mean,
  variance, minimum, first quartile, median, third quartile, and
  maximum, listing of the first ten tip labels, and (if available) of
  the first ten node labels. It is also printed whether some of these
  optional elements (branch lengths, node labels, and root edge) are not
  found in the tree.

  If the tree was estimated by maximum likelihood with
  \code{\link{mlphylo}}, a summary of the model fit and the parameter
  estimated is printed.

  \code{summary} simply prints its results on the standard output and is
  not meant for programming.
}
\value{
  A NULL value in the case of \code{summary}, a single numeric value for
  the three other functions.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link[base]{summary}} for the generic R
  function
}
\examples{
data(bird.families)
summary(bird.families)
Ntip(bird.families)
Nnode(bird.families)
Nedge(bird.families)
}
\keyword{manip}

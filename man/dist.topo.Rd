\name{dist.topo}
\alias{dist.topo}
\title{Topological Distances Between Two Trees}
\usage{
dist.topo(x, y, method = "PH85")
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{y}{an object of class \code{"phylo"}.}
  \item{method}{a character string giving the method to be used: either
    \code{"PH85"}, or \code{"score"}.}
}
\description{
  This function computes the topological distance between two
  phylogenetic trees using different methods.
}
\value{
  a single numeric value.
}
\details{
  Two methods are available: the one by Penny and Hendy (1985), and the
  branch length score by Kuhner and Felsenstein (1994). The trees are
  always considered as unrooted.

  The topological distance is defined as twice the number of internal
  branches defining different bipartitions of the tips (Penny and Hendy
  1985). Rzhetsky and Nei (1992) proposed a modification of the original
  formula to take multifurcations into account.

  The branch length score may be seen as similar to the previous
  distance but taking branch lengths into account. Kuhner and
  Felsenstein (1994) proposed to calculate the square root of the sum of
  the squared differences of the (internal) branch lengths defining
  similar bipartitions (or splits) in both trees.
}
\note{
  The geodesic distance of Billera et al. (2001) has been disabled.
}
\references{
  Billera, L. J., Holmes, S. P. and Vogtmann, K. (2001) Geometry of the
  space of phylogenetic trees. \emph{Advances in Applied Mathematics},
  \bold{27}, 733--767.

  Kuhner, M. K. and Felsenstein, J. (1994) Simulation comparison of
  phylogeny algorithms under equal and unequal evolutionary rates.
  \emph{Molecular Biology and Evolution}, \bold{11}, 459--468.

  Nei, M. and Kumar, S. (2000) \emph{Molecular evolution and
  phylogenetics}. Oxford: Oxford University Press.

  Penny, D. and Hendy, M. D. (1985) The use of tree comparison
  metrics. \emph{Systemetic Zoology}, \bold{34}, 75--82.

  Rzhetsky, A. and Nei, M. (1992) A simple method for estimating and
  testing minimum-evolution trees. \emph{Molecular Biology and
    Evolution}, \bold{9}, 945--967.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{read.tree}} to read tree files in Newick format,
  \code{\link{cophenetic.phylo}}, \code{\link{prop.part}}
}
\examples{
ta <- rtree(30)
tb <- rtree(30)
dist.topo(ta, ta) # = 0
dist.topo(ta, tb) # This is unlikely to be 0 !
}
\keyword{manip}

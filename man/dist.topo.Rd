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
    \code{"PH85"}, or \code{"BHV01"}.}
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
  one by Billera et al. (2001).

  The topological distance is defined as twice the number of internal
  branches defining different bipartitions of the tips (Penny and Hendy
  1985). Rzhetsky and Nei (1992) proposed a modification of the original
  formula to take multifurcations into account.

  Billera et al. (2001) developed a distance from the geometry of a tree
  space. The distance between two trees can be seen as the sum of the
  branch lengths that need be erased to have two similar trees.
}
\author{
  Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}
}
\references{
  Billera, L. J., Holmes, S. P. and Vogtmann, K. (2001) Geometry of the
  space of phylogenetic trees. \emph{Advances in Applied Mathematics},
  \bold{27}, 733--767.

  Nei, M. and Kumar, S. (2000) \emph{Molecular evolution and
  phylogenetics}. Oxford: Oxford University Press.

  Penny, D. and Hendy, M. D. (1985) The use of tree comparison
  metrics. \emph{Systemetic Zoology}, \bold{34}, 75--82.

  Rzhetsky, A. and Nei, M. (1992) A simple method for estimating and
  testing minimum-evolution trees. \emph{Molecular Biology and
    Evolution}, \bold{9}, 945--967.
}
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

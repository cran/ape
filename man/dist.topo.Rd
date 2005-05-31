\name{dist.topo}
\alias{dist.topo}
\title{Topological Distances Between Two Trees}
\usage{
dist.topo(x, y)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{y}{an object of class \code{"phylo"}.}
}
\description{
  This function computes the topological distance between two
  phylogenetic trees using the formula in Nei and Kumar (2000).
}
\value{
  a single numeric value.
}
\details{
  The topological distance is defined as twice the number of internal
  branches defining different bipartitions of the tips (Penny and Hendy
  1985). Rzhetsky and Nei (1992) proposed a modification of the original
  formula to take multifurcations into account.
}
\author{
  Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}
}
\references{
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
  \code{\link{dist.phylo}}, \code{\link{prop.part}}
}
\examples{
ta <- rtree(30)
tb <- rtree(30)
dist.topo(ta, ta) # = 0
dist.topo(ta, tb) # This is unlikely to be 0 !
}
\keyword{manip}

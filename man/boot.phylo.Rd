\name{boot.phylo}
\alias{boot.phylo}
\alias{prop.part}
\alias{prop.clades}
\title{Tree Bipartition and Bootstrapping Phylogenies}
\usage{
boot.phylo(phy, x, FUN, B = 100, block = 1)
prop.part(...)
prop.clades(phy, ..., part = NULL)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{x}{a taxa (rows) by characters (columns) matrix; this may be
    presented as a list.}
  \item{FUN}{the function used to estimate \code{phy} (see details).}
  \item{B}{the number of bootstrap replicates.}
  \item{block}{the number of columns in \code{x} that will be resampled
    together (see details).}
  \item{...}{either (i) a single object of class \code{"phylo"}, (ii) a
    series of such objects separated by commas, or (iii) a list
    containing such objects.}
  \item{part}{a list of partitions as returned by \code{prop.part}; if
    this is used then \code{...} is ignored.}
}
\description{
  These functions analyse bipartitions found in a series of trees.

  \code{prop.part} counts the number of bipartitions found in a series
  of trees given as \code{...}.

  \code{prop.clades} counts the number of times the bipartitions present
  in \code{phy} are present in a series of trees given as \code{...} or
  in the list previously computed and given with \code{part}.

  \code{boot.phylo} performs a bootstrap analysis.
}
\details{
  The argument \code{FUN} in \code{boot.phylo} must be the function used
  to estimate the tree from the original (resampled) data matrix. Thus,
  if the tree is estimated with neighbor-joining (see \code{nj}), one
  maybe wants something like \code{FUN = function(xx) nj(dist.dna(xx))}.

  \code{block} in \code{boot.phylo} specifies the number of columns to
  be resampled altogether. For instance, if one wants to resample at the
  codon-level, then \code{block = 3} must be used.
}
\value{
  \code{prop.part} returns a list with an attribute \code{"number"}. The
  elements of this list are the observed clades, and the attribute their
  respective numbers.

  \code{prop.clades} and \code{boot.phylo} returns a numeric vector
  which \emph{i}th element is the number associated to the \emph{i}th
  node of \code{phy}.
}
\references{
  Efron, B., Halloran, E. and Holmes, S. (1996) Bootstrap confidence
  levels for phylogenetic trees. \emph{Proceedings of the National
    Academy of Sciences USA}, \bold{93}, 13429--13434.

  Felsenstein, J. (1985) Confidence limits on phylogenies: an approach
  using the bootstrap. \emph{Evolution}, \bold{39}, 783--791.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{dist.topo}}
}
\keyword{manip}
\keyword{htest}

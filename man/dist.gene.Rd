\name{dist.gene}
\alias{dist.gene}
\alias{dist.gene.pairwise}
\alias{dist.gene.percentage}
\title{Pairwise Distances from Genetic Data}
\usage{
dist.gene(x, method = "pairwise", variance = FALSE)
dist.gene.pairwise(x, variance = FALSE)
dist.gene.percentage(x, variance = FALSE)
}
\arguments{
  \item{x}{a matrix or a data frame.}
  \item{method}{a character string specifying the method used to compute
    the distances; only two choices are available: \code{"pairwise"},
    and \code{"percentage"}.}
  \item{variance}{a logical, indicates whether the variance of the
    distances should be returned (default to FALSE).}
}
\description{
  These functions compute a matrix of distances between pairs of
  individuals from a matrix or a data frame of genetic data.
}
\details{
  This function is meant to be very general and accepts different kinds
  of data (alleles, haplotypes, DNA sequences, and so on). The rows of
  the data matrix represent the individuals, and the columns the loci.

  In the case of the pairwise method, the distance \eqn{d} between two
  individuals is the number of loci for which they differ, and the
  associated variance is \eqn{d(L - d)/L}, where \eqn{L} is the number
  of loci.

  In the case of the percentage method, this distance is divided by \eqn{L},
  and the associated variance is \eqn{d(1 - d)/L}.

  For more elaborate distances with DNA sequences, see the function
  \code{dist.dna}.
}
\value{
  either a numeric matrix with possibly the names of the individuals (as
  given by the rownames of the argument \code{x}) as colnames and rownames
  (if \code{variance = FALSE}, the default), or a list of two matrices names
  \code{distances} and \code{variance}, respectively (if \code{variance =
    TRUE}).
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{dist.dna}}, \code{\link{cophenetic.phylo}}
}
\keyword{manip}


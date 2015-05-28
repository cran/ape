\name{del.gaps}
\alias{del.gaps}
\alias{del.colgapsonly}
\title{Delete Alignment Gaps in DNA Sequences}
\description{
  These functions remove indel gaps (\code{"-"}) in a sample of DNA
  sequences.
}
\usage{
del.gaps(x)
del.colgapsonly(x)
}
\arguments{
  \item{x}{a matrix, a list, or a vector containing the DNA
    sequences; only matrices for \code{del.colgapsonly}.}
}
\details{
  \code{del.gaps} remove all gaps, so the returned sequences may not
  have all the same lengths and are therefore returned in a list.

  \code{del.colgapsonly} removes the columns that contain only gaps
  (useful when a small matrix is extracted from a large alignment).

  The sequences can be either in \code{"DNAbin"} or in another format,
  but the returned object is always of class \code{"DNAbin"}.
}
\value{
  \code{del.gaps} returns a vector (if there is only one input sequence)
  or a list of class \code{"DNAbin"}; \code{del.colgapsonly} returns a
  matrix.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{base.freq}}, \code{\link{GC.content}},
  \code{\link{seg.sites}}, \code{\link{image.DNAbin}}
}
\keyword{univar}

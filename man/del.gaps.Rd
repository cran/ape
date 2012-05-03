\name{del.gaps}
\alias{del.gaps}
\title{
  Delete Alignment Gaps in DNA Sequences
}
\usage{
del.gaps(x)
}
\arguments{
  \item{x}{a matrix, a list, or a vector containing the DNA sequences.}
}
\description{
  This function removes the insertion gaps (\code{"-"}) in a sample of
  DNA sequences.
}
\details{
  The sequences can be either in \code{"DNAbin"} or in character format,
  but the returned object is always of class \code{"DNAbin"}. If
  \code{x} is a vector, then a vector is returned; if it is a list or a
  matrix, then a list is returned.
}
\value{
  A vector (if there is only one input sequence) or a list of class
  \code{"DNAbin"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{base.freq}}, \code{\link{GC.content}},
  \code{\link{seg.sites}}, \code{\link{image.DNAbin}}
}
\keyword{univar}

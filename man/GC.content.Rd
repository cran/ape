\name{GC.content}
\alias{GC.content}
\title{Content in GC from DNA Sequences}
\usage{
GC.content(x)
}
\arguments{
  \item{x}{a vector, a matrix, a data frame, or a list which contains
    the DNA sequences.}
}
\description{
  This function computes the percentage of G+C in a sample of DNA sequences.
}
\details{
  The  percentage of G+C is computed over all sequences in the
  sample. All missing or unknown sites are discarded from the
  computations. The present function actually uses the function
  \code{base.freq}.
}
\value{
  A single numeric value is returned.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{base.freq}}, \code{\link{seg.sites}}
}
\keyword{univar}

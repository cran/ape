\name{base.freq}
\alias{base.freq}
\title{Base frequencies from DNA Sequences}
\usage{
base.freq(x)
}
\arguments{
  \item{x}{a vector, a matrix, a data frame, or a list which contains
    the DNA sequences.}
}
\description{
  This function computes the relative frequencies (i.e. percentages) of
  the four DNA bases (adenine, cytosine, guanine, and thymidine) from a
  sample of sequences.
}
\details{
  The base frequencies are computed over all sequences in the
  sample. All missing or unknown sites are discarded from the
  computations.
}
\value{
  A table of relative frequencies with the names "a", "c", "g", and "t".
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{GC.content}}
}
\keyword{univar}

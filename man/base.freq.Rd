\name{base.freq}
\alias{base.freq}
\title{Base frequencies from DNA Sequences}
\usage{
base.freq(x, freq = FALSE)
}
\arguments{
  \item{x}{a vector, a matrix, or a list which contains the DNA
    sequences.}
  \item{freq}{a logical specifying whether to return the proportions
    (the default) or the absolute frequencies (counts).}
}
\description{
  This function computes the relative frequencies (i.e. proportions) of
  the four DNA bases (adenine, cytosine, guanine, and thymidine) from a
  sample of sequences.
}
\details{
  The base frequencies are computed over all sequences in the
  sample. All missing or unknown sites are discarded from the
  computations.
}
\value{
  A numeric vector storing the relative frequencies with names
  \code{c("a", "c", "g", "t")}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{GC.content}}, \code{\link{seg.sites}},
  \code{\link{nuc.div}}, \code{\link{DNAbin}}
}
\keyword{univar}

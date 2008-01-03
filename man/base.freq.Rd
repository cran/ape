\name{base.freq}
\alias{base.freq}
\title{Base frequencies from DNA Sequences}
\usage{
base.freq(x)
}
\arguments{
  \item{x}{a vector, a matrix, or a list which contains the DNA
    sequences.}
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
  A numeric vector stoting the relative frequencies with names
  \code{c("a", "c", "g", "t")}.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{GC.content}}, \code{\link{seg.sites}},
  \code{\link{nuc.div}}, \code{\link{DNAbin}}
}
\keyword{univar}

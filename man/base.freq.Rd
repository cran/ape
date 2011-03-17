\name{base.freq}
\alias{base.freq}
\alias{Ftab}
\title{Base frequencies from DNA Sequences}
\description{
  \code{base.freq} computes the frequencies (absolute or relative) of
  the four DNA bases (adenine, cytosine, guanine, and thymidine) from a
  sample of sequences.

  \code{Ftab} computes the contingency table with the absolute
  frequencies of the DNA bases from a pair of sequences.
}
\usage{
base.freq(x, freq = FALSE, all = FALSE)
Ftab(x, y = NULL)
}
\arguments{
  \item{x}{a vector, a matrix, or a list which contains the DNA
    sequences.}
  \item{y}{a vector with a single DNA sequence.}
  \item{freq}{a logical specifying whether to return the proportions
    (the default) or the absolute frequencies (counts).}
  \item{all}{a logical; by default only the counts of A, C, G, and T are
    returned. If \code{all = TRUE}, all counts of bases, ambiguous codes,
    missing data, and alignment gaps are returned.}
}
\details{
  The base frequencies are computed over all sequences in the
  sample. All missing or unknown sites are discarded from the
  computations.

  For \code{Ftab}, if the argument \code{y} is given then both \code{x}
  and \code{y} are coerced as vectors and must be of equal length. If
  \code{y} is not given, \code{x} must be a matrix or a list and only
  the two first sequences are used.
}
\value{
  A numeric vector with names \code{c("a", "c", "g", "t")}, or a four by
  four matrix with similar dimnames.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{GC.content}}, \code{\link{seg.sites}},
  \code{\link{nuc.div}}, \code{\link{DNAbin}}
}
\examples{
data(woodmouse)
base.freq(woodmouse)
base.freq(woodmouse, TRUE)
Ftab(woodmouse)
Ftab(woodmouse[1, ], woodmouse[2, ]) # same than above
Ftab(woodmouse[14:15, ]) # between the last two
}
\keyword{univar}
\keyword{manip}

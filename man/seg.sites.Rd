\name{seg.sites}
\alias{seg.sites}
\title{
  Find Segregating Sites in DNA Sequences
}
\usage{
seg.sites(x)
}
\arguments{
  \item{x}{a matrix or a list which contains the DNA sequences.}
}
\description{
  This function gives the indices of segregating (polymorphic) sites in
  a sample of DNA sequences.
}
\details{
  If the sequences are in a list, all the sequences must be of the same
  length. Ambiguous nucleotides are ignored.
}
\value{
  A numeric (integer) vector giving the indices of the segregating
  sites.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{base.freq}}, \code{\link[pegas]{theta.s}},
  \code{\link[pegas]{nuc.div}}
}
\examples{
data(woodmouse)
y <- seg.sites(woodmouse)
y
length(y)
}
\keyword{univar}

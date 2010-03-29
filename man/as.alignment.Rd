\name{as.alignment}
\alias{as.alignment}
\alias{as.DNAbin}
\alias{as.DNAbin.character}
\alias{as.DNAbin.list}
\alias{as.DNAbin.alignment}
\alias{as.character.DNAbin}
\title{Conversion Among DNA Sequence Internal Formats}
\description{
  These functions transform a set of DNA sequences among various
  internal formats.
}
\usage{
as.alignment(x)
as.DNAbin(x, ...)

\method{as.DNAbin}{character}(x, ...)

\method{as.DNAbin}{list}(x, ...)

\method{as.DNAbin}{alignment}(x, ...)

\method{as.character}{DNAbin}(x, ...)
}
\arguments{
  \item{x}{a matrix or a list containing the DNA sequences, or an object
    of class \code{"alignment"}.}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
  For \code{as.alignment}, the sequences given as argument should be
  stored as matrices or lists of single-character strings (the format
  used in \pkg{ape} before version 1.10). The returned object is in the
  format used in the package \pkg{seqinr} to store aligned sequences.

  \code{as.DNAbin} is a generic function with methods so that it works
  with sequences stored into vectors, matrices, or lists.

  \code{as.character} is a generic function: the present method
  converts objects of class \code{"DNAbin"} into the format used
  before \pkg{ape} 1.10 (matrix of single characters, or list of vectors
  of single characters). This function must be used first to convert
  objects of class \code{"DNAbin"} into the class \code{"alignment"}.
}
\value{
  an object of class \code{"alignment"} in the case of
  \code{"as.alignment"}; an object of class \code{"DNAbin"} in the case
  of \code{"as.DNAbin"}; a matrix of mode character or a list containing
  vectors of mode character in the case of \code{"as.character"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{DNAbin}}, \code{\link{read.dna}},
  \code{\link{read.GenBank}}, \code{\link{write.dna}}
}
\examples{
data(woodmouse)
x <- as.character(woodmouse)
x[, 1:20]
str(as.alignment(x))
identical(as.DNAbin(x), woodmouse)
}
\keyword{manip}

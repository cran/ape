\name{write.dna}
\alias{write.dna}
\title{Write DNA Sequences in a File}
\usage{
write.dna(x, file = "", format = "interleaved", append = FALSE)
}
\arguments{
  \item{x}{a list a DNA sequences each made of a single vector of mode
    character where each element is a nucleotide.}
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{file = ""}(the default) then the
    tree is written on the standard output connection.}
  \item{format}{a character string specifying the format of the DNA
    sequences. Two choices are possible: \code{"interleaved"} or
    \code{"sequential"}.}
  \item{append}{a logical, if \code{TRUE} the tree is appended to the file
    without erasing the data possibly existing in the file, otherwise
    the file (if it exists) is overwritten (\code{FALSE} the default).}
}
\description{
  This function writes in a file a list of DNA sequences using either a
  sequential or an interleaved format. The names of the vectors of the
  list are used as taxa names.
}
\value{
  None (invisible `NULL').
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.dna}}, \code{\link{read.GenBank}}
}
\keyword{IO}

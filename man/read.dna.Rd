\name{read.dna}
\alias{read.dna}
\title{Read DNA Sequences in a File}
\usage{
read.dna(file, format = "interleaved", skip = 0,
         nlines = 0, comment.char = "#")
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string.}
  \item{format}{a character string specifying the format of the DNA
    sequences. Two choices are possible: \code{"interleaved"} or
    \code{"sequential"}.}
  \item{skip}{the number of lines of the input file to skip before
    beginning to read data.}
  \item{nlines}{the number of lines to be read (by default the file is
    read untill its end).}
  \item{comment.char}{a single character, the remaining of the line
    after this character is ignored.}
}
\description{
  This function reads DNA sequences in a file, and returns a list of DNA
  sequences with the names of the taxa read in the file as names of the
  vectors of the list.
}
\value{
  A list a DNA sequences each made of a single vector of mode character
  where each element is a nucleotide.
}
\details{
  This function tries to be as flexible as possible with respect to the
  presentation of the sequences. There is no restriction on the lengths
  of the taxa names. 
}
\seealso{
  \code{\link{read.GenBank}}, \code{\link{write.dna}}, \code{\link{dist.dna}}
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\keyword{IO}

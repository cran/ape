\name{read.GenBank}
\alias{read.GenBank}
\title{Read DNA Sequences from GenBank via Internet}
\usage{
read.GenBank(access.nb, seq.names = access.nb,
             species.names = TRUE, as.character = FALSE)
}
\arguments{
  \item{access.nb}{a vector of mode character giving the accession numbers.}
  \item{seq.names}{the names to give to each sequence; by default the
    accession numbers are used.}
  \item{species.names}{a logical indicating whether to attribute the
    species names to the returned object.}
  \item{as.character}{a logical controlling whether to return the
    sequences as an object of class \code{"DNAbin"} (the default).}
}
\description{
  This function connects to the GenBank database, and reads nucleotide
  sequences using accession numbers given as arguments.
}
\value{
  A list of DNA sequences made of vectors of class \code{"DNAbin"}, or
  of single characters (if \code{as.character = "TRUE"}).
}
\details{
  The function uses the site \url{http://www.ncbi.nlm.nih.gov/} from
  where the sequences are downloaded.

  If \code{species.names = TRUE}, the returned list has an attribute
  \code{"species"} containing the names of the species taken from the
  field ``ORGANISM'' in GenBank.
}
\seealso{
  \code{\link{read.dna}}, \code{\link{write.dna}},
  \code{\link{dist.dna}}, \code{\link{DNAbin}}
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\examples{
### This won't work if your computer is not connected
### to the Internet!!!
###
### Get the 8 sequences of tanagers (Ramphocelus)
### as used in Paradis (1997)
ref <- c("U15717", "U15718", "U15719", "U15720",
         "U15721", "U15722", "U15723", "U15724")
### Copy/paste or type the following commands if you
### want to try them.
\dontrun{
Rampho <- read.GenBank(ref)
### get the species names:
attr(Rampho, "species")
### build a matrix with the species names and the accession numbers:
cbind(attr(Rampho, "species"), names(Rampho))
### print the first sequence
### (can be done with `Rampho$U15717' as well)
Rampho[[1]]
### print the first sequence in a cleaner way
cat(Rampho[[1]], "\n", sep = "")
}
}
\keyword{IO}

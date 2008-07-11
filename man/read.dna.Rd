\name{read.dna}
\alias{read.dna}
\title{Read DNA Sequences in a File}
\usage{
read.dna(file, format = "interleaved", skip = 0,
         nlines = 0, comment.char = "#", seq.names = NULL,
         as.character = FALSE)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string.}
  \item{format}{a character string specifying the format of the DNA
    sequences. Four choices are possible: \code{"interleaved"},
    \code{"sequential"}, \code{"clustal"}, or \code{"fasta"}, or any
    unambiguous abbreviation of these.}
  \item{skip}{the number of lines of the input file to skip before
    beginning to read data.}
  \item{nlines}{the number of lines to be read (by default the file is
    read untill its end).}
  \item{comment.char}{a single character, the remaining of the line
    after this character is ignored.}
  \item{seq.names}{the names to give to each sequence; by default the
    names read in the file are used.}
  \item{as.character}{a logical controlling whether to return the
    sequences as an object of class \code{"DNAbin"} (the default).}
}
\description{
  This function reads DNA sequences in a file, and returns a matrix or a
  list of DNA sequences with the names of the taxa read in the file as
  rownames or names, respectively. By default, the sequences are stored
  in binary format, otherwise (if \code{as.character = "TRUE"}) in lower
  case.
}
\details{
  This function follows the interleaved and sequential formats defined
  in PHYLIP (Felsenstein, 1993) but with the original feature than there
  is no restriction on the lengths of the taxa names. For these two
  formats, the first line of the file must contain the dimensions of the
  data (the numbers of taxa and the numbers of nucleotides); the
  sequences are considered as aligned and thus must be of the same
  lengths for all taxa. For the FASTA format, the conventions defined in
  the URL below (see References) are followed; the sequences are taken as
  non-aligned. For all formats, the nucleotides can be arranged in any
  way with blanks and line-breaks inside (with the restriction that the
  first ten nucleotides must be contiguous for the interleaved and
  sequential formats, see below). The names of the sequences are read in
  the file unless the `seq.names' option is used. Particularities for
  each format are detailed below.

  \item{Interleaved:}{the function starts to read the sequences when it
    finds 10 contiguous characters belonging to the ambiguity code of
    the IUPAC (namely A, C, G, T, U, M, R, W, S, Y, K, V, H, D, B, and
    N, upper- or lowercase, so you might run into trouble if you have a
    taxa name with 10 contiguous letters among these!) All characters
    before the sequences are taken as the taxa names after removing the
    leading and trailing spaces (so spaces in a taxa name are
    allowed). It is assumed that the taxa names are not repeated in the
    subsequent blocks of nucleotides.}

  \item{Sequential:}{the same criterion than for the interleaved format
    is used to start reading the sequences and the taxa names; the
    sequences are then read until the number of nucleotides specified in
    the first line of the file is reached. This is repeated for each taxa.}

  \item{Clustal:}{this is the format output by the Clustal programs
    (.aln). It is somehow similar to the interleaved format: the
    differences being that the dimensions of the data are not indicated
    in the file, and the names of the sequences are repeated in each block.}

  \item{FASTA:}{This looks like the sequential format but the taxa names
    (or rather a description of the sequence) are on separate lines
    beginning with a `greater than' character ``>'' (there may be
    leading spaces before this character). These lines are taken as taxa
    names after removing the ``>'' and the possible leading and trailing
    spaces. All the data in the file before the first sequence is ignored.}
}
\value{
  a matrix or a list (if \code{format = "fasta"}) of DNA sequences
  stored in binary format, or of mode character (if \code{as.character =
    "TRUE"}).
}
\references{
  Anonymous. FASTA format description.
  \url{http://www.ncbi.nlm.nih.gov/BLAST/fasta.html}

  Anonymous. IUPAC ambiguity codes.
  \url{http://www.ncbi.nlm.nih.gov/SNP/iupac.html}

  Felsenstein, J. (1993) Phylip (Phylogeny Inference Package) version
  3.5c. Department of Genetics, University of Washington.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}
}
\seealso{
  \code{\link{read.GenBank}}, \code{\link{write.dna}},
  \code{\link{DNAbin}}, \code{\link{dist.dna}}, \code{\link{woodmouse}}
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\examples{
### a small extract from `data(woddmouse)'
cat("3 40",
"No305     NTTCGAAAAACACACCCACTACTAAAANTTATCAGTCACT",
"No304     ATTCGAAAAACACACCCACTACTAAAAATTATCAACCACT",
"No306     ATTCGAAAAACACACCCACTACTAAAAATTATCAATCACT",
file = "exdna.txt", sep = "\n")
ex.dna <- read.dna("exdna.txt", format = "sequential")
str(ex.dna)
ex.dna
### the same data in interleaved format...
cat("3 40",
"No305     NTTCGAAAAA CACACCCACT",
"No304     ATTCGAAAAA CACACCCACT",
"No306     ATTCGAAAAA CACACCCACT",
"          ACTAAAANTT ATCAGTCACT",
"          ACTAAAAATT ATCAACCACT",
"          ACTAAAAATT ATCAATCACT",
file = "exdna.txt", sep = "\n")
ex.dna2 <- read.dna("exdna.txt")
### ... in clustal format...
cat("CLUSTAL (ape) multiple sequence alignment", "",
"No305     NTTCGAAAAACACACCCACTACTAAAANTTATCAGTCACT",
"No304     ATTCGAAAAACACACCCACTACTAAAAATTATCAACCACT",
"No306     ATTCGAAAAACACACCCACTACTAAAAATTATCAATCACT",
"           ************************** ******  ****",
file = "exdna.txt", sep = "\n")
ex.dna3 <- read.dna("exdna.txt", format = "clustal")
### ... and in FASTA format
cat("> No305",
"NTTCGAAAAACACACCCACTACTAAAANTTATCAGTCACT",
"> No304",
"ATTCGAAAAACACACCCACTACTAAAAATTATCAACCACT",
"> No306",
"ATTCGAAAAACACACCCACTACTAAAAATTATCAATCACT",
file = "exdna.txt", sep = "\n")
ex.dna4 <- read.dna("exdna.txt", format = "fasta")
### The first three are the same!
identical(ex.dna, ex.dna2)
identical(ex.dna, ex.dna3)
identical(ex.dna, ex.dna4)
unlink("exdna.txt") # clean-up
}
\keyword{IO}

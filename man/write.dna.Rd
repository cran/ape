\name{write.dna}
\alias{write.dna}
\title{Write DNA Sequences in a File}
\usage{
write.dna(x, file, format = "interleaved", append = FALSE,
          nbcol = 6, colsep = " ", colw = 10, indent = NULL,
          blocksep = 1)
}
\arguments{
  \item{x}{a list a DNA sequences each made of a single vector of mode
    character where each element is a nucleotide.}
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string.}
  \item{format}{a character string specifying the format of the DNA
    sequences. Three choices are possible: \code{"interleaved"},
    \code{"sequential"}, or \code{"fasta"}, or any unambiguous
    abbreviation of these.}
  \item{append}{a logical, if \code{TRUE} the data are appended to the
    file without erasing the data possibly existing in the file,
    otherwise the file (if it exists) is overwritten (\code{FALSE} the
    default).}
  \item{nbcol}{a numeric specifying the number of columns per row (6 by
    default); may be negative implying that the nucleotides are printed
    on a single line.}
  \item{colsep}{a character used to separate the columns (a single
    space by default).}
  \item{colw}{a numeric specifying the number of nucleotides per column
    (10 by default).}
  \item{indent}{a numeric or a character specifying how the blocks of
    nucleotides are indented (see details).}
  \item{blocksep}{a numeric specifying the number of lines between the
    blocks of nucleotides (this has an effect only if `format =
    "interleaved"').}
}
\description{
  This function writes in a file a list of DNA sequences in sequential,
  interleaved, or FASTA format. The names of the vectors of the list are
  used as taxa names.
}
\details{
  The same three formats are supported in the present function than in
  \code{\link{read.dna}}: see its help page and the references below for
  a description of these formats.
  
  If the sequences have no names, then they are given "1", "2", ... as
  names in the file.

  With the interleaved and sequential formats, the sequences must be all
  of the same length; if the taxon names are longer than 10 characters,
  they are truncated and a warning message is issued.
  
  The argument `indent' specifies how the rows of nucleotides are
  indented. In the interleaved and sequential formats, the rows with
  the taxon names are never indented; the subsequent rows are indented
  with 10 spaces by default (i.e. if `indent = NULL)'. In the FASTA
  format, the rows are not indented by default. This default behaviour
  can be modified by specifying a value to `indent': the rows are then
  indented with `indent' (if it is a character) or `indent' spaces (if
  it is a numeric). For example, specifying `indent = "   "' or `indent
  = 3' will have exactly the same effect (use `indent = "\t"' for a
  tabulation).
}
\note{
  Specifying a negative value for `nbcol' (meaning that the nucleotides
  are printed on a single line) gives the same result for the
  interleaved and sequential formats.

  The different options are intended to give flexibility in formatting
  the sequences. For instance, if the sequences are very long it may be
  judicious to remove all the spaces beween columns `(colsep = ""), in
  the margins (indent = 0), and between the blocks (blocksep = 0) to
  produce a smaller file.
}
\value{
  None (invisible `NULL').
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\references{
  Felsenstein, J. (1993) Phylip (Phylogeny Inference Package) version
  3.5c. Department of Genetics, University of Washington.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}

  \url{http://www.ncbi.nlm.nih.gov/SNP/iupac.html}

  \url{http://www.ncbi.nlm.nih.gov/BLAST/fasta.html}
}
\seealso{
  \code{\link{read.dna}}, \code{\link{read.GenBank}}
}
\keyword{IO}

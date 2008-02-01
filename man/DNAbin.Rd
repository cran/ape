\name{DNAbin}
\alias{DNAbin}
\alias{print.DNAbin}
\alias{summary.DNAbin}
\alias{[.DNAbin}
\alias{rbind.DNAbin}
\alias{cbind.DNAbin}
\alias{as.matrix.DNAbin}
\title{Manipulate DNA Sequences in Bit-Level Format}
\description{
  These functions help to manipulate DNA sequences coded in the
  bit-level coding scheme.
}
\usage{
\method{print}{DNAbin}(x, \dots)
\method{summary}{DNAbin}(object, printlen = 6, digits = 3, \dots)
\method{rbind}{DNAbin}(\dots)
\method{cbind}{DNAbin}(\dots, check.names = TRUE)
\method{[}{DNAbin}(x, i, j, drop = TRUE)
\method{as.matrix}{DNAbin}(x, \dots)
}
\arguments{
  \item{x, object}{an object of class \code{"DNAbin"}.}
  \item{\dots}{either further arguments to be passed to or from other
    methods in the case of \code{print}, \code{summary}, and
    \code{as.matrix}, or a series of objects of class \code{"DNAbin"} in
    the case of \code{rbind} and \code{cbind}.}
  \item{printlen}{the number of labels to print (6 by default).}
  \item{digits}{the number of digits to print (3 by default).}
  \item{check.names}{a logical specifying whether to check the rownames
    before binding the columns (see details).}
  \item{i, j}{indices of the rows and/or columns to select or to drop.
    They may be numeric, logical, or character (in the same way than for
    standard R objects).}
  \item{drop}{logical; if \code{TRUE} (the default), the returned object
    is of the lowest possible dimension.}
}
\details{
  These are all `methods' of generic functions which are here applied to
  DNA sequences stored as objects of class \code{"DNAbin"}. They are
  used in the same way than the standard R functions to manipulate
  vectors, matrices, and lists. Additionally, the operators \code{[[}
  and \code{$} may be used to extract a vector from a list.

  These functions are provided to manipulate easily DNA sequences coded
  with the bit-level coding scheme. The latter allows much faster
  comparisons of sequences, as well as storing them in less memory
  compared to the format used before \pkg{ape} 1.10.

  For \code{cbind}, if \code{"check.names = TRUE"}, the rownames of each
  matrix are checked, and the rows are reordered if necessary. If the
  rownames differ among matrices, an error occurs. If
  \code{"check.names = FALSE"}, the matrices are simply binded and the
  rownames of the first matrix are used.

  \code{as.matrix} may be used to convert DNA sequences (of the same
  length) stored in a list into a matrix while keeping the names and the
  class.
}
\value{
  an object of class \code{"DNAbin"} in the case of \code{rbind},
  \code{cbind}, and \code{[}.
}
\references{
  Paradis, E. (2007) A Bit-Level Coding Scheme for Nucleotides.
  \url{http://pbil.univ-lyon1.fr/R/ape/misc/BitLevelCodingScheme_20April2007.pdf}
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{as.DNAbin}}, \code{\link{read.dna}},
  \code{\link{read.GenBank}}, \code{\link{write.dna}}

  The corresponding generic functions are documented in the package
  \pkg{base}.
}
\examples{
data(woodmouse)
woodmouse
summary(woodmouse)
summary(woodmouse, 15, 6)
summary(woodmouse[1:5, 1:300], 15, 6)
### Just to show how distances could be influenced by sampling:
dist.dna(woodmouse[1:2, ])
dist.dna(woodmouse[1:3, ])
}
\keyword{manip}

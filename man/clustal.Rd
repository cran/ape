\name{clustal}
\alias{clustal}
\alias{muscle}
\alias{tcoffee}
\title{Multiple Sequence Alignment with External Applications}
\description{
  These functions call their respective program from \R to align a set of
  nucleotide sequences of class \code{"DNAbin"}.
}
\usage{
clustal(x, pw.gapopen = 10, pw.gapext = 0.1,
        gapopen = 10, gapext = 0.2, exec = NULL,
        MoreArgs = "", quiet = TRUE)
muscle(x, exec = "muscle", MoreArgs = "", quiet = TRUE)
tcoffee(x, exec = "t_coffee", MoreArgs = "", quiet = TRUE)
}
\arguments{
  \item{x}{an object of class \code{"DNAbin"}.}
  \item{pw.gapopen, pw.gapext}{gap opening and gap extension penalties
    used by Clustal during pairwise alignments.}
  \item{gapopen, gapext}{idem for global alignment.}
  \item{exec}{a character string giving the name of the program, with
    its path if necessary. \code{clustal} tries to guess it depending on
    the operating system (see details).}
  \item{MoreArgs}{a character string giving additional options.}
  \item{quiet}{a logical: the default is to not print on \R's console the
    messages from the external program.}
}
\details{
  \code{clustal} tries to guess the name of the executable program
  depending on the operating system. Specifically, the followings are
  used: ``clustalw'' under Linux, ``clustalw2'' under MacOS, or
  ``C:/Program Files/ClustalW2/clustalw2'' under Windows.

  The calculations are done in a temporary directory which is deleted
  when \R is quit. So it is possible to find the files created by the
  last call in the directory printed by \code{tempdir()}.

  When called without arguments (i.e., \code{clustal()}, \dots), the
  function prints the options of the program which may be passed to
  \code{MoreArgs}.
}
\value{
  an object of class \code{"DNAbin"} with the aligned sequences.
}
\references{
  Chenna, R., Sugawara, H., Koike, T., Lopez, R., Gibson, T. J.,
  Higgins, D. G. and Thompson, J. D. (2003) Multiple sequence alignment
  with the Clustal series of programs. \emph{Nucleic Acids Research}
  \bold{31}, 3497--3500.
  \url{http://www.clustal.org/}

  Edgar, R. C. (2004) MUSCLE: Multiple sequence alignment with high
  accuracy and high throughput. \emph{Nucleic Acids Research},
  \bold{32}, 1792--1797.
  \url{http://www.drive5.com/muscle/muscle_userguide3.8.html}

  Notredame, C., Higgins, D. and Heringa, J. (2000) T-Coffee: A novel
  method for multiple sequence alignments. \emph{Journal of Molecular
  Biology}, \bold{302}, 205--217.
  \url{http://www.tcoffee.org/Documentation/t_coffee/t_coffee_technical.htm}
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{image.DNAbin}}, \code{\link{del.gaps}}

  The package \pkg{phyloch} which has similar functions for the MAFFT
  and Prank.
}
\examples{
\dontrun{
### display the options:
clustal()
muscle()
tcoffee()

data(woodmouse)
### open gaps more easily:
clustal(woodmouse, pw.gapopen = 1, pw.gapext = 1)
### T-Coffee requires negative values (quite slow; muscle is much faster):
tcoffee(woodmouse,  MoreArgs = "-gapopen=-10 -gapext=-2")
}}
\keyword{manip}

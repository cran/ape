\name{write.tree}
\alias{write.tree}
\title{Write Tree File in Parenthetic Format}
\usage{
write.tree(phy, file = "", append = FALSE,
           multi.line = TRUE, digits = 10)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{file = ""} (the default) then the
    tree is written on the standard output connection (i.e. the console).}
  \item{append}{a logical, if \code{TRUE} the tree is appended to the file
    without erasing the data possibly existing in the file, otherwise
    the file (if it exists) is overwritten (\code{FALSE} the default).}
  \item{multi.line}{a logical, if \code{TRUE} (the default) line breaks
    are introduced in the tree in parenthetic form.}
  \item{digits}{a numeric giving the number of digits used for printing
    branch lengths.}
}
\description{
  This function writes in a file a tree in parenthetic format using the
  Newick (also known as New Hampshire) format.
}
\value{
  a vector of mode character if \code{file = ""}, none (invisible
  `NULL') otherwise.
}
\details{
  The node labels and the root edge length, if available, are written in
  the file.
}
\references{
  Felsenstein, J. The Newick tree format.
  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}

  Olsen, G. Interpretation of the "Newick's 8:45" tree format standard.
  \url{http://evolution.genetics.washington.edu/phylip/newick_doc.html}
}

\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link{read.nexus}},
  \code{\link{write.nexus}}
}
\keyword{manip}
\keyword{IO}

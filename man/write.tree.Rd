\name{write.tree}
\alias{write.tree}
\title{Write Tree File in Parenthetic Format}
\usage{
write.tree(phy, file = "", append = FALSE, format = "Newick")
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{file = ""}(the default) then the
    tree is written on the standard output connection.}
  \item{append}{a logical, if \code{TRUE} the tree is appended to the file
    without erasing the data possibly existing in the file, otherwise
    the file (if it exists) is overwritten (\code{FALSE} the default).}
  \item{format}{a character giving the format of the tree file; by
    default \code{"Newick"}. Currently only \code{"Newick"} is available.}
}
\description{
  This function writes in a file a tree in parenthetic format using the
  Newick or New Hampshire format.
}
\value{
  None (invisible `NULL').
}
\details{
  The node labels and the root edge length, if available, are written in
  the file.
}
\references{
  \url{http://evolution.genetics.washington.edu/phylip/newick_doc.html}

  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
}

\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link{read.nexus}},
  \code{\link{write.nexus}}
}
\keyword{manip}
\keyword{IO}

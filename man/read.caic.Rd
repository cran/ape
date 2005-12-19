\name{read.caic}
\alias{read.caic}
\title{Read Tree File in CAIC Format}
\description{
This function reads one tree from a CAIC file.
A second file containing branch lengths values may also be passed (experimental).
}
\usage{
read.caic(file, brlen = NULL, skip = 0, comment.char = "#", ...)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character, or a double-quoted string.}
  \item{brlen}{a file name for the branch lengths file.}
  \item{skip}{the number of lines of the input file to skip before beginning to read data (this is passed directly to scan()).}
  \item{comment.char}{a single character, the remaining of the line after this character is ignored (this is passed directly to scan()).}
  \item{\dots}{Further arguments to be passed to scan().}
}
\details{
  Read a tree from a file in the format used by the CAIC and MacroCAIc program.
}
\value{
an object of class "phylo" with the following components:
\item{edge}{a two-column matrix of mode character where each row represents an edge of the tree; the nodes and the tips are symbolized with numbers (these numbers are not treated as numeric, hence the mode character); the nodes are represented with negative numbers (the root being "-1"), and the tips are represented with positive numbers. For each row, the first column gives the ancestor. This representation allows an easy manipulation of the tree, particularly if it is rooted.}
\item{edge.length}{a numeric vector giving the lengths of the branches given by edge.}
\item{tip.label}{a vector of mode character giving the names of the tips; the order of the names in this vector corresponds to the (positive) number in edge.}
\item{node.label}{(optional) a vector of mode character giving the names of the nodes (set to NULL if not available in the file).}
\item{root.edge}{(optional) a numeric value giving the length of the branch at the root is it exists (NULL otherwise).}
}
\references{
  Purvis, A. and Rambaut, A. (1995) Comparative analysis by independent
  contrasts (CAIC): an Apple Macintosh application for analysing
  comparative data. \emph{CABIOS}, \bold{11} :241--251.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\section{Warning }{The branch length support is still experimental and was not fully tested.}
\seealso{ \code{\link{read.tree}}, \code{\link{read.nexus}} }
\examples{
### The same example than in read.tree, without branch lengths.
### An extract from Sibley and Ahlquist (1990)
cat("AAA","Strix_aluco","AAB","Asio_otus",
   "AB","Athene_noctua","B","Tyto_alba",
   file = "ex.tre", sep = "\n")
tree.owls <- read.caic("ex.tre")
plot(tree.owls)
tree.owls
unlink("ex.tre") # delete the file "ex.tre"
}
\keyword{hplot}

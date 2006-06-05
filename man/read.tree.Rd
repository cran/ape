\name{read.tree}
\alias{read.tree}
\title{Read Tree File in Parenthetic Format}
\usage{
read.tree(file = "", format = "Newick", rooted = TRUE, text = NULL,
          tree.names = NULL, skip = 0, comment.char = "#", ...)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{file = ""} (the default) then the
    tree is input on the keyboard, the entry being terminated with a
    blank line.}
  \item{format}{a character giving the format of the tree file; by
    default \code{"Newick"}. Currently only the Newick format is
    available (see \code{\link{read.nexus}} to read trees in Nexus
    files).}
  \item{rooted}{a logical indicating whether the tree is rooted or not;
    defaults to \code{TRUE}. This is currently not implemented, meaning
    that all trees are read as rooted.}
  \item{text}{alternatively, the name of a variable of mode character
    which contains the tree(s) in parenthetic format. By default, this
    is ignored (set to \code{NULL}, meaning that the tree is read in a
    file); if \code{text} is not \code{NULL}, then the argument
    \code{file} is ignored.}
  \item{tree.names}{if there are several trees to be read, a vector of
    mode character that gives names to the individual trees; if
    \code{NULL} (the default), the trees are named \code{"tree1"},
    \code{"tree2"}, ...}
  \item{skip}{the number of lines of the input file to skip before
    beginning to read data (this is passed directly to\code{ scan()}).}
  \item{comment.char}{a single character, the remaining of the line
    after this character is ignored (this is passed directly to
    \code{scan()}).}
  \item{...}{Further arguments to be passed to \code{scan()}.}
}
\description{
  This function reads a file which contains one or several trees in
  parenthetic format known as the Newick or New Hampshire format.
}
\details{
  The default option for \code{file} allows to type directly the tree on
  the keyboard (or possibly to copy from an editor and paste in R's
  console) with, e.g., \code{mytree <- read.tree()}.

  `read.tree' tries to represent correctly trees with a badly
  represented root edge (i.e. with an extra pair of parentheses). For
  instance, the tree "((A:1,B:1):10);" will be read like "(A:1,B:1):10;"
  but a warning message will be issued in the former case as this is
  apparently not a valid Newick format. If there are two root edges
  (e.g., "(((A:1,B:1):10):10);"), then the tree is not read and an error
  message is issued.
}
\value{
  an object of class \code{"phylo"} with the following components:
  \item{edge}{a two-column matrix of mode character where each row
    represents an edge of the tree; the nodes and the tips are
    symbolized with numbers (these numbers are not treated as numeric,
    hence the mode character); the nodes are represented with negative
    numbers (the root being \code{"-1"}), and the tips are represented with
    positive numbers. For each row, the first column gives the
    ancestor. This representation allows an easy manipulation of the
    tree, particularly if it is rooted.}
  \item{edge.length}{(optional) a numeric vector giving the lengths of the
    branches given by \code{edge}.}
  \item{tip.label}{a vector of mode character giving the names of the
    tips; the order of the names in this vector corresponds to the
    (positive) number in \code{edge}.}
  \item{node.label}{(optional) a vector of mode character giving the
    names of the nodes (set to \code{NULL} if not available in the file).}
  \item{root.edge}{(optional) a numeric value giving the length of the
    branch at the root if it exists (\code{NULL} otherwise).}

  If several trees are read in the file, the returned object is of class
  \code{c("multi.tree", "phylo")}, and is a list of objects of class
  \code{"phylo"}.
}
\references{
  Felsenstein, J. The Newick tree format.
  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}

  Olsen, G. Interpretation of the "Newick's 8:45" tree format standard.
  \url{http://evolution.genetics.washington.edu/phylip/newick_doc.html}
}

\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{write.tree}}, \code{\link{read.nexus}},
  \code{\link{write.nexus}}, \code{\link[base]{scan}} for the basic R
  function to read data in a file
}
\examples{
### An extract from Sibley and Ahlquist (1990)
cat("(((Strix_aluco:4.2,Asio_otus:4.2):3.1,",
   "Athene_noctua:7.3):6.3,Tyto_alba:13.5);",
   file = "ex.tre", sep = "\n")
tree.owls <- read.tree("ex.tre")
str(tree.owls)
tree.owls
unlink("ex.tre") # delete the file "ex.tre"
### Only the first three species using the option `text'
TREE <- "((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);"
TREE
tree.owls.bis <- read.tree(text = TREE)
str(tree.owls.bis)
tree.owls.bis
}
\keyword{manip}
\keyword{IO}

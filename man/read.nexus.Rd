\name{read.nexus}
\alias{read.nexus}
\title{Read Tree File in Nexus Format}
\usage{
read.nexus(file, tree.names = NULL)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string.}
  \item{tree.names}{if there are several trees to be read, a vector of
    mode character that gives names to the individual trees; if
    \code{NULL} (the default), the trees are named \code{"tree1"},
    \code{"tree2"}, ...}
}
\description{
  This function reads one or several trees in a NEXUS file.
}
\details{
  The present implementation tries to follow as much as possible the
  NEXUS standard. Only the block "TREES" is read; the other data can be
  read with other functions (e.g., \code{\link{read.dna}},
  \code{\link[base]{read.table}}, ...). A trace of the original data is
  kept with the attribute \code{"origin"} (see below).

  `read.nexus tries to represent correctly trees with a badly
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
  \item{edge.length}{a numeric vector giving the lengths of the
    branches given by \code{edge}.}
  \item{tip.label}{a vector of mode character giving the names of the
    tips; the order of the names in this vector corresponds to the
    (positive) number in \code{edge}.}
  \item{node.label}{(optional) a vector of mode character giving the
    names of the nodes (set to \code{NULL} if not available in the file).}
  \item{root.edge}{(optional) a numeric value giving the length of the
    branch at the root is it exists (\code{NULL} otherwise).}

  If several trees are read in the file, the returned object is of class
  \code{c("phylo", "multi.tree")}, and is a list of objects of class
  \code{"phylo"}.

  An attribute \code{"origin"} is further given to the returned object
  which gives the name of the source file (with its path). This is used
  to write a tree in a NEXUS file where all the original data must be
  written (not only the tree) in accordance to the specifications of
  Maddison et al. (1997).
}
\references{
  Maddison, D. R., Swofford, D. L. and Maddison, W. P. (1997) NEXUS: an
  extensible file format for systemaric information. \emph{Systematic
    Biology}, \bold{46}, 590--621.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link{write.nexus}},
  \code{\link{write.tree}}
}
\keyword{manip}
\keyword{IO}

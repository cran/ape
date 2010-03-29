\name{print.phylo}
\alias{print.phylo}
\alias{print.multiPhylo}
\alias{[.multiPhylo}
\alias{[[.multiPhylo}
\alias{$.multiPhylo}
\alias{str.multiPhylo}
\title{Compact Display of a Phylogeny}
\usage{
\method{print}{phylo}(x, printlen = 6 ,...)
\method{print}{multiPhylo}(x, details = FALSE ,...)
\method{[}{multiPhylo}(x, i)
\method{[[}{multiPhylo}(x, i)
\method{$}{multiPhylo}(x, name)
\method{str}{multiPhylo}(object, ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{object}{an object of class \code{"multiPhylo"}.}
  \item{printlen}{the number of labels to print (6 by default).}
  \item{details}{a logical indicating whether to print information on
    all trees.}
  \item{i}{indices of the tree(s) to select from a list; this may be a
    vector of integers, logicals, or names.}
  \item{name}{a character string specifying the tree to be extracted.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  These functions prints a compact summary of a phylogeny, or a list of
  phylogenies, on the console.

  The operators \code{[}, \code{[[}, and \code{$} propagate the class
  correctly.
}
\value{
  An object of class \code{"phylo"} (\code{[[}, \code{$}) or of class
  \code{"multiPhylo"} (\code{[}), or NULL.
}
\author{Ben Bolker \email{bolker@zoo.ufl.edu} and Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}}, \code{\link{summary.phylo}},
  \code{\link[base]{print}} for the generic R function
}
\keyword{manip}

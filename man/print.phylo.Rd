\name{print.phylo}
\alias{print.phylo}
\alias{print.multiPhylo}
\alias{[.multiPhylo}
\title{Compact Display of a Phylogeny}
\usage{
\method{print}{phylo}(x, printlen = 6 ,...)
\method{print}{multiPhylo}(x, details = FALSE ,...)
\method{[}{multiPhylo}(x, i)
}
\arguments{
  \item{x}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{printlen}{the number of labels to print (6 by default).}
  \item{details}{a logical indicating whether to print information on
    all trees.}
  \item{i}{indices of the trees to select from a list; this may be a
    vector of integers, logicals, or names.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  These functions prints a compact summary of a phylogeny, or a list of,
  on the console.
}
\value{
  An object of class \code{"multiPhylo"} or NULL.
}
\author{Ben Bolker \email{bolker@zoo.ufl.edu} and Emmanuel Paradis
  \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link{summary.phylo}},
  \code{\link[base]{print}} for the generic R function
}
\keyword{manip}

\name{print.phylo}
\alias{print.phylo}
\alias{print.multi.tree}
\title{Compact Display of a Phylogeny}
\usage{
\method{print}{phylo}(x, printlen = 6 ,...)
\method{print}{multi.tree}(x, details = FALSE ,...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{printlen}{the number of labels to print (6 by default).}
  \item{details}{a logical indicating whether to print information on
    all trees.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  These functions prints a compact summary of a phylogeny, or a list of,
  on the console.
}
\value{
  A NULL value is returned, the results are simply printed.
}
\author{Ben Bolker \email{bolker@zoo.ufl.edu} and Emmanuel Paradis
  \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link{summary.phylo}},
  \code{\link[base]{print}} for the generic R function
}
\keyword{manip}

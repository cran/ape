\name{mrca}
\alias{mrca}
\title{Find Most Recent Common Ancestors Between Pairs}
\description{
  This function returns for each pair of tips (and nodes) its most
  recent common ancestor (MRCA).
}
\usage{
mrca(phy, full = FALSE, as.numeric = FALSE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{full}{a logical indicating whether to return the MRCAs among
    all tips and nodes (if \code{TRUE}); the default is to return only
    the MRCAs among tips.}
  \item{as.numeric}{a logical indicating whether to return the matrix of
    tip and node numbers as numeric.}
}
\details{
  The diagonal is set to the number of the tips (and nodes if \code{full
    = TRUE}).

  If \code{full = FALSE}, the colnames and rownames are set with the tip
  labels of the tree; otherwise the numbers are given as names.
}
\value{
  a matrix of mode character (by default) or numeric if \code{as.numeric
  = TRUE}.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\keyword{manip}

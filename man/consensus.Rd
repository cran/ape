\name{consensus}
\alias{consensus}
\title{Concensus Trees}
\usage{
consensus(..., p = 1)
}
\arguments{
  \item{...}{either (i) a single object of class \code{"phylo"}, (ii) a
    series of such objects separated by commas, or (iii) a list
    containing such objects.}
  \item{p}{a numeric value between 0.5 and 1 giving the proportion for a
    clade to be represented in the consensus tree.}
}
\description{
  Given a series of trees, this function returns the consensus tree. By
  default, the strict-consensus tree is computed. To get the
  majority-rule consensus tree, use \code{p = 0.5}. Any value between
  0.5 and 1 can be used.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{prop.part}}, \code{\link{dist.topo}}
}
\keyword{manip}

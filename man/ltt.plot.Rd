\name{ltt.plot}
\alias{ltt.plot}
\title{Lineages Through Time Plot}
\usage{
ltt.plot(phy, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{...}{further arguments to be passed to \code{plot()} (see
    \code{Details:} on how to transform the axes).}
}
\description{
  This plots, on the current graphical device, the minimum numbers of
  lineages through time from a phylogenetic tree.
}
\details{
  The axes can be log-transformed by adding one of the following
  options: \code{log = "x"}, \code{log = "y"}, \code{log = "xy"}, to
  transform the x-axis, the y-axis, or both, respectively.
}
\references{
  Harvey, P. H., May, R. M. and Nee, S. (1994) Phylogenies without
  fossils. \emph{Evolution}, \bold{48}, 523--529.

  Nee, S., Holmes, E. C., Rambaut, A. and Harvey, P. H. (1995) Inferring
  population history from molecular phylogenies. \emph{Philosophical
    Transactions of the Royal Society of London. Series B. Biological
    Sciences}, \bold{349}, 25--31.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link[base]{plot}} for the basic plotting function in R,
  \code{\link{skyline.plot}}, \code{\link{branching.times}}
}
\keyword{hplot}

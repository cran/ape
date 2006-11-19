\name{as.phylo.formula}
\alias{as.phylo.formula}
\title{ Conversion from Taxonomy Variables to Phylogenetic Trees }
\description{
  The function \code{as.phylo.formula} (short form \code{as.phylo})
  builds a phylogenetic tree (an object of class \code{phylo}) from
  a set of nested taxonomic variables.
}
\usage{
\method{as.phylo}{formula}(x, data = parent.frame(), ...)
}
\arguments{
  \item{x}{ a right-side formula describing the taxonomic relationship: \code{~C1/C2/.../Cn}. }
  \item{data}{ the data.frame where to look for the variables (default to environment). }
  \item{\dots}{ further arguments to be passed from other methods.}
}
\details{
  Taxonomic variables must be nested and passed in the correct order: the higher clade must be on the left of the formula, for instance \code{~Order/Family/Genus/Species}.
	In most cases, the resulting tree will be unresolved and contain polytomies.
}
\value{
  An object of class \code{phylo}.
}
\author{Julien Dutheil \email{Julien.Dutheil@univ-montp2.fr}}
\seealso{ \code{\link{as.phylo}}, \code{\link{read.tree}} for a description of \code{phylo} objects, \code{\link{multi2di}} }
\examples{
data(carnivora)
plot(as.phylo(~SuperFamily/Family/Genus/Species, data=carnivora))
}
\keyword{manip}

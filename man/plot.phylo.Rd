\name{plot.phylo}
\alias{plot.phylo}
\title{Plot Phylogenies}
\usage{
plot.phylo(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{...}{further arguments to be passed to \code{plot()}.}
}
\description{
  This plots a phylogenetic tree on the current graphical device.
}
\note{
  For the moment, trees are drawn only horizontally with tip label in
  italics. Many other functionalities will be added soon.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link[base]{plot}} for the basic plotting function in R,
  \code{\link{read.tree}}
}
\examples{
### An extract from Sibley and Ahlquist (1990)
cat("(((Strix aluco:4.2,Asio otus:4.2):3.1,",
   "Athene noctua:7.3):6.3,Tyto alba:13.5);",
   file = "ex.tre", sep = "\n")
tree.owls <- read.tree("ex.tre")
plot(tree.owls)
unlink("ex.tre") # delete the file "ex.tre"
}
\keyword{hplot}

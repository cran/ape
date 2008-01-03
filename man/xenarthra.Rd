\name{xenarthra}
\alias{xenarthra}
\title{Molecular Phylogeny of Living Xenarthrans}
\description{
  This phylogeny was inferred by maximum likelihood analysis of the
  nuclear gene BRCA1 (breast cancer susceptibility, 2788 sites)
  sequences for 47 placental and 3 marsupial taxa.
}
\usage{
data(xenarthra)
}
\format{
  The data are stored as an object of class \code{"phylo"} which
  structure is described in the help page of the function
  \code{\link{read.tree}}.
}
\source{
  Delsuc, F., Scally, M., Madsen, O., Stanhope, M. J., de Jong, W. W.,
  Catzeflis, F. M., Springer, M. S. and Douzery, E. J. P. (2002)
  Molecular phylogeny of living xenarthrans and the impact of character
  and taxon sampling on the placental tree rooting. \emph{Molecular
    Biology and Evolution}, \bold{19}, 1656--1671.
}
\seealso{
  \code{\link{read.tree}}
}
\examples{
data(xenarthra)
plot(xenarthra)
### remove the margins...
plot(xenarthra, no.margin = TRUE)
### ... and use a smaller font size
plot(xenarthra, no.margin = TRUE, cex = 0.8)
}
\keyword{datasets}

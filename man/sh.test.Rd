\name{sh.test}
\alias{sh.test}
\title{Shimodaira-Hasegawa Test}
\usage{
sh.test(..., x, model = DNAmodel(), B = 100)
}
\arguments{
  \item{...}{either a series of objects of class \code{"phylo"}
    separated by commas, or a list containing such objects.}
  \item{x}{a list, a matrix, or a data frame containing the (aligned)
    DNA sequences.}
  \item{model}{the model to be fitted to each tree (as an object of
    \code{"DNAmodel"}).}
  \item{B}{the number of bootstrap replicates.}
}
\description{
  This function computes the Shimodaira--Hasegawa test for a set of
  trees.
}
\details{
  The present implementation follows the original formulation of
  Shimodaira and Hasegawa (1999). A difference is that the bootstrap
  resampling is done on the original sequence data rather than the RELL
  method as sugested by Shimodaira and Hasegawa.
}
\value{
  a numeric vector with the P-value associated with each tree given in
  \code{...}.
}
\references{
  Shimodaira, H. and Hasegawa, M. (1999) Multiple comparisons of
  log-likelihoods with applications to phylogenetic
  inference. \emph{Molecular Biology and Evolution}, \bold{16},
  1114--1116.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{mlphylo}}, \code{\link{DNAmodel}}
}
\examples{
data(woodmouse)
t1 <- nj(dist.dna(woodmouse))
t2 <- rtree(15, tip.label = t1$tip.label)
t3 <- rtree(15, tip.label = t1$tip.label)
### Are the NJ tree and two random tress significantly different?
sh.test(t1, t2, t3, x = woodmouse, B = 100)
}
\keyword{models}

\name{vcv.phylo}
\alias{vcv.phylo}
\title{Phylogenetic Variance-covariance or Correlation Matrix}
\usage{
vcv.phylo(phy, model = "Brownian", cor = FALSE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{model}{a character giving the model used to compute the
    variances and covariances of the phynotype; by default
    \code{"Brownian"}. Currently only Brownian is available.}
  \item{cor}{a logical indicating whether the correlation matrix should
    be returned (\code{TRUE}); by default the variance-covariance matrix
    is returned (\code{FALSE}).}
}
\description{
  This function computes the expected variances and covariances of a
  continuous phenotype assuming it evolves under a given model
  (currently only the model of Brownian motion is available).
}
\value{
  a numeric matrix with the names of the tips (as given by the tip.label
  of the argument `phy') as colnames and rownames.
}
\references{
  Garland, T. Jr. and Ives, A. R. (2000) Using the past to predict the
  present: confidence intervals for regression equations in phylogenetic
  comparative methods. \emph{American Naturalist}, \bold{155}, 346--364.
}

\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}} to read tree files in Newick format
}
\keyword{manip}
\keyword{multivariate}

\name{dist.taxo}
\alias{dist.taxo}
\title{Pairwise Distances from a Taxonomic Level}
\usage{
  dist.taxo(x)
}
\arguments{
  \item{x}{a vector describing a taxonomic level}
}
\description{
  This function computes the pairwise distances between the pairs of
  data in an alpha numeric vector.
  This vector contains taxonomic data, either as a full description 
  (for instance 'animal' and 'vegetal') or as a code (for instance '0' and '1').
  This function is useful for drawing correlograms.
}
\details{
  The matrix \eqn{m_{i,j}}{mij} will be:
  \enumerate{
    \item 0 if \eqn{i=j}{i=j};
    \item 0 if \eqn{i \neq j}{i!=j} and \eqn{x_i \neq x_j}{x[i]!=x[j]};
    \item 1 if \eqn{i \neq j}{i!=j} and \eqn{x_i = x_j}{x[i]==x[j]}.
  }
}
\value{
  a square matrix of size \code{length(x)}.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
	\code{\link{Moran.I}} for example of use
}
\keyword{regression}

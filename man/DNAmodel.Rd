\name{DNAmodel}
\alias{DNAmodel}
\title{Defines Models of DNA Evolution}
\usage{
DNAmodel(model = "K80", partition = 1,
         ncat.isv = 1, invar = FALSE,
         equal.isv = TRUE, equal.invar = 1)
}
\arguments{
  \item{model}{a vector of mode character giving the substition model.}
  \item{partition}{a vector of integers defining the partitions for the
    substition models (eventually recycled).}
  \item{ncat.isv}{the number of categories in each partition.}
  \item{invar}{a logical value specifying whether there are invariants.}
  \item{equal.isv}{a logical value specifying whether the `alpha'
    parameter is the same in all partitions; has no effet if \code{ncat
      = 1} or if \code{partition = 1}.}
  \item{equal.invar}{similar to the argument above but for the
    proportion of invariants.}
}
\description{
  This function defines a model of evolution for a set of DNA sequences
  with possible partitions.
}
\details{
  \code{partition} is recycled along the sequence: thus by default there
  is a single partition. For instance, to partition a sequence of 1000
  sites into two partitions of equal length, one will use
  \code{partition = c(rep(1, 500), rep(2, 500))}. The partitions must be
  numbered with a series of integers (1, 2, 3, ...). To partition the
  codon positions, one could do \code{partition = c(1, 1, 2)}.

  The substition models are the same in all partitions. Branch lengths
  are the same in all partitions up to a multiplying coefficient (the
  contrast parameter, denoted 'xi').

  The substitution models must be among the followings: \code{"JC69"}
  \code{"K80"}, \code{"F81"}, \code{"F84"}, \code{"HKY85"},
  \code{"T92"}, \code{"TN93"}, and \code{"GTR"}. These models (except
  HKY85 and GTR) are described in the help page of
  \code{\link{dist.dna}}.

  Inter-sites variation in substitution rates (ISV) is allowed by
  specifying \code{ncat.isv} greater than one.
}
\note{
  The result of this function is not intended to be used by the user,
  but rather to be passed to \code{\link{mlphylo}}.
}
\value{
  an object of class \code{"DNAmodel"} with components defined by the
  arguments of the function call.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{mlphylo}}, \code{\link{dist.dna}}
}
\examples{
### the K80 model:
mod <- DNAmodel()
### the simplest substitution model:
mod <- DNAmodel("JC69")
### the classical GTR + G4 + I:
mod <- DNAmodel("GTR", ncat.isv = 4, invar = TRUE)
### codon-partitioning (with K80):
mod <- DNAmodel(partition = c(1, 1, 2))
### the same but adding inter-sites variation (the alpha parameter
### is the same for both partitions):
mod <- DNAmodel(partition = c(1, 1, 2), ncat.isv = 4)
### ... and with different `alpha' for each partition:
mod <- DNAmodel(partition = c(1, 1, 2), ncat.isv = 4, equal.isv = FALSE)
}
\keyword{models}

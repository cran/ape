\name{DNAmodel}
\alias{DNAmodel}
\title{Defines Models of DNA Evolution}
\usage{
DNAmodel(model = "K80", part.model = 1,
         ncat = 1, part.gamma = 1,
         invar = FALSE, part.invar = 1)
}
\arguments{
  \item{model}{a vector of mode character giving the substition model
    for each partition (must be of the same length than
    \code{unique(part.model)}).}
  \item{part.model}{a vector of integers defining the partitions for the
    substition models (eventually recycled).}
  \item{ncat}{the number of categories in each partition defined by
    \code{part.gamma} (must be of the same length than
    \code{unique(part.gamma)}).}
  \item{part.gamma}{a vector of integers defining the partitions for the
    inter-sites variation in substitution parameters (eventually
    recycled).}
  \item{invar}{a vector of logical values specifying whether there are
    invariants in each partition defined by \code{part.invar} (must be
    of the same length than \code{unique(part.invar)}).}
  \item{part.invar}{a vector of integers defining the partitions for the
    proportion of invariants (eventually recycled).}
}
\description{
  This function defines a model of evolution for a set of DNA sequences
  with possible partitions.
}
\details{
  The six arguments of this function work in pairs: \code{model} defines
  the substitution model for each partition defined by
  \code{part.model}. The latter is recycled along the sequence: thus by
  default there is a single partition. For instance, to partition a
  sequence of 1000 sites into two partitions of equal length, one will
  use \code{part.model = c(rep(1, 500), rep(2, 500))}. The partitions
  must be numbered with a series of integers (1, 2, 3, ...). To
  partition the codon positions, one could do \code{part.model = c(1, 1,
    2)}.

  There must be as many models than partitions, but the former may be
  the same. For instance, \code{model = c("K80", "K80")} means that
  partitions 1 and 2 evolve under Kimura's model but its parameters will
  be estimated separately in each partition. In addition, branch lengths
  are the same in all partitions up to a multiplying coefficient (the
  contrast parameter, denoted 'xi').

  The substitution models must be among the followings: \code{"JC69"}
  \code{"K80"}, \code{"F81"}, \code{"F84"}, \code{"HKY85"},
  \code{"T92"}, \code{"TN93"}, and \code{"GTR"}. These models (except
  HKY85 and GTR) are described in the help page of \code{\link{dist.dna}}.

  The two other pairs of functions work in the same way. This decoupling
  allows one to share some parameters accross partitions (see examples).
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
### the homogeneous K80 model:
mod <- DNAmodel()
### the simplest substitution model:
mod <- DNAmodel("JC69")
### the classical GTR + G4 + I:
mod <- DNAmodel("GTR", ncat = 4, invar = TRUE)
### assuming K80 for the two first codon-positions, and
### JC69 for the third position:
mod <- DNAmodel(c("K80", "JC69"), part.model = c(1, 1, 2))
### the same but adding inter-sites variation (the alpha parameter
### is the same for both partitions):
mod <- DNAmodel(c("K80", "JC69"), part.model = c(1, 1, 2),
                ncat = 4)
### Only the alpha parameter of ISV is different between the
### partitions (K80 is assumed with equal Ts/Tv):
mod <- DNAmodel(ncat = c(4, 4), part.gamma = c(1, 1, 2))
### A much more complete model resulting in three
### overall partitions:
mod <- DNAmodel(c("K80", "K80"), part.model = c(rep(1, 600), rep(2, 200)),
                ncat = c(4, 4), part.gamma = c(rep(1, 200), rep(2, 600)))
}
\keyword{models}

\name{chronopl}
\alias{chronopl}
\title{Molecular Dating With Penalized Likelihood}
\usage{
chronopl(phy, lambda, node.age = 1, node = "root", CV = FALSE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{lambda}{value of the smoothng parameter.}
  \item{node.age}{numeric values specifying the fixed node ages.}
  \item{node}{the numbers of the nodes whose ages are given by
    \code{node.age}; \code{"root"} is a short-cut the number of the
    node.}
  \item{CV}{whether to perform cross-validation.}
}
\description{
  This function estimates the node ages of a tree using semi-parametric
  method based on penalized likelihood (Sanderson 2002). The branch
  lengths of the input tree are interpreted as (mean) numbers of
  substitutions.
}
\details{
  The idea of this method is to use a trade-off between a parametric
  formulation where each branch has its own rate, and a nonparametric
  term where changes in rates are minimized between contiguous
  branches. A smoothing parameter (lambda) controls this trade-off. If
  lambda = 0, then the parametric component dominates and rates vary as
  much as possible among branches, whereas for increasing values of
  lambda, the variation are smoother to tend to a clock-like model (same
  rate for all branches).

  \code{lambda} must be given. The known ages are given in
  \code{node.age}, and the correponding node numbers in \code{node}.
  These two arguments must obviously be of the same length. By default,
  an age of 1 is assumed for the root, and the ages of the other nodes
  are estimated.

  The cross-validation used here is different from the one proposed by
  Sanderson (2002). Here, each tip is dropped successively and the
  analysis is repeated with the reduced tree: the estimated dates for
  the remaining nodes are compared with the estimates from the full
  data. For the \eqn{i}{i}th tip the following is calculated:

  \deqn{\sum_{j=1}^{n-2}{\frac{(t_j - t_j^{-i})^2}{t_j}}}{SUM[j = 1, ..., n-2] (tj - tj[-i])^2/tj},

  where \eqn{t_j}{tj} is the estimated date for the \eqn{j}{j}th node
  with the full phylogeny, \eqn{t_j^{-i}}{tj[-i]} is the estimated date
  for the \eqn{j}{j}th node after removing tip \eqn{i}{i} from the tree,
  and \eqn{n}{n} is the number of tips.
}
\value{
  an object of class \code{"phylo"} with branch lengths as estimated by
  the function. There are two or three further attributes:

  \item{ploglik}{the maximum penalized log-likelihood.}
  \item{rates}{the estimated rates for each branch.}
  \item{D2}{the influence of each observation on overall date
    estimates (if \code{CV = TRUE}).}
}
\references{
  Sanderson, M. J. (2002) Estimating absolute rates of molecular
  evolution and divergence times: a penalized likelihood
  approach. \emph{Molecular Biology and Evolution}, \bold{19},
  101--109.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}}
\seealso{
  \code{\link{chronogram}}, \code{\link{ratogram}},
  \code{\link{NPRS.criterion}}, \code{\link{chronoMPL}}
}
\keyword{models}

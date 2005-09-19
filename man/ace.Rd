\name{ace}
\alias{ace}
\title{Ancestral Character Estimation}
\usage{
ace(x, phy, type = "continuous", method = "ML", CI = TRUE,
    model = if (type == "continuous") "BM" else "ER",
    scaled = TRUE, kappa = 1, corStruct = NULL)
}
\arguments{
  \item{x}{a vector or a factor.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{type}{the type of variable; either \code{"continuous"} or
    \code{"discrete"} (or an abbreviation of these).}
  \item{method}{a character specifying the method used for
    estimation. Three choices are possible: \code{"ML"}, \code{"pic"},
    or \code{"GLS"}.}
  \item{CI}{a logical specifying whether to return the 95\% confidence
    intervals of the ancestral state estimates (for continuous
    characters) or the likelihood of the different states (for discrete
    ones).}
  \item{model}{a character specifying the model (ignored if \code{method
      = "GLS"}), or a numeric matrix if \code{type = "discrete"} (see
    details).}
  \item{scaled}{a logical specifying whether to scale the contrast
    estimate (used only if \code{method = "pic"}).}
  \item{kappa}{a positive value giving the exponent transformation of
    the branch lengths (see details).}
  \item{corStruct}{if \code{method = "GLS"}, specifies the correlation
    structure to be used (this also gives the assumed model).}
}
\description{
  This function estimates ancestral character states, and the associated
  uncertainty, for continuous and discrete characters.
}
\details{
  If \code{type = "continuous"}, the default model is Brownian motion
  where characters evolve randomly following a random walk. This model
  can be fitted by maximum likelihood (the default, Schluter et
  al. 1997), least squares (\code{method = "pic"}, Felsenstein 1985), or
  generalized least squares (\code{method = "GLS"}, Martins and Hansen
  1997). In the latter case, the specification of \code{phy} and
  \code{model} are actually ignored: it is instead given through a
  correlation structure with the option \code{corStruct}.

  For discrete characters (\code{type = "discrete"}), only maximum
  likelihood estimation is available (Pagel 1994). The model is
  specified through a numeric matrix with integer values taken as
  indices of the parameters. The numbers of rows and of columns of this
  matrix must be equal, and are taken to give the number of states of
  the character. For instance, \code{matrix(c(0, 1, 1, 0), 2)} will
  represent a model with two character states and equal rates of
  transition, \code{matrix(c(0, 1, 2, 0), 2)} a model with unequal
  rates, \code{matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3)} a model with
  three states and equal rates of transition (the diagonal is always
  ignored). There are short-cuts to specify these models: \code{"ER"} is
  an equal-rates model (e.g., the first and third examples above),
  \code{"ARD"} is an all-rates-different model (the second example), and
  \code{"SYM"} is a symmetrical model (e.g., \code{matrix(c(0, 1, 2, 1,
    0, 3, 2, 3, 0), 3)}). If a short-cut is used, the number of states
  if determined from the data.
}
\note{
  The current version (as in ape 1.7) is somehow experimental, and the
  code has not been thoroughly tested. Any feedback will be welcome!
}
\value{
  a list with the following elements:

  \item{ace}{if \code{type = "continuous"}, the estimates of the
    ancestral character values.}
  \item{CI95}{if \code{type = "continuous"}, the estimated 95\%
    confidence intervals.}
  \item{sigma2}{if \code{type = "continuous"}, \code{model = "BM"}, and
    \code{method = "ML"}, the maximum likelihood estimate of the
    Brownian parameter.}
  \item{rates}{if \code{type = "discrete"}, the maximum likelihood
    estimates of the transition rates.}
  \item{loglik}{if \code{method = "ML"}, the maximum log-likelihood.}
  \item{lik.anc}{if \code{type = "discrete"}, the scaled likelihoods of
    each ancestral state.}
  \item{call}{the function call.}
}
\references{
  Felsenstein, J. (1985) Phylogenies and the comparative
  method. \emph{American Naturalist}, \bold{125}, 1--15.

  Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the
  comparative method: a general approach to incorporating phylogenetic
  information into the analysis of interspecific data. \emph{American
    Naturalist}, \bold{149}, 646--667.

  Pagel, M (1994) Detecting correlated evolution on phylogenies: a
  general method for the comparative analysis of discrete
  characters. \emph{Proceedings of the Royal Society of London. Series
    B. Biological Sciences}, \bold{255}, 37--445.

  Schluter, D., Price, T., Mooers, A. O. and Ludwig, D. (1997)
  Likelihood of ancestor states in adaptive radiation. \emph{Evolution},
  \bold{51}, 1699--1711.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{corBrownian}}, \code{\link{corGrafen}},
  \code{\link{corMartins}}, \code{\link{compar.ou}}
}
\examples{
### Just some random data...
data(bird.orders)
x <- rnorm(23)
### Compare the three methods for continuous characters:
ace(x, bird.orders)
ace(x, bird.orders, method = "pic")
ace(x, bird.orders, method = "GLS",
    corStruct = corBrownian(1, bird.orders))
### For discrete characters:
ans <- ace(c(rep(0, 5), rep(1, 18)), bird.orders, type = "d")
#### Showing he likelihood of state "1" on each node:
plot(bird.orders, type = "c", FALSE)
nodelabels(thermo = ans$lik.anc[, 2], cex = 0.75)
}
\keyword{models}

\name{ace}
\alias{ace}
\alias{logLik.ace}
\alias{deviance.ace}
\alias{AIC.ace}
\alias{anova.ace}
\title{Ancestral Character Estimation}
\usage{
ace(x, phy, type = "continuous", method = "ML", CI = TRUE,
    model = if (type == "continuous") "BM" else "ER",
    scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1)
\method{logLik}{ace}(object, ...)
\method{deviance}{ace}(object, ...)
\method{AIC}{ace}(object, ..., k = 2)
\method{anova}{ace}(object, ...)
}
\arguments{
  \item{x}{a vector or a factor.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{type}{the variable type; either \code{"continuous"} or
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
  \item{ip}{the initial value(s) used for the ML estimation procedure
    when \code{type == "discrete"} (possibly recycled).}
  \item{object}{an object of class \code{"ace"}.}
  \item{k}{a numeric value giving the penalty per estimated parameter;
    the default is \code{k = 2} which is the classical Akaike
    information criterion.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  This function estimates ancestral character states, and the associated
  uncertainty, for continuous and discrete characters.

  \code{logLik}, \code{deviance}, and \code{AIC} are generic functions
  used to extract the log-likelihood, the deviance (-2*logLik), or the
  Akaike information criterion of a tree. If no such values are
  available, \code{NULL} is returned.

  \code{anova} is another generic function that is used to compare
  nested models: the significance of the additional parameter(s) is
  tested with likelihood ratio tests. You must ensure that the models
  are effectively nested (if they are not, the results will be
  meaningless). It is better to list the models from the smallest to the
  largest.
}
\details{
  If \code{type = "continuous"}, the default model is Brownian motion
  where characters evolve randomly following a random walk. This model
  can be fitted by maximum likelihood (the default, Schluter et
  al. 1997), least squares (\code{method = "pic"}, Felsenstein 1985), or
  generalized least squares (\code{method = "GLS"}, Martins and Hansen
  1997, Cunningham et al. 1998). In the latter case, the specification
  of \code{phy} and \code{model} are actually ignored: it is instead
  given through a correlation structure with the option \code{corStruct}.

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
  is determined from the data.
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
  \item{se}{if \code{type = "discrete"}, the standard-errors of
    estimated rates.}
  \item{index.matrix}{if \code{type = "discrete"}, gives the indices of
    the \code{rates} in the rate matrix.}
  \item{loglik}{if \code{method = "ML"}, the maximum log-likelihood.}
  \item{lik.anc}{if \code{type = "discrete"}, the scaled likelihoods of
    each ancestral state.}
  \item{call}{the function call.}
}
\references{
  Cunningham, C. W., Omland, K. E. and Oakley, T. H. (1998)
  Reconstructing ancestral character states: a critical
  reappraisal. \emph{Trends in Ecology & Evolution}, \bold{13},
  361--366.

  Felsenstein, J. (1985) Phylogenies and the comparative
  method. \emph{American Naturalist}, \bold{125}, 1--15.

  Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the
  comparative method: a general approach to incorporating phylogenetic
  information into the analysis of interspecific data. \emph{American
    Naturalist}, \bold{149}, 646--667.

  Pagel, M. (1994) Detecting correlated evolution on phylogenies: a
  general method for the comparative analysis of discrete
  characters. \emph{Proceedings of the Royal Society of London. Series
    B. Biological Sciences}, \bold{255}, 37--45.

  Schluter, D., Price, T., Mooers, A. O. and Ludwig, D. (1997)
  Likelihood of ancestor states in adaptive radiation. \emph{Evolution},
  \bold{51}, 1699--1711.
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}, Ben Bolker
\email{bolker@zoo.ufl.edu}}
\seealso{
  \code{\link{corBrownian}}, \code{\link{corGrafen}},
  \code{\link{corMartins}}, \code{\link{compar.ou}},
  \code{\link[stats]{anova}}
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
x <- factor(c(rep(0, 5), rep(1, 18)))
ans <- ace(x, bird.orders, type = "d")
#### Showing the likelihoods on each node:
plot(bird.orders, type = "c", FALSE, label.offset = 1)
co <- c("blue", "yellow")
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 2, adj = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.75)
### An example of the use of the argument `ip':
tr <- character(4)
tr[1] <- "((((t10:5.03,t2:5.03):2.74,(t9:4.17,"
tr[2] <- "t5:4.17):3.60):2.80,(t3:4.05,t7:"
tr[3] <- "4.05):6.53):2.32,((t6:4.38,t1:4.38):"
tr[4] <- "2.18,(t8:2.17,t4:2.17):4.39):6.33);"
tr <- read.tree(text = paste(tr, collapse = ""))
y <- c(rep(1, 6), rep(2, 4))
### The default `ip = 0.1' makes ace fails:
ace(y, tr, type = "d")
ace(y, tr, type = "d", ip = 0.01)
### Surprisingly, using an initial value farther to the
### MLE than the default one works:
ace(y, tr, type = "d", ip = 0.3)
}
\keyword{models}

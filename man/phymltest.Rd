\name{phymltest}
\alias{phymltest}
\alias{print.phymltest}
\alias{summary.phymltest}
\alias{plot.phymltest}
\title{Fits a Bunch of Models with PHYML}
\usage{
phymltest(seqfile, format = "interleaved", itree = NULL,
          exclude = NULL, execname = NULL, path2exec = NULL)
\method{print}{phymltest}(x, ...)
\method{summary}{phymltest}(object, ...)
\method{plot}{phymltest}(x, ...)
}
\arguments{
  \item{seqfile}{a character string giving the name of the file that
    contains the DNA sequences to be analysed by PHYML.}
  \item{format}{a character string specifying the format of the DNA
    sequences: either \code{"interleaved"} (the default), or
    \code{"sequential"}.}
  \item{itree}{a character string giving the name of a file with a tree
    in Newick format to be used as an initial tree by PHYML. If
    \code{NULL} (the default), PHYML uses a ``BIONJ'' tree.}
  \item{exclude}{a vector of mode character giving the models to be
    excluded from the analysis. These must be among those below, and
    follow the same syntax.}
  \item{execname}{a character string specifying the name of the PHYML
    binary file. If \code{NULL} and the operating system is ``windows'',
    the default name \code{"phyml_w32"} is used, otherwise the name must
    be given.}
  \item{path2exec}{a character string giving the path to the PHYML
    binary file. If \code{NULL} the file must be accessible to R (either
    it is in the computer path, or it is in R's working directory).}
  \item{x}{an object of class \code{"phymltest"}.}
  \item{object}{an object of class \code{"phymltest"}.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  This function calls the software PHYML and fits successively 28 models
  of DNA evolution. The results are saved on disk, as PHYML usually
  does, and returned in R as a vector with the log-likelihood value of
  each model.
}
\details{
  The present function has been tested with version 2.4 of PHYML; it
  should also work with version 2.3, but it won't work with version 2.1.

  Under unix-like systems, it seems necessary to run R from a csh or a
  similar shell (sh might not work).

  The user must take care to set correctly the three diffferent paths
  involved here: the path to PHYML's binary, the path to the sequence
  file, and the path to R's working directory. The function should work
  if all three paths are different. Obviously, there should be no problem
  if they are all the same.

  If the usual output files of PHYML already exist, they are not
  deleted and PHYML's results are appended.

  The following syntax is used for the models:

  "X[Y][Z]00[+I][+G]"

  where "X" is the first letter of the author of the model, "Y" and "Z"
  are possibly other co-authors of the model, "00" is the year of the
  publication of the model, and "+I" and "+G" indicates whether the
  presence of invariant sites and/or a gamma distribution of
  substitution rates have been specified. Thus, Kimura's model is
  denoted "K80" and not "K2P". The exception to this rule is the general
  time-reversible model which is simple denoted "GTR" model.

  The seven substitution models used are: "JC69", "K80", "F81", "F84",
  "HKY85", "TN93", and "GTR". These models are then altered by adding
  the "+I" and/or "+G", resulting thus in four variants for each of them
  (e.g., "JC69", "JC69+I", "JC69+G", "JC69+I+G"). Some of these models
  are described in the help page of \code{\link{dist.dna}}.

  When a gamma distribution of substitution rates is specified, four
  categories are used (which is PHYML's default behaviour), and the
  ``alpha'' parameter is estimated from the data.

  For the models with a different substition rate for transitions and
  transversions, these rates are left free and estimated from the data
  (and not constrained with a ratio of 4 as in PHYML's default).
}
\note{
  It is important to note that the models fitted by this function is
  only a small fraction of the models possible with PHYML. For instance,
  it is possible to vary the number of categories in the (discretized)
  gamma distribution of substitution rates, and many parameters can be
  fixed by the user. The results from the present function should rather
  be taken as indicative of a best model.
}
\value{
  \code{phymltest} returns an object of class \code{"phymltest"}: a
  numeric vector with the models as names.

  The \code{print} method prints an object of class \code{"phymltest"}
  as matrix with the name of the models, the number of free parameters,
  the log-likelihood value, and the value of the Akaike information
  criterion (AIC = -2 * loglik + 2 * number of free parameters)

  The \code{summary} method prints all the possible likelihood ratio
  tests for an object of class \code{"phymltest"}.

  The \code{plot} method plots the values of AIC of an object of class
  \code{"phymltest"} on a vertical scale.
}
\references{
  Posada, D. and Crandall, K. A. (2001) Selecting the best-fit model of
  nucleotide substitution. \emph{Systematic Biology}, \bold{50},
  580--601.

  Guindon, S. and Gascuel, O. (2003) A simple, fast, and accurate
  algorithm to estimate large phylogenies by maximum likelihood.
  \emph{Systematic Biology}, \bold{52}, 696--704.
  \url{http://atgc.lirmm.fr/phyml/}
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.tree}}, \code{\link{write.tree}},
  \code{\link{dist.dna}}
}
\examples{
### A "fake" example with random likelihood values: it does not
### make sense, but does not need PHYML and gives you a flavour
### of how the output looks like:
x <- runif(28, -100, -50)
names(x) <- .phymltest.model
class(x) <- "phymltest"
x
summary(x)
plot(x)
### This example needs PHYML, copy/paste or type the
### following commands if you want to try them, eventually
### changing setwd() and the options of phymltest()
\dontrun{
setwd("D:/phyml_v2.4/exe") # under Windows
data(woodmouse)
write.dna(woodmouse, "woodmouse.txt")
X <- phymltest("woodmouse.txt")
X
summary(X)
plot(X)
}
}
\keyword{models}

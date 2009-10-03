\name{ape-defunct}
\alias{klastorin}
\alias{mlphylo}
\alias{DNAmodel}
\alias{sh.test}
\alias{heterozygosity}
\alias{H}
\alias{nuc.div}
\alias{theta.h}
\alias{theta.k}
\alias{theta.s}
\title{Defunct Ape Functions}
\description{
  These functions have been removed from \pkg{ape} or moved to another
  package.
}
\usage{
klastorin(phy)
mlphylo(x, phy, model = DNAmodel(), search.tree = FALSE,
        quiet = FALSE, value = NULL, fixed = FALSE)
DNAmodel(model = "K80", partition = 1,
         ncat.isv = 1, invar = FALSE,
         equal.isv = TRUE, equal.invar = 1)
sh.test(..., x, model = DNAmodel(), B = 100)
heterozygosity(x, variance = FALSE)
H(x, variance = FALSE)
nuc.div(x, variance = FALSE, pairwise.deletion = FALSE)
theta.h(x, standard.error = FALSE)
theta.k(x, n = NULL, k = NULL)
theta.s(s, n, variance = FALSE)
}
\details{
  \code{klastorin} has been removed because it does not seem to be used
  and this helped to clear some internal C code (this function may be
  put back with a different coding).

  \code{mlphylo}, \code{DNAmodel} and \code{sh.test} have been removed:
  see the package \pkg{phangorn} for a much better implementation of
  these methods (and others).

  \code{heterozygosity}, \code{nuc.div}, \code{theta.h}, \code{theta.k}
  and \code{theta.s} have been moved to \pkg{pegas}.
}
\keyword{internal}

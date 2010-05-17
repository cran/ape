\name{ape-defunct}
\alias{ape-defunct}
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
\alias{evolve.phylo}
\alias{plot.ancestral}
\alias{chronogram}
\alias{ratogram}
\alias{NPRS.criterion}
\title{Defunct Ape Functions}
\description{
  These functions have been removed from \pkg{ape} or moved to another
  package.
}
\usage{
klastorin(phy)
mlphylo(...)
DNAmodel(...)
sh.test(...)
heterozygosity(x, variance = FALSE)
H(x, variance = FALSE)
nuc.div(x, variance = FALSE, pairwise.deletion = FALSE)
theta.h(x, standard.error = FALSE)
theta.k(x, n = NULL, k = NULL)
theta.s(s, n, variance = FALSE)
evolve.phylo(phy, value, var)
plot.ancestral(...)
chronogram(...)
ratogram(...)
NPRS.criterion(...)
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

  \code{evolve.phylo} and \code{plot.ancestral} have been deprecated by
  the new function \code{\link{rTraitCont}}.

  \code{chronogram}, \code{ratogram}, and \code{NPRS.criterion} have
  ceased to be maintained: consider using \code{\link{chronopl}}.
}
\keyword{internal}

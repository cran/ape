\name{dist.dna}
\alias{dist.dna}
\alias{dist.dna.JukesCantor}
\alias{dist.dna.TajimaNei}
\alias{dist.dna.Kimura}
\alias{dist.dna.Tamura}
\title{Pairwise Distances from DNA Sequences}
\usage{
dist.dna(x, y = NULL, variance = FALSE, gamma = NULL,
         method = "Kimura", basefreq = NULL, GCcontent = NULL)
dist.dna.JukesCantor(x, y, variance = FALSE, gamma = NULL)
dist.dna.TajimaNei(x, y, variance = FALSE, basefreq = NULL)
dist.dna.Kimura(x, y, variance = FALSE, gamma = NULL)
dist.dna.Tamura(x, y, variance = FALSE, GCcontent = NULL)
}
\arguments{
  \item{x}{either, a vector with a single DNA sequence, or a matrix of
    DNA sequences, or a list of DNA sequences (the latrer can be taken
    from, e.g., \code{read.GenBank}).}
  \item{y}{a vector with a single DNA sequence.}
  \item{gamma}{a value for the gamma parameter which is possibly used to
    apply a gamma correction to the distances (by default \code{gamma =
      NULL} so no correction is applied).}
  \item{variance}{a logical indicating whether to compute the variances
    of the distances; defaults to \code{FALSE} so the variances are not
    computed.}
  \item{method}{a character string specifying the method used to compute
    the distance. Currently four choices are possible: "JukesCantor",
    "TajimaNei", "Kimura" (the default), and "Tamura".}
  \item{basefreq}{the base frequencies to be used in the computations
    (if applicable, i.e. if \code{method = "TajimaNei"}). By default, the
    base frequencies are computed from the whole sample of sequences.}
  \item{GCcontent}{the content in G+C to be used in the computations
    (if applicable, i.e. if \code{method = "Tamura"}). By default, this
    percentage is computed from the whole sample of sequences.}
}
\description{
  These functions compute a matrix of pairwise distances from DNA
  sequences using a model of DNA evolution. Four models are currently
  available.
}
\details{
  For the function \code{dist.dna}, if the argument \code{y} is specified,
  then it is binded to \code{x}, and the distances between all columns
  of the resulting matrix are computed; otherwise, \code{x} must be a
  matrix or a list. The four other functions take two single sequences
  as arguments.

  The function \code{dist.dna} actually calls one of the other function
  depending on the argument \code{method} (by default "Kimura")
  eventually passing the relevant arguments. For instance, specifying a
  value for the option \code{basefreq} has no effect if the option
  \code{method} is set to "Kimura" or "JukesCantor" (the base
  frequencies are assumed to be equal to 0.25 un both models).

  The molecular evolutionary models available through the option
  \code{method} have been extensively described in the literature. A
  brief description is given; more details can be found in the
  References.

  \item{"JukesCantor"}{This model was developed by Jukes and Cantor
  (1969). It assumes that all substitutions (i.e. a change of a base by
  another one) have the same probability. This probability is the same
  for all sites along the DNA sequence. This last assumption can be
  relaxed by assuming that the substition rate varies among site
  following a gamma distribution which parameter must be given by the
  user. By default, no gamma correction is applied. Another assumption
  is that the base frequencies are balanced and thus equal to 0.25.}

  \item{"TajimaNei"}{Tajima and Nei (1984) developed an extension of the
  Jukes--Cantor model which relaxes the assumption of balanced base
  frequencies. The latter are estimated from the data. In the present
  function, the base frequencies are either given by the user, or
  estimated from the data. This allows the user to compute the base
  frequencies from a different (possibly much larger) data set than the
  one (s)he is interested in computing the distances. If the Tajima--Nei
  distances are computed with the function \code{dist.dna} and no base
  frequencies are given (\code{basefreq =  NULL}), then they are
  computed from the whole vectors, matrix, or list given as argument. If
  the distances are computed with the function \code{dist.dna.TajimaNei}
  and no base frequencies are given, then they are computed from both
  vectors given as argument.}

  \item{"Kimura"}{The distance derived by Kimura (1980), sometimes
  referred to as ``Kimura's 2-parameters distance'', has the same underlying
  assumptions than the Jukes--Cantor distance except that two kinds of
  substitutions are considered: transitions (A <-> G, C <-> T), and
  transversions (A <-> C, A <-> T, C <-> G, G <-> T). They are assumed
  to have different probabilities. A transition is the substitution of a
  purine (C, T) by another one, or the substitution of a pyrimidine (A,
  G) by another one. A transversion is the substitution of a purine by a
  pyrimidine, or vice-versa. Both transition and transversion rates are
  the same for all sites along the DNA sequence. Jin and Nei (1990)
  modified the Kimura model to allow for variation among sites following
  a gamma distribution. Like for the Jukes--Cantor model, the gamma parameter
  must be given by the user. By default, no gamma correction is applied.}

  \item{"Tamura"}{Tamura (1992) generalized the Kimura model by relaxing
  the assumption of equal base frequencies. This is done by taking into
  account the bias in G+C content in the sequences. The substitution
  rates are assumed to be the same for all sites along the DNA sequence.}
}
\value{
  a numeric matrix with possibly the names of the individuals (as given
  by the rownames of the argument `x') as colnames and rownames (if
  \code{variance = FALSE}, the default).
}
\note{
  This is a first attempt to implement the computation of DNA
  evolutionary distances in `ape'. Four models are currently available,
  and this is supposed to encompass the most widely used models of DNA
  evolution (particularly Kimura's one with the gamma correction).

  The models of DNA evolution available in `ape' follow somewhat those
  available in the software MEGA (Kumar et al. 2001) with the exception
  that the more complex model by Tamura and Nei (1993) is not available
  in `ape'. This model assumes distinct rates for both kinds of
  transition (A <-> G versus C <-> T), and transversions. The base
  frequencies are not assumed to be equal and are estimated from the
  data. A gamma correction of the inter-site variation in substitution
  rates is also possible (Kumar et al. 2001).

  The package Phylip (Felsenstein 1993) has also a model not considered
  here. This model (called the maximum likelihood, or ML, model by
  Felsenstein) assumes different rates of transition and transition, and
  the base frequencies are allowed to be unequal. This seems to be a
  particular case of the Tamura--Nei model (Tamura Nei 1993) where both
  kinds of transition are equiprobable.

  It is planned to include more complex models of DNA evolution similar
  to those just described in future versions of `ape'.
}
\references{
  Felsenstein, J. (1993) Phylip (Phylogeny Inference Package) version
  3.5c. Department of Genetics, University of Washington.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}

  Jukes, T. H. and C. R. Cantor. (1969) Evolution of protein
  molecules. in \emph{Mammalian Protein Metabolism}, ed. Munro, H. N.,
  pp. 21--132, New York: Academic Press.

  Kimura, M. (1980) A simple method for estimating evolutionary rates of
  base substitutions through comparative studies of nucleotide
  sequences. \emph{Journal of Molecular Evolution}, \bold{16}, 111--120.

  Kumar, S., Tamura, K., Jakobsen, I. B. and Nei, M. (2001) MEGA2:
  Molecular Evolutionary Genetics Analysis software.
  \emph{Bioinformatics}, \bold{17}, 1244--1245.
  \url{http://www.megasoftware.net/}

  Jin, L. and M. Nei (1990) Limitations of the evolutionary parsimony
  method of phylogenetic analysis. \emph{Molecular Biology and
    Evolution}, \bold{7}, 82--102.

  Tajima, F. and Nei., M. (1984) Estimation of evolutionary distance
  between nucleotide sequences. \emph{Molecular Biology and Evolution},
  \bold{1}, 269--285.

  Tamura, K. 1992. Estimation of the number of nucleotide substitutions
  when there are strong transition-transversion and G + C-content
  biases. \emph{Molecular Biology and Evolution}, \bold{9}, 678--687.

  Tamura, K. and M. Nei. 1993. Estimation of the number of nucleotide
  substitutions in the control region of mitochondrial DNA in humans and
  chimpanzees. \emph{Molecular Biology and Evolution}, \bold{10}, 512--526. 
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{dist.gene}}, \code{\link{dist.phylo}}
}
\keyword{manip}
\keyword{multivariate}

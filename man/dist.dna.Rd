\name{dist.dna}
\alias{dist.dna}
\alias{dist.dna.JukesCantor}
\alias{dist.dna.TajimaNei}
\alias{dist.dna.Kimura}
\alias{dist.dna.Tamura}
\alias{dist.dna.TamuraNei}
\title{Pairwise Distances from DNA Sequences}
\usage{
dist.dna(x, y = NULL, variance = FALSE, gamma = NULL,
         method = "Kimura", basefreq = NULL, GCcontent = NULL)
dist.dna.JukesCantor(x, y, variance = FALSE, gamma = NULL)
dist.dna.TajimaNei(x, y, variance = FALSE, basefreq = NULL)
dist.dna.Kimura(x, y, variance = FALSE, gamma = NULL)
dist.dna.Tamura(x, y, variance = FALSE, GCcontent = NULL)
dist.dna.TamuraNei(x, y, variance = FALSE, basefreq = NULL,
                   gamma = NULL)
}
\arguments{
  \item{x}{either, a vector with a single DNA sequence, or a matrix of
    DNA sequences, or a list of DNA sequences (the latter can be taken
    from, e.g., \code{read.GenBank}).}
  \item{y}{a vector with a single DNA sequence.}
  \item{gamma}{a value for the gamma parameter which is possibly used to
    apply a gamma correction to the distances (by default \code{gamma =
      NULL} so no correction is applied).}
  \item{variance}{a logical indicating whether to compute the variances
    of the distances; defaults to \code{FALSE} so the variances are not
    computed.}
  \item{method}{a character string specifying the method used to compute
    the distance. Currently four choices are possible: \code{"JukesCantor"},
    \code{"TajimaNei"}, \code{"Kimura"} (the default), \code{"Tamura"},
    and \code{"TamuraNei"}.}
  \item{basefreq}{the base frequencies to be used in the computations
    (if applicable, i.e. if \code{method = "TajimaNei"}). By default, the
    base frequencies are computed from the whole sample of sequences.}
  \item{GCcontent}{the content in G+C to be used in the computations
    (if applicable, i.e. if \code{method = "Tamura"}). By default, this
    percentage is computed from the whole sample of sequences.}
}
\description{
  These functions compute a matrix of pairwise distances from DNA
  sequences using a model of DNA evolution. Five models are currently
  available.
}
\details{
  For the function \code{dist.dna}, if the argument \code{y} is specified,
  then it is binded to \code{x}, and the distances between all columns
  of the resulting matrix are computed; otherwise, \code{x} must be a
  matrix or a list. The four other functions take two single sequences
  as arguments.

  The function \code{dist.dna} actually calls one of the other function
  depending on the argument \code{method} (by default \code{"Kimura"})
  eventually passing the relevant arguments. For instance, specifying a
  value for the option \code{basefreq} has no effect if the option
  \code{method} is set to "Kimura" or "JukesCantor" (the base
  frequencies are assumed to be equal to 0.25 in both models).

  The molecular evolutionary models available through the option
  \code{method} have been extensively described in the literature. A
  brief description is given below; more details can be found in the
  References.

  \item{``JukesCantor''}{This model was developed by Jukes and Cantor
  (1969). It assumes that all substitutions (i.e. a change of a base by
  another one) have the same probability. This probability is the same
  for all sites along the DNA sequence. This last assumption can be
  relaxed by assuming that the substition rate varies among site
  following a gamma distribution which parameter must be given by the
  user. By default, no gamma correction is applied. Another assumption
  is that the base frequencies are balanced and thus equal to 0.25.}

  \item{``TajimaNei''}{Tajima and Nei (1984) developed an extension of the
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

  \item{``Kimura''}{The distance derived by Kimura (1980), sometimes
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

  \item{``Tamura''}{Tamura (1992) generalized the Kimura model by relaxing
  the assumption of equal base frequencies. This is done by taking into
  account the bias in G+C content in the sequences. The substitution
  rates are assumed to be the same for all sites along the DNA
  sequence.}

  \item{``TamuraNei''}{Tamura and Nei (1993) developed a model which
  assumes distinct rates for both kinds of transition (A <-> G versus C
  <-> T), and transversions. The base frequencies are not assumed to be
  equal and are estimated from the data. A gamma correction of the
  inter-site variation in substitution rates is possible.}
}
\value{
  a numeric matrix with possibly the names of the individuals (as given
  by the rownames of the argument \code{x}) as colnames and rownames (if
  \code{variance = FALSE}, the default), or a list of two matrices names
  \code{distances} and \code{variance}, respectively (if \code{variance =
    TRUE}).
}
\note{
  The models of DNA evolution available in `ape' follow somewhat those
  available in the software MEGA (Kumar et al. 2001).
}
\references{
  Felsenstein, J. (1993) Phylip (Phylogeny Inference Package) version
  3.5c. Department of Genetics, University of Washington.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}

  Jukes, T. H. and Cantor, C. R. (1969) Evolution of protein
  molecules. in \emph{Mammalian Protein Metabolism}, ed. Munro, H. N.,
  pp. 21--132, New York: Academic Press.

  Kimura, M. (1980) A simple method for estimating evolutionary rates of
  base substitutions through comparative studies of nucleotide
  sequences. \emph{Journal of Molecular Evolution}, \bold{16}, 111--120.

  Kumar, S., Tamura, K., Jakobsen, I. B. and Nei, M. (2001) MEGA2:
  Molecular Evolutionary Genetics Analysis software.
  \emph{Bioinformatics}, \bold{17}, 1244--1245.
  \url{http://www.megasoftware.net/}

  Jin, L. and Nei, M. (1990) Limitations of the evolutionary parsimony
  method of phylogenetic analysis. \emph{Molecular Biology and
    Evolution}, \bold{7}, 82--102.

  Tajima, F. and Nei, M. (1984) Estimation of evolutionary distance
  between nucleotide sequences. \emph{Molecular Biology and Evolution},
  \bold{1}, 269--285.

  Tamura, K. (1992) Estimation of the number of nucleotide substitutions
  when there are strong transition-transversion and G + C-content
  biases. \emph{Molecular Biology and Evolution}, \bold{9}, 678--687.

  Tamura, K. and Nei, M. (1993) Estimation of the number of nucleotide
  substitutions in the control region of mitochondrial DNA in humans and
  chimpanzees. \emph{Molecular Biology and Evolution}, \bold{10}, 512--526. 
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{read.GenBank}}, \code{\link{read.dna}}, \code{\link{write.dna}},
  \code{\link{dist.gene}}, \code{\link{dist.phylo}}
}
\keyword{manip}
\keyword{multivariate}

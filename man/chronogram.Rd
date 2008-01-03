\name{chronogram}
\alias{chronogram}

\title{Chronogram Computed by Nonparametric Rate Smoothing}
\usage{
chronogram(phy, scale = 1, expo = 2, minEdgeLength = 1e-06)
}
\arguments{
  \item{phy}{A phylogenetic tree (i.e. an object of class
    \code{"phylo"}), where the branch lengths are measured in substitutions.}
  \item{scale}{Age of the root in the inferred chronogram (default value: 0). }
  \item{expo}{Exponent in the objective function (default value: 2)}
  \item{minEdgeLength}{Minimum edge length in the phylogram (default
    value: 1e-06). If any branch lengths are smaller then they will be
    set to this value.}
}
\description{
 \code{chronogram} computes a chronogram from a phylogram by applying the NPRS
 (nonparametric rate smoothing) algorithm described in Sanderson (1997).
}
\details{
  Please refer to Sanderson (1997) for mathematical details
}
\value{
\code{chronogram} returns an object of class \code{"phylo"}. The branch lengths of this
tree will be clock-like and scaled so that the root node has age 1 (or the value
set by the option \code{scale}
}
\author{
  Gangolf Jobb (\url{http://www.treefinder.de}) and
Korbinian Strimmer (\url{http://www.stat.uni-muenchen.de/~strimmer/})
}
\seealso{
\code{\link{ratogram}}, \code{\link{NPRS.criterion}}.
}
\references{
  Sanderson, M. J. (1997) A nonparametric approach to estimating
    divergence times in the absence of rate constancy. \emph{Molecular
    Biology and Evolution}, \bold{14}, 1218--1231.
}
\examples{
# get tree
data("landplants.newick") # example tree in NH format
tree.landplants <- read.tree(text = landplants.newick)

# plot tree
tree.landplants
plot(tree.landplants, label.offset = 0.001)

# estimate chronogram
chrono.plants <- chronogram(tree.landplants)

# plot
plot(chrono.plants, label.offset = 0.001)

# value of NPRS function for our estimated chronogram
NPRS.criterion(tree.landplants, chrono.plants)
}
\keyword{manip}

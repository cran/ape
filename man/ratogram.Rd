\name{ratogram}
\alias{ratogram}

\title{Ratogram Computed by Nonparametric Rate Smoothing}
\usage{
ratogram(phy, scale = 1, \dots)
}
\arguments{
  \item{phy}{A phylogenetic tree (i.e. an object of class \code{"phylo"}), where
    the branch lengths are measured in substitutions.}
  
  \item{scale}{Age of the root in the chronogram corresponding to the inferred ratogram(default value: 0). }
    
  \item{minEdgeLength}{Minimum edge length in the phylogram (default value: 1e-06). If any branch lengths are
    smaller then they will be set to this value. }
  
  \item{expo = 2}{Exponent in the objective function (default value: 2)}
    
  \item{...}{Any of the above parameters.}
  
}
\description{

 \code{ratogram} computes a ratogram from a phylogram by applying the NPRS
 (nonparametric rate smoothing) algorithm described in Sanderson (1997).
}

\details{
  Please refer to Sanderson (1997) for mathematical details
}

\value{
\code{chronogram} returns an object of class \code{"phylo"}. The branch lengths of this
tree will be the absolute rates estimated for each branch.}

}

\author{Gangolf Jobb (\url{http://www.treefinder.de}) and
Korbinian Strimmer (\url{http://www.stat.uni-muenchen.de/~strimmer/})
}

\seealso{
\code{\link{chronogram}}, 
\code{\link{NPRS.criterion}}.
}


\references{
  Sanderson, M. J. (1997) A nonparametric approach to estimating
    divergence times in the absence of rate constancy. \emph{Molecular
    Biology and Evolution}, \bold{14}, 1218--1231.
}

\examples{
library(ape)

# get tree
data("landplants.newick") # example tree in NH format
tree.landplants <- read.tree(text = landplants.newick)

# plot tree
tree.landplants
plot(tree.landplants, label.offset = 0.001)

# estimate ratogram
rato.plants <- ratogram(tree.landplants)

# plot and write to file
plot(rato.plants, label.offset = 0.001)
write.tree(rato.plants, file = "ratogram.phy")
}
\keyword{manip}

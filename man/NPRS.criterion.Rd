\name{NPRS.criterion}
\alias{NPRS.criterion}

\title{Objective Function Employed in Nonparametric Rate Smoothing}
\usage{
NPRS.criterion(phy, chrono, expo = 2, minEdgeLength = 1e-06)
}
\arguments{
  \item{phy}{A non-clock-like phylogenetic tree (i.e. an object of class \code{"phylo"}), where
    the branch lengths are measured in substitutions.}

  \item{chrono}{A chronogram, i.e. a clock-like tree (i.e. an object of class \code{"phylo"}), where
    the branch lengths are measured in absolute time.}

  \item{expo}{Exponent in the objective function (default value: 2)}
    
  \item{minEdgeLength}{Minimum edge length in the phylogram (default value: 1e-06). If any branch lengths are
    smaller then they will be set to this value. }  
}
\description{

 \code{NPRS.criterion} computes the objective function to be minimized in the NPRS
 (nonparametric rate smoothing) algorithm described in Sanderson (1997).
}

\details{
  Please refer to Sanderson (1997) for mathematical details. Note that
  is is not computationally efficient to optimize the branch lengths in
  a chronogram by using \code{NPRS.criterion} - please use
  \code{\link{chronogram}} instead.
}

\value{
\code{NPRS.criterion} returns the value of the objective function given a phylogram and a chronogram.}

}

\author{Gangolf Jobb (\url{http://www.treefinder.de}) and
Korbinian Strimmer (\url{http://www.stat.uni-muenchen.de/~strimmer/})
}

\seealso{
\code{\link{ratogram}},
\code{\link{chronogram}}.
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

# estimate chronogram
chrono.plants <- chronogram(tree.landplants)

# plot and write to file
plot(chrono.plants, label.offset = 0.001)
write.tree(chrono.plants, file = "chronogram.phy")

# value of NPRS function for our estimated chronogram
NPRS.criterion(tree.landplants, chrono.plants)
}
\keyword{manip}

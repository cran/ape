\name{as.phylo}
\alias{as.phylo}
\alias{as.phylo.hclust}
\alias{as.hclust.phylo}

\title{Conversion between "phylo" and "hclust" trees}
\usage{
as.phylo(hc)
as.phylo.hclust(hc)
as.hclust.phylo(phy)
}
\arguments{
  \item{hc}{an object of class \code{"hclust"}), i.e. a clustering tree 
            computed using  \code{\link{hclust}}.}
  \item{phy}{a phylogenetic tree, i.e. an object of class \code{"phy"}. To be converted
            into a \code{"hclust"} object the phylogenetic tree needs to be  ultra-metric and binary.}
   
}
\description{
 The function \code{as.phylo.hclust} (short form \code{as.phylo}) converts a clustering 
 tree computed using \code{\link{hclust}}
 into a phylogenetic tree, i.e. into an object of class \code{"phylo"}.
 
 The function \code{as.hclust.phylo} (short form \code{as.hclust}) converts a phylogenetic tree into a 
 \code{"hclust"} object.
}
\value{
An object of class \code{"hclust"} or \code{"phylo"}, respectively.

}
\seealso{
\code{\link{hclust}}, \code{\link{as.hclust}}, \code{\link{is.binary.tree}},
\code{\link{is.ultrametric}}.
}

\author{Gangolf Jobb (\url{http://www.treefinder.de})
}

\examples{
library(ape)
library(mva) # for hclust and example data set

# create hclust object
data(USArrests)
hclust.tree<-hclust(dist(USArrests[1:6,]))
plot(hclust.tree)

# convert into phylo object
phylo.tree <- as.phylo(hclust.tree)
phylo.tree
plot(phylo.tree)

# phylogenetic tree is ultrametric and binary 
is.ultrametric(phylo.tree)
is.binary.tree(phylo.tree)

# convert back into hclust object
converted.tree <- as.hclust(phylo.tree)
plot(converted.tree)
}
\keyword{manip}

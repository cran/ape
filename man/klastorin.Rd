\name{klastorin}
\alias{klastorin}
\title{Klastorin's (1982) method for classifying genes as suggested by Misawa
and Tajima (2000)}
\usage{
klastorin(phy)
}
\arguments{
  \item{phy}{a phylogenetic tree, i.e. an object of class \code{"phy"}. The root of the tree should make
   sense biologically.
  }
}
\description{
  The function \code{klastorin} uses the method by Klastorin's (1982) as
  suggested by Misawa and Tajima (2000) for identifying groups within
  gene trees.
}
\value{
A vector indication the class affiliation for each sequence/taxon in the tree.
}
\seealso{
\code{\link{opsin}}.
}
\references{
   Klastorin T.D. (1982) An alternative method for hospital partition
   determination using hierarchical cluster analysis. \emph{Operations
     Research} \bold{30},1134--1147.

   Misawa, K. (2000) A simple method for classifying genes and a bootstrap
   test for classifications. \emph{Molecular Biology and Evolution},
   \bold{17}, 1879--1884.
}
\author{Gangolf Jobb (\url{http://www.treefinder.de})}
\examples{
library(ape)

# find groups in landplant tree
data("landplants.newick")
tree1 <- read.tree(text = landplants.newick)
plot(tree1, label.offset = 0.001)
klastorin(tree1)
tree1$tip.label

# find groups in opsin tree
data("opsin.newick")
tree2 <- read.tree(text = opsin.newick)
plot(tree2,label.offset = 0.01)
groups <- klastorin(tree2)
groups
tree2$tip.label[groups==1]
tree2$tip.label[groups==2]
tree2$tip.label[groups==3]
tree2$tip.label[groups==4]
tree2$tip.label[groups==5]
}
\keyword{manip}

\name{bind.tree}
\alias{bind.tree}
\title{Binds Trees}
\usage{
bind.tree(x, y, node = -1, branch = NULL, position = NULL)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{y}{an object of class \code{"phylo"}.}
  \item{node}{a negative integer giving the node number of the tree
    \code{x} where the tree \code{y} is binded (by default -1, i.e. the
    root).}
  \item{branch}{an integer; if specified this gives the number of a tip
    or of a label of the tree \code{x} below which the tree \code{y} is
    binded.}
  \item{position}{a numeric value giving the position from the tip or
    node given by \code{branch} where the tree \code{y} is binded.}
}
\description{
  This function binds together two phylogenetic trees to result in a
  single object of class \code{"phylo"}.
}
\details{
  The argument \code{node} is ignored if \code{branch} is given. The
  argument \code{x} can be seen as the receptor tree, whereas \code{y}
  is the donor tree. The root of \code{y} is then sticked on a location
  of \code{x} specified by either the argument \code{node}, or the
  arguments \code{branch} and \code{position}. If \code{y} has a root
  edge, this is added as in internal branch in the resulting tree.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{drop.tip}}
}
\examples{
### binds the two clades of bird orders
cat("((Struthioniformes:21.8,Tinamiformes:21.8):4.1,",
    "((Craciformes:21.6,Galliformes:21.6):1.3,Anseriformes:22.9):3.0):2.1;",
    file = "ex1.tre", sep = "\n")
cat("(Turniciformes:27.0,(Piciformes:26.3,((Galbuliformes:24.4,",
    "((Bucerotiformes:20.8,Upupiformes:20.8):2.6,",
    "(Trogoniformes:22.1,Coraciiformes:22.1):1.3):1.0):0.6,",
    "(Coliiformes:24.5,(Cuculiformes:23.7,(Psittaciformes:23.1,",
    "(((Apodiformes:21.3,Trochiliformes:21.3):0.6,",
    "(Musophagiformes:20.4,Strigiformes:20.4):1.5):0.6,",
    "((Columbiformes:20.8,(Gruiformes:20.1,Ciconiiformes:20.1):0.7):0.8,",
    "Passeriformes:21.6):0.9):0.6):0.6):0.8):0.5):1.3):0.7):1.0;",
    file = "ex2.tre", sep = "\n")
tree.bird1 <- read.tree("ex1.tre")
tree.bird2 <- read.tree("ex2.tre")
unlink(c("ex1.tre", "ex2.tre")) # clean-up
birds <- bind.tree(tree.bird1, tree.bird2, branch = -1,
                   position = tree.bird1$root.edge)
birds
layout(matrix(c(1, 2, 3, 3), 2, 2))
plot(tree.bird1)
plot(tree.bird2)
plot(birds)
}
\keyword{manip}

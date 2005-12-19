\name{evolve.phylo}
\alias{evolve.phylo}
\title{Ancestral Character Simulation}
\description{
  Simulate the (independent) evolution of one or several continuous
  characters along a given phylogenetic tree under a homogeneous
  Brownian model.
}
\usage{
evolve.phylo(phy, value, var)
}
\arguments{
  \item{phy}{an object of class 'phylo' with branch lengths.}
  \item{value}{ancestral states, one by character. The (optional) names
    of this vector will be used as character names.}
  \item{var}{the variance of each character.}
}
\details{
  Let x be a continuous character.
  If it evolves according to a Brownian model, its value at time t follows a normal law with mean x0 and variance t*sigma\_x,
  where x0 is the value of the character at time 0, and sigma\_x is the 'inner' variance of the character.
  The evolution of a continuous character is performed by letting the character evolve on each branch, from its ancestral root state.
  The final state of a branch is the ancestral states of the daughter branches, and so on.
}
\value{
  An object of class 'ancestral', inheriting from the 'phylo' class. The
  following components are added:

  \item{node.character}{a data.frame with node ids as rownames and one
    column by character, containing all the inner node values for each
    character.}
  \item{tip.character}{a data.frame with tip ids as rownames and one
    column by character, containing all the tip values for each
    character.}
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
  \code{\link{plot.ancestral}}, \code{\link{ace}}
}
\examples{
data(bird.orders)
x <- rep(0, 5)
names(x) <- c("A", "B", "C", "D", "E")
anc1 <- evolve.phylo(bird.orders, x, 1)
anc2 <- evolve.phylo(bird.orders, x, 1)
cor(anc1$tip.character, anc2$tip.character)
}
\keyword{models}

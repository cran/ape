\name{compute.brlen}
\alias{compute.brlen}
\title{Branch lengths Computation}
\usage{
	compute.brlen(phy, method="Grafen", power=1)
}
\arguments{
	\item{phy}{An object of class \code{phylo} representing the tree.}
	\item{method}{The method to use to compute lengths. Only Grafen's 1989 method is implemented for now.}
	\item{power}{The power at which heights must be raised (see below).}
}
\description{
	Compute branch lengths.
}
\details{
	Grafen's (1989) computation of branch lengths:
	Each node is gviven a 'height', namely the number of leaves of the subtree minus one, 0 for leaves.
	Each height is scaled so that root height is 1, and then raised at power 'rho' (>0).
	Branch lengths are then computed as the difference between height of lower node and height of upper node.
}
\value{
	An object of class \code{phylo} with branch lengths.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\references{
	Grafen, A. (1989) The Phylogenetic Regression. \emph{Philosophical Transactions of the Royal society of London. Series B, Biological Sciences}, \bold{326}(1233), 119--157.
}
\seealso{
  \code{\link{read.tree}} for a description of \code{phylo} objects.
}
\examples{
	data(bird.orders)
	split.screen(c(2,2))
	screen(1)
	plot(compute.brlen(bird.orders, power=1), main=expression(rho==1))
	screen(2)
	plot(compute.brlen(bird.orders, power=3), main=expression(rho==3))
	screen(3)
	plot(compute.brlen(bird.orders, power=0.5), main=expression(rho==0.5))
	screen(4)
	plot(compute.brlen(bird.orders, power=0.1), main=expression(rho==0.1))
	close.screen(all=TRUE)
}
\keyword{manip}

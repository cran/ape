\name{rotate}
\alias{rotate}
\title{Swopping sister clades}
\description{
For a given node, rotate exchanges the position of two clades descending from this node. It can handle dichotomies as well as polytomies. In the latter case, two clades from the polytomy are selected for swopping.}
\usage{
rotate(phy, node, polytom = c(1, 2))
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{node}{a vector of mode numeric or character specifying the number of the node}
  \item{polytom}{a vector of mode numeric and length two specifying the two clades that should be exchanged in a polytomy}
}
\details{
phy can be either rooted or unrooted, contain polytomies and lack branch lengths. In the presence of very short branch lengths it is convenient to plot the phylogenetic tree without branch lengths in order to identify the number of the node in question.

node can be any of the interior nodes of a phylogenetic tree including the root node. Number of the nodes can be identified by the nodelabels function. Alternatively, you can specify a vector of length two that contains either the number or the names of two tips that coalesce in the node of interest.

If the node subtends a polytomy, any two clades of the the polytomy can be chosen by polytom. On a plotted phylogeny, the clades are numbered from bottom to top and polytom is used to index the two clades one likes to swop.
}
\value{
an object of class \code{"phylo"}.
}
\author{Christoph Heibl \email{heibl@lmu.de}}
\seealso{
	\code{\link{plot.phylo}},
	\code{\link{nodelabels}},
	\code{\link{root}},
	\code{\link{drop.tip}}}
\examples{
# create a random tree:
tre <- rtree(25)

# visualize labels of internal nodes:
plot.phylo(tre, use.edge.length=FALSE)
nodelabels()

# rotate clades around node 30:
tre.new <- rotate(tre, 30)

# compare the results:
X11() # open new graphical device
par(mfrow=c(1,2)) # devide graphical device
plot(tre) # plot old tre
plot(tre.new) # plot new tree

# visualize labels of terminal nodes:
X11() # open new graphical device
plot.phylo(tre)
tiplabels()

# rotate clades containing nodes 12 and 20:
tre.new <- rotate(tre, c(12, 21))

# compare the results:
X11() # open new graphical device
par(mfrow=c(1,2)) # devide graphical device
plot(tre) # plot old tre
plot(tre.new) # plot new tree

# or you migth just specify tiplabel names:
tre.new <- rotate(tre, c("t3", "t14"))

# compare the results:
X11() # open new graphical device
par(mfrow=c(1,2)) # devide graphical device
plot(tre) # plot old tre
plot(tre.new) # plot new tree
}
\keyword{manip}

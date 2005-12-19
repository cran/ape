### evolve.tree.R  (2005-12-04)
###
###     Character Simulation under a Brownian Model
###
### Copyright 2005 Julien Dutheil <julien.dutheil@univ-montp2.fr>
###
### This file is part of the `ape' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

evolve.phylo <- function(phy, value, var) {
  if (!("phylo" %in% class(phy)))
      stop("object \"phy\" is not of class \"phylo\"")
  if (is.null(phy$edge.length))
      stop("tree \" phy\" must have branch lengths.")
  nchar <- max(length(value), length(var))
  value <- rep(value, length=nchar)
  var   <- rep(var,   length=nchar)
  char.names <- names(value);

  edges <- phy$edge
  nodes <- unique(as.vector(edges))
  n <- length(nodes) # Number of nodes
  root <- match("-1", nodes)
  states<-list();
  states[["-1"]] <- value
  for(node in nodes[-root]) {
    edge.index <- match(node, edges[,2])
    edge.length <- phy$edge.length[edge.index]
    ancestor <- edges[edge.index, 1]
    ancestor.node.index <- match(ancestor, nodes)
    ancestor.states <- states[[ancestor.node.index]]
    index <- match(node, nodes)
    x <- numeric(nchar)
    for(i in 1:nchar) {
      x[i] <- rnorm(1, mean=ancestor.states[i], sd=sqrt(var[i]*edge.length))
    }
    states[[index]] <- x;
  }
  nodes.states <- as.data.frame(matrix(ncol=nchar, nrow=0))
  if(!is.null(char.names)) names(nodes.states) <- char.names
  count <- 1
  for(i in unique(edges[,1])) {
    nodes.states[i,] <- states[[match(i, nodes)]]
    count <- count + 1
  }

  nl <- length(phy$tip.label) #Number of leaves
  leaves.states <- as.data.frame(matrix(ncol=nchar, nrow=0))
  if(!is.null(char.names)) names(leaves.states) <- char.names
  count <- 1
  for(i in 1:nl) {
    leaves.states[as.character(count),] <- states[[match(as.character(i), nodes)]]
    count <- count + 1
  }

  phy[["node.character"]] <- nodes.states;
  phy[["tip.character"]]  <- leaves.states;
  if(! "ancestral" %in% class(phy)) class(phy) <- c("ancestral", class(phy));
  return(phy)
}

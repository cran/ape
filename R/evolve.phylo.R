### evolve.tree.R (2005-12-04)
###
###   Character Simulation under a Brownian Model
###
### Copyright 2005 Julien Dutheil
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

evolve.phylo <- function(phy, value, var) {
  if (!("phylo" %in% class(phy)))
      stop("object \"phy\" is not of class \"phylo\"")
  if (is.null(phy$edge.length))
      stop("tree \" phy\" must have branch lengths.")
  nchar <- max(length(value), length(var))
  value <- rep(value, length=nchar)
  var   <- rep(var,   length=nchar)
  char.names <- names(value);

  ## added by EP for the new coding of "phylo" (2006-10-04):
  phy <- new2old.phylo(phy)
  ## End
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
  ## added by EP for the new coding of "phylo" (2006-10-04):
  phy <- old2new.phylo(phy)
  ## End
  return(phy)
}

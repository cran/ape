## branching.times.R (2006-10-04)

##    Branching Times of a Phylogenetic Tree

## Copyright 2002-2006 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

branching.times <- function(phy)
{
### the tree must be in cladewise order
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    xx <- numeric(phy$Nnode)
    interns <- which(phy$edge[, 2] > n)
    ## we loop only on the internal edges, this assumes
    ## that `xx' is already set with 0
    for (i in interns)
      xx[phy$edge[i, 2] - n] <- xx[phy$edge[i, 1] - n] + phy$edge.length[i]
    depth <- xx[phy$edge[N, 1] - n] + phy$edge.length[N]
    xx <- depth - xx
    names(xx) <-
      if (is.null(phy$node.label)) (n + 1):(n + phy$Nnode) else phy$node.label
    xx
}

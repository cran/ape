### rotate.R (2006-10-05)
###
###     Rotate an Internal Branch of a Tree
###
### Copyright 2004-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

rotate <- function(phy, group)
{
    if (class(phy) != "phylo")
      stop("object \"phy\" is not of class \"phylo\"")
    nb.tip <- length(phy$tip.label)
    if (length(group) == 1) {
        if (group == "all") {
            ind <- which(phy$edge[, 2] <= nb.tip)
            phy$edge[ind, 2] <- phy$edge[rev(ind), 2]
            phy$tip.label <- rev(phy$tip.label)
        }
        return(phy)
    }
    group <-
      if (is.numeric(group)) sort(group) else which(phy$tip.label %in% group)
    ## Check that the group is monophyletic
    msg <- "the specified group is not monophyletic!"
    if (!all(diff(group) == 1)) stop(msg)
    ## Find the MRCA of the tips given as group
    ## (see `root.R' for comments on the code)
    seq.nod <- .Call("seq_root2tip", phy$edge[, 1], phy$edge[, 2],
                     nb.tip, phy$Nnode, PACKAGE = "ape")
    sn <- seq.nod[group]
    MRCA <- nb.tip + 1
    i <- 2
    repeat {
        x <- unique(unlist(lapply(sn, "[", i)))
        if (length(x) != 1) break
        MRCA <- x
        i <- i + 1
    }
    ## Check that all descendants of this node
    ## are included in the outgroup
    desc <- which(unlist(lapply(seq.nod,
                                function(x) any(x %in% MRCA))))
    if (length(group) != length(desc)) stop(msg)
    if (!all(sort(group) == sort(desc))) stop(msg)
    ind <- which(phy$edge[, 2] %in% group)
    phy$edge[ind, 2] <- phy$edge[rev(ind), 2]
    phy$tip.label[group] <- phy$tip.label[rev(group)]
    phy
}

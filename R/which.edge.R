## which.edge.R (2013-05-10)

##   Identifies Edges of a Tree

## Copyright 2004-2013 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

getMRCA <- function(phy, tip)
### Find the MRCA of the tips given as `tip'
### (see `root.R' for comments on the code)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (is.character(tip)) tip <- which(phy$tip.label %in% tip)
    if (length(tip) < 2) return(NULL)
    Ntip <- length(phy$tip.label)
    seq.nod <- .Call(seq_root2tip, phy$edge, Ntip, phy$Nnode)
    sn <- seq.nod[tip]
    MRCA <- Ntip + 1
    i <- 2
    repeat {
        x <- unique(unlist(lapply(sn, "[", i)))
        if (length(x) != 1) break
        MRCA <- x
        i <- i + 1
    }
    MRCA
}

which.edge <- function(phy, group)
{
    if (!inherits(phy, "phylo"))
      stop('object "phy" is not of class "phylo"')
    if (is.character(group))
      group <- which(phy$tip.label %in% group)
    if (length(group) == 1)
      return(match(group, phy$edge[, 2]))
    nb.tip <- length(phy$tip.label)
    MRCA <- getMRCA(phy, group)
    if (MRCA == nb.tip + 1) {
        from <- 1
        to <- dim(phy$edge)[1]
    } else {
        from <- which(phy$edge[, 2] == MRCA) + 1
        to <- max(which(phy$edge[, 2] %in% group))
    }
    wh <- from:to
    tmp <- phy$edge[wh, 2]
    ## check that there are no extra tips:
    ## (we do this by selecting the tips in `group' and the nodes
    ##  i.e., the internal edges)
    test <- tmp %in% group | tmp > nb.tip
    if (any(!test)) {
        wh <- wh[test] # drop the extra tips
        ## see if there are no extra internal edges:
        tmp <- phy$edge[wh, ]
        test <- !(tmp[, 2] %in% tmp[, 1]) & tmp[, 2] > nb.tip
        while (any(test)){
            wh <- wh[!test]
            tmp <- phy$edge[wh, ]
            test <- !(tmp[, 2] %in% tmp[, 1]) & tmp[, 2] > nb.tip
        }
    }
    wh
}

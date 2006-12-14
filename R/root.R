### root.R (2006-11-29)
###
###      Root of Phylogenetic Trees
###
### Copyright 2004-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

is.rooted <- function(phy)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    if (!is.null(phy$root.edge)) return(TRUE)
    else
      if (tabulate(phy$edge[, 1])[length(phy$tip.label) + 1] > 2)
        return(FALSE)
      else return(TRUE)
}

unroot <- function(phy)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    if (dim(phy$edge)[1] < 3)
      stop("cannot unroot a tree with two edges.")
    ## delete FIRST the root.edge (in case this is sufficient to
    ## unroot the tree, i.e. there is a multichotomy at the root)
    if (!is.null(phy$root.edge)) phy$root.edge <- NULL
    if (!is.rooted(phy)) return(phy)
    ## We remove one of the edges coming from the root, and
    ## eventually adding the branch length to the other one
    ## also coming from the root.
    ## In all cases, the node deleted is the 2nd one (numbered
    ## nb.tip+2 in `edge'), so we simply need to renumber the
    ## nodes by adding 1, except the root (this remains the
    ## origin of the tree).
    nb.tip <- length(phy$tip.label)
    if (phy$edge[1, 2] <= nb.tip) {
        ## If the 1st edge is terminal, we delete the second which
        ## is necessarily internal and coming from the root too:
        j <- 1 # the target where to stick the edge
        i <- 2 # the edge to delete
    } else { # remove the 1st edge
        j <- which(phy$edge[, 1] == -1)[2]
        i <- 1
    }
    phy$edge <- phy$edge[-i, ]
    nodes <- phy$edge > nb.tip + 1 # renumber all nodes except the root
    phy$edge[nodes] <- phy$edge[nodes] - 1
    if (!is.null(phy$edge.length)) {
        phy$edge.length[j] <- phy$edge.length[j] + phy$edge.length[i]
        phy$edge.length <- phy$edge.length[-i]
    }
    phy$Nnode <- phy$Nnode - 1
    if (!is.null(phy$node.label))
      phy$node.label <- phy$node.label[-2]
    phy
}

root <- function(phy, outgroup)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    if (is.character(outgroup))
      outgroup <- which(phy$tip.label %in% outgroup)
    nb.tip <- length(phy$tip.label)
    if (length(outgroup) == nb.tip) return(phy)

    ## First check that the outgroup is monophyletic--
    ## unless there's only one tip specified of course
    ROOT <- nb.tip + 1
    if (length(outgroup) > 1) {
        msg <- "the specified outgroup is not monophyletic!"
        ## If all tips in `tip' are not contiguous, then
        ## no need to go further:
        if (!all(diff(outgroup) == 1)) stop(msg)
        seq.nod <- .Call("seq_root2tip", phy$edge[, 1], phy$edge[, 2],
                         nb.tip, phy$Nnode, PACKAGE = "ape")
        sn <- seq.nod[outgroup]
        ## We go from the root to the tips: the sequence of nodes
        ## is identical until the MRCA:
        newroot <- ROOT
        i <- 2 # we start at the 2nd position since the root
               # of the tree is a common ancestor to all tips
        repeat {
            x <- unique(unlist(lapply(sn, "[", i)))
            if (length(x) != 1) break
            newroot <- x
            i <- i + 1
        }
        ## Check that all descendants of this node
        ## are included in the outgroup
        ## (1st solution... there may be something smarter)
        desc <- which(unlist(lapply(seq.nod,
                                    function(x) any(x %in% newroot))))
        if (length(outgroup) != length(desc)) stop(msg)
        if (!all(sort(outgroup) == sort(desc))) stop(msg)

    } else newroot <- phy$edge[which(phy$edge[, 2] == outgroup), 1]

    if (newroot == ROOT) return(phy)

### <FIXME>
### The remaining part of the code has not been improved; this
### does not seem obvious. This is delayed...     (2006-09-23)
### </FIXME>

    ## Invert all branches from the new root to the old one
    i <- which(phy$edge[, 2] == newroot)
    nod <- phy$edge[i, 1]
    while (nod != ROOT) {
        j <- which(phy$edge[, 2] == nod)
        phy$edge[i, 1] <- phy$edge[i, 2]
        phy$edge[i, 2] <- nod
        i <- j
        nod <- phy$edge[i, 1]
    }

    i.oroot <- which(phy$edge[, 1] == ROOT)
    ## Unroot the tree if there's a basal dichotomy...
    if (length(i.oroot) == 2) {
        j <- i.oroot[which(i.oroot != i)]
        phy$edge[j, 1] <- phy$edge[i, 2]
        phy$edge <- phy$edge[-i, ]
        if (!is.null(phy$edge.length)) {
            phy$edge.length[j] <- phy$edge.length[j] + phy$edge.length[i]
            phy$edge.length <- phy$edge.length[-i]
        }
        phy$edge[which(phy$edge == newroot)] <- ROOT
    } else {
        ## ... otherwise just invert the root with the newroot
        phy$edge[which(phy$edge == newroot)] <- ROOT
        phy$edge[i.oroot] <- newroot
        ## ... and invert finally! (fixed 2005-11-07)
        phy$edge[i, ] <- rev(phy$edge[i, ])
    }
    if (!is.null(phy$node.label)) {
        tmp <- phy$node.label[1]
        phy$node.label[1] <- phy$node.label[-as.numeric(newroot)]
        phy$node.label[-as.numeric(newroot)] <- tmp
    }
    read.tree(text = write.tree(phy, multi.line = FALSE))
}

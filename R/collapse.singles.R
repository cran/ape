## collapse.singles.R (2006-07-15)

##    Collapse "Single" Nodes

## Copyright 2006 Ben Bolker

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

collapse.singles <- function(tree)
{
    elen <- tree$edge.length
    xmat <- tree$edge
    singles <- NA
    while (length(singles) > 0) {
        ## changed by EP to make it slightly more efficient:
        ## tx <- table(xmat[xmat < 0])
        ## singles <- as.numeric(names(tx)[tx < 3])
        tx <- tabulate(xmat[, 1])
        singles <- which(tx == 1)
        ## END
        if (length(singles) > 0) {
            i <- singles[1]
            prev.node <- which(xmat[, 2] == i)
            next.node <- which(xmat[, 1] == i)
            xmat[prev.node, 2] <- xmat[next.node, 2]
            xmat <- xmat[xmat[, 1] != i, ] ## drop
            ## changed by EP for the new coding of "phylo" (2006-10-05):
            ## xmat[xmat < i] <- xmat[xmat < i] + 1 ## adjust indices
            xmat[xmat > i] <- xmat[xmat > i] - 1 ## adjust indices
            ## END
            elen[prev.node] <- elen[prev.node] + elen[next.node]
            elen <- elen[-next.node]
        }
    }
    tree$edge <- xmat
    tree$edge.length <- elen
    tree
}

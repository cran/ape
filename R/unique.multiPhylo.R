## unique.multiPhylo.R (2007-11-16)

##   Revomes Duplicate Trees from a List

## Copyright 2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

unique.multiPhylo <- function(x, incomparables = FALSE,
                              use.edge.length = FALSE,
                              use.tip.label = TRUE, ...)
{
    n <- length(x)
    keep <- !logical(n)
    for (i in 2:n) {
        j <- 1
        while (j < i) {
            if (all.equal(x[[j]], x[[i]],
                          use.edge.length = use.edge.length,
                          use.tip.label = use.tip.label)) {
                keep[i] <- FALSE
                break
            }
            j <- j + 1
        }
    }
    structure(x[keep], class = "multiPhylo")
}

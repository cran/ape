## unique.multiPhylo.R (2008-06-09)

##   Revomes Duplicate Trees from a List

## Copyright 2007-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

unique.multiPhylo <- function(x, incomparables = FALSE,
                              use.edge.length = FALSE,
                              use.tip.label = TRUE, ...)
{
    n <- length(x)
    keep <- 1L
    for (i in 2:n) {
        already.seen <- FALSE
        for (s in x[keep]) {
            if (all.equal(s, x[[i]],
                          use.edge.length = use.edge.length,
                          use.tip.label = use.tip.label)) {
                already.seen <- TRUE
                break
            }
        }
        if (!already.seen) keep <- c(keep, i)
    }
    x[keep]
}

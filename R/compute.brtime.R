## compute.brtime.R (2011-07-26)

##   Compute and Set Branching Times

## Copyright 2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

compute.brtime <-
    function(phy, method = "coalescent", force.positive = NULL)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    n <- length(phy$tip.label)
    m <- phy$Nnode
    N <- Nedge(phy)

    ## x: branching times (aka, node ages or heights)

    if (identical(method, "coalescent")) { # the default
        x <- 2 * rexp(m)/(as.double((m + 1):2) * as.double(m:1))
        ## x <- 2 * rexp(n - 1)/(as.double(n:2) * as.double((n - 1):1))
        if (is.null(force.positive)) force.positive <- TRUE
    } else if (is.numeric(method)) {
        x <- as.vector(method)
        if (length(x) != m)
            stop("number of branching times given is not equal to the number of nodes")
        if (is.null(force.positive))
            force.positive <- FALSE
    }

    y <- c(rep(0, n), x) # for all nodes (terminal and internal)

    phy <- reorder(phy, "pruningwise")
    e1 <- phy$edge[, 1L] # local copies of the pointer
    e2 <- phy$edge[, 2L] #
    el <- numeric(N)
    if (force.positive) y[unique(e1)] <- sort(x)

    for (i in 1:N) el[i] <- y[e1[i]] - y[e2[i]]

    phy$edge.length <- el
    reorder(phy)
}


## plot.phyloExtra.R (2016-02-18)

##   Extra Functions for Plotting and Annotating

## Copyright 2016 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plotBreakLongEdges <- function(phy, n = 1, ...) {
    o <- order(phy$edge.length, decreasing = TRUE)
    i <- o[seq_len(n)]
    phy$edge.length[i] <- max(phy$edge.length[-i])
    plot.phylo(phy, ...)
    edgelabels(edge = i, pch = 19, col = "white")
    edgelabels("//", i, frame = "n")
}

drawSupportOnEdges <- function(value, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    n <- lastPP$Ntip
    m <- lastPP$Nnode
    if (length(value) == m) value <- value[-1]
    else if (length(value) != m - 1)
        stop("incorrect number of support values")
    nodes <- 2:m + n
    i <- match(nodes, lastPP$edge[, 2])
    edgelabels(value, i, ...)
}

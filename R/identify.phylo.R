## identify.phylo.R (2008-02-28)

##   Graphical Identification of Nodes and Tips

## Copyright 2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

identify.phylo <- function(x, nodes = TRUE, tips = FALSE,
                           labels = FALSE, ...)
{
    cat("Click close to a node of the tree...\n")
    xy <- locator(1)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    d <- sqrt((xy$x - lastPP$xx)^2 + (xy$y - lastPP$yy)^2)
    NODE <- which.min(d)
    res <- list()
    if (NODE <= lastPP$Ntip) {
        res$tips <- if (labels) x$tip.label[NODE] else NODE
        return(res)
    }
    if (tips) {
        TIPS <- prop.part(x)[[NODE - lastPP$Ntip]]
        res$tips <- if (labels) x$tip.label[TIPS] else TIPS
    }
    if (nodes) {
        if (is.null(x$node.label)) labels <- FALSE
        res$nodes <- if (labels) x$node.label[NODE - lastPP$Ntip] else NODE
    }
    res
}

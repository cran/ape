## identify.phylo.R (2007-12-14)

##   Graphical Identification of Nodes and Tips

## Copyright 2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

identify.phylo <- function(x, nodes = TRUE, tips = FALSE,
                           labels = FALSE, ...)
{
    xy <- locator(1)
    Ntip <- .last_plot.phylo$Ntip
    d <- sqrt((xy$x - .last_plot.phylo$xx)^2 +
              (xy$y - .last_plot.phylo$yy)^2)
    NODE <- which.min(d)
    print(NODE)
    res <- list()
    if (NODE <= Ntip) {
        res$tips <- if (labels) x$tip.label[NODE] else NODE
        return(res)
    }
    if (tips) {
        TIPS <- prop.part(x)[[NODE - Ntip]]
        res$tips <- if (labels) x$tip.label[TIPS] else TIPS
    }
    if (nodes) {
        if (is.null(x$node.label)) labels <- FALSE
        res$nodes <- if (labels) x$node.label[NODE - Ntip] else NODE
    }
    res
}

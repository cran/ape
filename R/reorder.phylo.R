## reorder.phylo.R (2012-10-12)

##   Internal Reordering of Trees

## Copyright 2006-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

reorder.phylo <- function(x, order = "cladewise", index.only = FALSE, ...)
{
    ORDER <- c("cladewise", "postorder", "pruningwise")
    io <- pmatch(order, ORDER)
    if (is.na(io)) stop("ambiguous order")
    order <- ORDER[io]
    nb.edge <- dim(x$edge)[1]
    if (!is.null(attr(x, "order")))
        if (attr(x, "order") == order)
            if (index.only) return(1:nb.edge) else return(x)
    nb.node <- x$Nnode
    if (nb.node == 1) return(x)
    nb.tip <- length(x$tip.label)
    if (io == 3) {
        x <- reorder(x)
        neworder <-
            .C("neworder_pruningwise", as.integer(nb.tip),
               as.integer(nb.node), as.integer(x$edge[, 1]),
               as.integer(x$edge[, 2]), as.integer(nb.edge),
               integer(nb.edge), PACKAGE = "ape")[[6]]
    } else {
        neworder <-
            .C("neworder_phylo", as.integer(nb.tip),
               as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
               as.integer(nb.edge), integer(nb.edge), io,
               DUP = FALSE, NAOK = TRUE, PACKAGE = "ape")[[5]]
    }
    if (index.only) return(neworder)
    x$edge <- x$edge[neworder, ]
    if (!is.null(x$edge.length))
        x$edge.length <- x$edge.length[neworder]
    attr(x, "order") <- order
    x
}


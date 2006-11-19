### reorder.phylo.R (2006-11-08)
###
###     Internal Reordering of Trees
###
### Copyright 2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

reorder.phylo <- function(x, order = "cladewise", ...)
{
    order <- match.arg(order, c("cladewise", "pruningwise"))
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    nb.edge <- dim(x$edge)[1]
    neworder <- if (order == "cladewise")
      .C("neworder_cladewise", as.integer(nb.tip),
         as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
         as.integer(nb.edge), integer(nb.edge),
         PACKAGE = "ape")[[5]]
    else
      .C("neworder_pruningwise", as.integer(nb.tip),
         as.integer(nb.node), as.integer(x$edge[, 1]),
         as.integer(x$edge[, 2]), as.integer(nb.edge),
         integer(nb.edge), PACKAGE = "ape")[[6]]
    x$edge <- x$edge[neworder, ]
    if (!is.null(x$edge.length))
      x$edge.length <- x$edge.length[neworder]
    attr(x, "order") <- order
    x
}

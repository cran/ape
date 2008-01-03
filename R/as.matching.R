## as.matching.R (2007-12-23)

##    Conversion Between Phylo and Matching Objects

## Copyright 2005-2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

as.matching <- function(x, ...) UseMethod("as.matching")

as.matching.phylo <- function(x, labels = TRUE, ...)
{
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    if (nb.tip != nb.node + 1)
      stop("the tree must be dichotomous AND rooted.")
    x <- reorder(x, "pruningwise")
    mat <- matrix(x$edge[, 2], ncol = 2, byrow = TRUE)
    nodes <- x$edge[seq(by = 2, length.out = nb.node), 1]
    ## we can use match() becoz each node appears once in `mat'
    O <- match(mat, nodes)
    new.nodes <- 1:nb.node + nb.tip
    sel <- !is.na(O)
    mat[sel] <- new.nodes[O[sel]]
    mat <- cbind(t(apply(mat, 1, sort)), new.nodes, deparse.level = 0)

    obj <- list(matching = mat)
    if (!is.null(x$edge.length))
        warning("branch lengths have been ignored")
    if (labels) {
        obj$tip.label <- x$tip.label
        if (!is.null(x$node.label))
          obj$node.label <- x$node.label[match(new.nodes, nodes)]
    }
    class(obj) <- "matching"
    obj
}

as.phylo.matching <- function(x, ...)
{
    N <- 2*dim(x$matching)[1]
    edge <- matrix(NA, N, 2)
    nb.tip <- (N + 2)/2
    nb.node <- nb.tip - 1
    new.nodes <- numeric(N + 1)
    new.nodes[N + 1] <- nb.tip + 1
    nextnode <- nb.tip + 2
    j <- 1
    for (i in nb.node:1) {
        edge[j:(j + 1), 1] <- new.nodes[x$matching[i, 3]]
        for (k in 1:2) {
            if (x$matching[i, k] > nb.tip) {
                edge[j + k - 1, 2] <- new.nodes[x$matching[i, k]] <- nextnode
                nextnode <- nextnode + 1
            } else edge[j + k - 1, 2] <- x$matching[i, k]
        }
        j <- j + 2
    }
    obj <- list(edge = edge)
    if (!is.null(x$tip.label)) obj$tip.label <- x$tip.label
    else obj$tip.label <- as.character(1:nb.tip)
    class(obj) <- "phylo"
    read.tree(text = write.tree(obj))
}

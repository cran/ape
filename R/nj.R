### nj.R (2006-09-15)
###
###      Neighbor-Joining Tree Estimation
###
### Copyright 2004-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

nj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    edge1 <- edge2 <- integer(2*N - 3)
    edge.length <- numeric(2*N - 3)
    ans <- .C("nj", as.double(X), as.integer(N), as.integer(edge1),
              as.integer(edge2), as.double(edge.length), PACKAGE = "ape")
    obj <- list(edge = cbind(ans[[3]], ans[[4]]),
                edge.length = ans[[5]], tip.label = labels)
    obj$Nnode <- N - 2
    class(obj) <- "phylo"
    reorder(obj)
}

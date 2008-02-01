## me.R (2008-01-12)

##   Tree Estimation Based on Minimum Evolution Algorithm

## Copyright 2007 Vincent Lefort with modifications by
##                 Emmanuel Paradis (2008)

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.fastme <- function(X, nni, lib)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- as.integer(attr(X, "Size"))
    labels <- sprintf("%6s", 1:N)
    edge1 <- edge2 <- integer(2*N - 3)
    ans <- .C(lib, as.double(X), N, labels, as.integer(nni),
              edge1, edge2, double(2*N - 3), character(N),
              PACKAGE = "ape")
    labels <- substr(ans[[8]], 1, 6)
    LABS <- attr(X, "Labels")
    labels <- if (!is.null(LABS)) LABS[as.numeric(labels)]
    else gsub("^ ", "", labels)
    structure(list(edge =  cbind(ans[[5]], ans[[6]]), edge.length = ans[[7]],
                   tip.label = labels, Nnode = N - 2),
              class = "phylo")
}

fastme.bal <- function(X, nni = TRUE) .fastme(X, nni, "me_b")

fastme.ols <- function(X, nni = TRUE) .fastme(X, nni, "me_o")

bionj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- as.integer(attr(X, "Size"))
    labels <- sprintf("%6s", 1:N)
    edge1 <- edge2 <- integer(2*N - 3)
    ans <- .C("bionj", as.double(X), N, labels, edge1, edge2,
              double(2*N - 3), character(N), PACKAGE = "ape")
    labels <- substr(ans[[7]], 1, 6)
    LABS <- attr(X, "Labels")
    labels <- if (!is.null(LABS)) LABS[as.numeric(labels)]
    else gsub("^ ", "", labels)
    structure(list(edge =  cbind(ans[[4]], ans[[5]]), edge.length = ans[[6]],
                   tip.label = labels, Nnode = N - 2),
              class = "phylo")
}

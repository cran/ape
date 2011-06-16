## me.R (2011-05-12)

##   Tree Estimation Based on Minimum Evolution Algorithm

## Copyright 2007 Vincent Lefort with modifications by
##                 Emmanuel Paradis (2008-2011)

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

fastme.bal <- function(X, nni = TRUE, spr = TRUE, tbr = TRUE)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- as.integer(attr(X, "Size"))
    labels <- sprintf("%6s", 1:N)
    edge1 <- edge2 <- integer(2*N - 3)

    ans <- .C("me_b", as.double(X), N, labels, as.integer(nni),
              as.integer(spr), as.integer(tbr), edge1, edge2,
              double(2*N - 3), character(N), PACKAGE = "ape")

    labels <- substr(ans[[10]], 1, 6)
    LABS <- attr(X, "Labels")
    labels <- if (!is.null(LABS)) LABS[as.numeric(labels)]
    else gsub("^ ", "", labels)
    structure(list(edge =  cbind(ans[[7]], ans[[8]]), edge.length = ans[[9]],
                   tip.label = labels, Nnode = N - 2L),
              class = "phylo")
}

fastme.ols <- function(X, nni = TRUE)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- as.integer(attr(X, "Size"))
    labels <- sprintf("%6s", 1:N)
    edge1 <- edge2 <- integer(2*N - 3)
    ans <- .C("me_o", as.double(X), N, labels, as.integer(nni),
              edge1, edge2, double(2*N - 3), character(N),
              PACKAGE = "ape")
    labels <- substr(ans[[8]], 1, 6)
    LABS <- attr(X, "Labels")
    labels <- if (!is.null(LABS)) LABS[as.numeric(labels)]
    else gsub("^ ", "", labels)
    structure(list(edge =  cbind(ans[[5]], ans[[6]]), edge.length = ans[[7]],
                   tip.label = labels, Nnode = N - 2L),
              class = "phylo")
}

bionj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    if (any(X > 100))
        stop("at least one distance was greater than 100")
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
                   tip.label = labels, Nnode = N - 2L),
              class = "phylo")
}

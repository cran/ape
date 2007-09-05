## me.R (2007-07-12)

##      Tree Estimation Based on Minimum Evolution Algorithm

## Copyright 2007 Vincent Lefort

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

fastme.bal <- function(X, nni = TRUE)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C("me_b", as.double(X), as.integer(N), as.character(labels),
              "", as.integer(nni), PACKAGE = "ape")
    read.tree(text = ans[[4]])
}

fastme.ols <- function(X, nni = TRUE)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C("me_o", as.double(X), as.integer(N), as.character(labels),
              "", as.integer(nni), PACKAGE = "ape")
    read.tree(text = ans[[4]])
}

bionj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C("bionj", as.double(X), as.integer(N),
              as.character(labels), "", PACKAGE = "ape")
    read.tree(text = ans[[4]])
}

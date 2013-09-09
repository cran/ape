## additive.R (2011-10-11)

##   Incomplete Distance Matrix Filling

## Copyright 2011 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.addit_ultra <- function(libs, X)
{
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)] <- -1
    X[X < 0] <- -1
    X[is.nan(X)] <- -1
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    m <- sum(X == -1)
    ans <- .C(libs, as.double(X), as.integer(N),
              as.integer(m), double(N*N), PACKAGE = "ape")
    matrix(ans[[4]], N, N)
}

additive <- function(X) .addit_ultra(C_additive, X)

ultrametric <- function(X) .addit_ultra(C_ultrametric, X)

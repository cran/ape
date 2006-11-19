### dist.gene.R (2002-08-28)
###
###    Pairwise Distances from Genetic Data
###
### Copyright 2002 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

dist.gene.pairwise <- function(x, variance = FALSE)
{
    if (is.data.frame(x)) x <- as.matrix(x)
    L <- ncol(x)
    n <- nrow(x)
    D <- matrix(NA, n, n)
    diag(D) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            D[i, j] <- D[j, i] <- L - sum(x[i, ] == x[j, ])
        }
    }
    if (!is.null(rownames(x))) rownames(D) <- colnames(D) <- rownames(x)
    if (variance) {
        var.D <- D * (L - D) / L
        return(list(distances = D, variance = var.D))
    }
    else return(D)
}

dist.gene.percentage <- function(x, variance = FALSE)
{
    L <- ncol(x)
    D <- dist.gene.pairwise(x) / L
    if (variance) {
        var.D <- D * (1 - D) / L
        return(list(pairwise.distances = D, variance = var.D))
    }
    else return(D)
}

dist.gene <- function(x, method = "pairwise", variance = FALSE)
{
    if (method == "pairwise")
      return(dist.gene.pairwise(x, variance = variance))
    if (method == "percentage")
      return(dist.gene.percentage(x, variance = variance))
}

## vcv.phylo.R (2012-02-09)

##   Phylogenetic Variance-Covariance or Correlation Matrix

## Copyright 2002-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

vcv <- function(phy, ...) UseMethod("vcv")

vcv.phylo <- function(phy, model = "Brownian", corr = FALSE, ...)
{
    if (is.null(phy$edge.length))
        stop("the tree has no branch lengths")

    pp <- prop.part(phy)

    phy <- reorder(phy, "pruningwise")

    n <- length(phy$tip.label)
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    EL <- phy$edge.length

    ## xx: vecteur donnant la distance d'un noeud
    ##     ou d'un tip à partir de la racine
    ## (same than in is.ultrametric)
    xx <- numeric(n + phy$Nnode)

    vcv <- matrix(0, n, n)

    ## the loop below starts from the bottom of the edge matrix, so
    ## from the root

    for (i in length(e1):1) {
        var.cur.node <- xx[e1[i]]
        xx[e2[i]] <- var.cur.node + EL[i] # sets the variance
        j <- i - 1L
        while (e1[j] == e1[i] && j > 0) {
            left <- if (e2[j] > n) pp[[e2[j] - n]] else e2[j]
            right <- if (e2[i] > n) pp[[e2[i] - n]] else e2[i]
            vcv[left, right] <- vcv[right, left] <- var.cur.node # sets the covariance
            j <- j - 1L
        }
    }

    diag(vcv) <- xx[1:n]

    if (corr) {
        ## This is inspired from the code of `cov2cor' (2005-09-08):
        M <- vcv
        Is <- sqrt(1/M[1 + 0:(n - 1)*(n + 1)])
        vcv[] <- Is * M * rep(Is, each = n)
        vcv[1 + 0:(n - 1)*(n + 1)] <- 1
    }

    dimnames(vcv)[1:2] <- list(phy$tip.label)

    vcv
}

vcv.corPhyl <- function(phy, corr = FALSE, ...)
{
    labels <- attr(phy, "tree")$tip.label
    dummy.df <- data.frame(1:length(labels), row.names = labels)
    res <- nlme::corMatrix(Initialize.corPhyl(phy, dummy.df), corr = corr)
    dimnames(res) <- list(labels, labels)
    res
}

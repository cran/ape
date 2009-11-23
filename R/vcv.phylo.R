## vcv.phylo.R (2009-11-19)

##   Phylogenetic Variance-Covariance or Correlation Matrix

## Copyright 2002-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

vcv <- function(phy, ...) UseMethod("vcv")

vcv.phylo <- function(phy, model = "Brownian", corr = FALSE, ...)
{
    if (is.null(phy$edge.length))
      stop("the tree has no branch lengths")

    foo <- function(node, var, endofclade) {
        ## First, get the extent of clade descending
        ## from `node' in the matrix `edge':
        from <- which(phy$edge[, 1] == node)
        to <- c(from[-1] - 1, endofclade)
        ## Get the #'s of the descendants of `node':
        desc <- phy$edge[from, 2]
        ## The variance of each of these is easy:
        vcv[desc, desc] <<- var + phy$edge.length[from]
        nd <- length(desc)
        ## The variance of `node' is equal to the covariance of
        ## each possible pair among its descendant clades.
        for (i in 1:(nd - 1))
          for (j in (i + 1):nd)
            for (k in phy$edge[from[i]:to[i], 2])
              for (l in phy$edge[from[j]:to[j], 2])
                vcv[k, l] <<- vcv[l, k] <<- var
        for (i in 1:nd) {
            if (desc[i] <= n) next
            foo(desc[i], vcv[desc[i], desc[i]], to[i])
        }
    }

    phy <- reorder(phy)
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    N <- n.node + n
    vcv <- matrix(0, N, N)
    foo(n + 1, 0, dim(phy$edge)[1])

    vcv <- vcv[1:n, 1:n]
    if (corr) {
        ## This is inspired from the code of `cov2cor' (2005-09-08):
        M <- vcv
        Is <- sqrt(1/M[1 + 0:(n - 1)*(n + 1)])
        vcv[] <- Is * M * rep(Is, each = n)
        vcv[1 + 0:(n - 1)*(n + 1)] <- 1
    }
    rownames(vcv) <- colnames(vcv) <- phy$tip.label
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

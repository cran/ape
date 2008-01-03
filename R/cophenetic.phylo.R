## cophenetic.phylo.R (2007-01-23)

##   Pairwise Distances from a Phylogenetic Tree

## Copyright 2006-2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

dist.nodes <- function(x)
{
    if (is.null(x$edge.length))
      stop("your tree has no branch lengths")

    if (!is.binary.tree(x) || !is.rooted(x))
      x <- multi2di(x, random = FALSE)
    n <- length(x$tip.label)
    n.node <- x$Nnode
    N <- n + n.node
    x <- reorder(x, order = "pruningwise")

    res <- matrix(NA, N, N)
    res[cbind(1:N, 1:N)] <- 0 # implicit mode conversion

    ## I like the simplicity of this one:
    res[x$edge] <- res[x$edge[, 2:1]] <- x$edge.length

    ## compute the distances ...
    for (i in seq(from = 1, by = 2, length.out = n.node)) {
        j <- i + 1
        anc <- x$edge[i, 1]
        des1 <- x$edge[i, 2]
        des2 <- x$edge[j, 2]

        ## If `des1' is a node, we look for the nodes and tips for
        ## which the distance up to `des1' has already been
        ## computed, including `des1' itself. For all these, we can
        ## compute the distance up to `anc' and all node(s) and
        ## tip(s) in `des2'.
        if (des1 > n) des1 <- which(!is.na(res[des1, ]))

        ## id. for `des2'
        if (des2 > n) des2 <- which(!is.na(res[des2, ]))

        ## The following expression is vectorized only on `des2' and
        ## not on `des1' because they may be of different lengths.
        for (y in des1)
          res[y, des2] <- res[des2, y] <- res[anc, y] + res[anc, des2]
        ## compute the distances between the tip(s) and node(s)
        ## in `des2' and the ancestor of `anc'; id. for `des2'
        ## (only if it is not the root)
        if (anc != n + 1) {
            ind <- which(x$edge[, 2] == anc)
            nod <- x$edge[ind, 1] # the ancestor of `anc'
            l <- x$edge.length[ind]
            res[des2, nod] <- res[nod, des2] <- res[anc, des2] + l
            res[des1, nod] <- res[nod, des1] <- res[anc, des1] + l
        }
    }
    dimnames(res)[1:2] <- list(1:N)
    res
}

cophenetic.phylo <- function(x)
{
    n <- length(x$tip.label)
    ans <- dist.nodes(x)[1:n, 1:n]
    dimnames(ans)[1:2] <- list(x$tip.label)
    ans
}

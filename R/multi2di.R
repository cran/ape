## multi2di.R (2008-03-17)

##   Collapse and Resolve Multichotomies

## Copyright 2005-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

multi2di <- function(phy, random = TRUE)
{
    degree <- tabulate(phy$edge[, 1])
    target <- which(degree > 2)
    if (!length(target)) return(phy)
    nb.edge <- dim(phy$edge)[1]
    nextnode <- length(phy$tip.label) + phy$Nnode + 1
    new.edge <- edge2delete <- NULL
    wbl <- FALSE
    if (!is.null(phy$edge.length)) {
        wbl <- TRUE
        new.edge.length <- NULL
    }

    for (node in target) {
        ind <- which(phy$edge[, 1] == node)
        N <- length(ind)
        desc <- phy$edge[ind, 2]
        if (random) {
          ## if we shuffle the descendants, we need to eventually
          ## reorder the corresponding branch lenghts (see below)
          ## so we store the result of sample()
            tmp <- sample(length(desc))
            desc <- desc[tmp]
            res <- rtree(N)$edge
        } else {
            res <- matrix(0, 2*N - 2, 2)
            res[, 1] <- N + rep(1:(N - 1), each = 2)
            res[, 2] <- N + rep(2:N, each = 2)
            res[seq(1, by = 2, length.out = N - 1), 2] <- 1:(N - 1)
            res[length(res)] <- N
        }
        if (wbl) {
            ## keep the branch lengths coming from `node'
            el <- numeric(dim(res)[1]) # initialized with 0's
            el[res[, 2] <= N] <-
              if (random) phy$edge.length[ind][tmp] else phy$edge.length[ind]
        }
        ## now substitute the nodes in `res'
        ## `node' stays at the "root" of these new
        ## edges whereas their "tips" are `desc'
        Nodes <- c(node, seq(from = nextnode, length.out = N - 2))
        res[, 1] <- Nodes[res[, 1] - N]
        tmp <- res[, 2] > N
        res[tmp, 2] <- Nodes[res[tmp, 2] - N]
        res[!tmp, 2] <- desc[res[!tmp, 2]]
        new.edge <- rbind(new.edge, res)
        edge2delete <- c(edge2delete, ind)
        if (wbl) new.edge.length <- c(new.edge.length, el)
        nextnode <- nextnode + N - 2
        phy$Nnode <- phy$Nnode + N - 2
    }
    phy$edge <- rbind(phy$edge[-edge2delete, ], new.edge)
    if (wbl)
      phy$edge.length <- c(phy$edge.length[-edge2delete], new.edge.length)
    if (!is.null(attr(phy, "order"))) attr(phy, "order") <- NULL
    print(phy$node.label)
    if (!is.null(phy$node.label))
        phy$node.label <-
            c(phy$node.label, rep("", phy$Nnode - length(phy$node.label)))
    print(phy$node.label)
    reorder(phy)
    ##read.tree(text = write.tree(phy))
}

di2multi <- function(phy, tol = 1e-8)
{
    if (is.null(phy$edge.length)) stop("the tree has no branch length")
    ## We select only the internal branches which are
    ## significantly small:
    ind <- which(phy$edge.length < tol & phy$edge[, 2] > length(phy$tip.label))
    n <- length(ind)
    if (!n) return(phy)
    ## recursive function to `propagate' node #'s in case
    ## there is a series of consecutive edges to remove
    foo <- function(ancestor, des2del) {
        wh <- which(phy$edge[, 1] == des2del)
        for (k in wh) {
            if (phy$edge[k, 2] %in% node2del) foo(ancestor, phy$edge[k, 2])
            else phy$edge[k, 1] <<- ancestor
        }
    }
    node2del <- phy$edge[ind, 2]
    anc <- phy$edge[ind, 1]
    for (i in 1:n) {
        if (anc[i] %in% node2del) next
        foo(anc[i], node2del[i])
    }
    phy$edge <- phy$edge[-ind, ]
    phy$edge.length <- phy$edge.length[-ind]
    phy$Nnode <- phy$Nnode - n
    ## Now we renumber the nodes that need to be:
    sel <- phy$edge > min(node2del)
    for (i in which(sel))
      phy$edge[i] <- phy$edge[i] - sum(node2del < phy$edge[i])
    if (!is.null(phy$node.label))
        phy$node.label <- phy$node.label[-(node2del - length(phy$tip.label))]
    phy
}

## rtree.R (2007-12-22)

##   Generates Random Trees

## Copyright 2004-2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

rtree <- function(n, rooted = TRUE, tip.label = NULL, br = runif, ...)
{
    foo <- function(n, pos) {
        n1 <- .Internal(sample(n - 1, 1, FALSE, NULL))
        n2 <- n - n1
        po2 <- pos + 2*n1 - 1
        edge[c(pos, po2), 1] <<- nod
        nod <<- nod + 1
        if (n1 > 2) {
            edge[pos, 2] <<- nod
            foo(n1, pos + 1)
        } else if (n1 == 2) {
            edge[c(pos + 1, pos + 2), 1] <<- edge[pos, 2] <<- nod
            nod <<- nod + 1
        }
        if (n2 > 2) {
            edge[po2, 2] <<- nod
            foo(n2, po2 + 1)
        } else if (n2 == 2) {
            edge[c(po2 + 1, po2 + 2), 1] <<- edge[po2, 2] <<- nod
            nod <<- nod + 1
        }
    }

    if (n < 2) stop("a tree must have at least 2 tips.")
    nbr <- 2 * n - 2
    if (!rooted) nbr <- nbr - 1
    edge <- matrix(NA, nbr, 2)

    if (n == 2) {
        if (rooted) edge[] <- c(3, 3, 1, 2)
        else stop("an unrooted tree must have at least 3 tips.")
    } else if (n == 3) {
        edge[] <-
          if (rooted) c(4, 5, 5, 4, 5, 1:3)
          else c(4, 4, 4, 1:3)
    } else if (n == 4 && !rooted) {
        edge[] <- c(5, 6, 6, 5, 5, 6, 1:4)
    } else {
        nod <- n + 1
        if (rooted) { # n > 3
            foo(n, 1)
            ## The following is slightly more efficient than affecting the
            ## tip numbers in foo(): the gain is 0.006 s for n = 1000.
            i <- which(is.na(edge[, 2]))
            edge[i, 2] <- 1:n
        } else { # n > 4
            n1 <- .Internal(sample(n - 2, 1, FALSE, NULL))
            if (n1 == n - 2) {
                n2 <- n3 <- 1
            } else {
                n2 <- .Internal(sample(n - n1 - 1, 1, FALSE, NULL))
                n3 <- n - n1 - n2
            }
            po2 <- 2*n1
            po3 <- 2*(n1 + n2) - 1
            edge[c(1, po2, po3), 1] <- nod
            nod <- nod + 1
            if (n1 > 2) {
                edge[1, 2] <- nod
                foo(n1, 2)
            } else if (n1 == 2) {
                edge[2:3, 1] <- edge[1, 2] <- nod
                nod <- nod + 1
            }
            if (n2 > 2) {
                edge[po2, 2] <- nod
                foo(n2, po2 + 1)
            } else if (n2 == 2) {
                edge[c(po2 + 1, po2 + 2), 1] <- edge[po2, 2] <- nod
                nod <- nod + 1
            }
            if (n3 > 2) {
                edge[po3, 2] <- nod
                foo(n3, po3 + 1)
            } else if (n3 == 2) {
                edge[c(po3 + 1, po3 + 2), 1] <- edge[po3, 2] <- nod
                ## nod <- nod + 1
            }
            i <- which(is.na(edge[, 2]))
            edge[i, 2] <- 1:n
        }
    }
    phy <- list(edge = edge)
    phy$tip.label <-
      if (is.null(tip.label)) paste("t", sample(n), sep = "")
      else sample(tip.label)
    if (is.function(br)) phy$edge.length <- br(nbr, ...)
    phy$Nnode <- if (rooted) n - 1 else n - 2
    class(phy) <- "phylo"
    phy
}

rcoal <- function(n, tip.label = NULL, br = rexp, ...)
{
    nbr <- 2*n - 2
    edge <- matrix(NA, nbr, 2)
    x <- br(n - 1, ...) # coalescence times
    if (n == 2) {
        edge[] <- c(3, 3, 1:2)
        edge.length <- rep(x, 2)
    } else if (n == 3) {
        edge[] <- c(4, 5, 5, 4, 5, 1:3)
        edge.length <- c(x[2], x[1], x[1], sum(x))
    } else {
        edge.length <- numeric(nbr)
        h <- numeric(2*n - 1) # initialized with 0's
        node.height <- cumsum(x)
        pool <- 1:n
        nextnode <- 2*n - 1
        for (i in 1:(n - 1)) {
            y <- sample(pool, size = 2)
            ind <- (i - 1)*2 + 1:2
            edge[ind, 2] <- y
            edge[ind, 1] <- nextnode
            edge.length[ind] <- node.height[i] - h[y]
            h[nextnode] <- node.height[i]
            pool <- c(pool[! pool %in% y], nextnode)
            nextnode <- nextnode - 1
        }
    }
    phy <- list(edge = edge, edge.length = edge.length)
    phy$tip.label <-
      if (is.null(tip.label)) paste("t", 1:n, sep = "")
      else tip.label
    phy$Nnode <- n - 1
    class(phy) <- "phylo"
    ##reorder(phy)
    ## to avoid crossings when converting with as.hclust:
    read.tree(text = write.tree(phy))
}

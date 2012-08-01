## root.R (2011-08-05)

##   Root of Phylogenetic Trees

## Copyright 2004-2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.rooted <- function(phy)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (!is.null(phy$root.edge)) TRUE
    else
        if (tabulate(phy$edge[, 1])[length(phy$tip.label) + 1] > 2)
            FALSE else TRUE
}

unroot <- function(phy)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    N <- dim(phy$edge)[1]
    if (N < 3)
        stop("cannot unroot a tree with less than three edges.")

    ## delete FIRST the root.edge (in case this is sufficient to
    ## unroot the tree, i.e. there is a multichotomy at the root)
    if (!is.null(phy$root.edge)) phy$root.edge <- NULL
    if (!is.rooted(phy)) return(phy)

    n <- length(phy$tip.label)
    ROOT <- n + 1L

### EDGEROOT[1]: the edge to delete
### EDGEROOT[2]: the target where to stick the edge to delete

### If the tree is in pruningwise order, then the last two edges
### are those connected to the root; the node situated in
### phy$edge[N - 2L, 1L] will be the new root...

    ophy <- attr(phy, "order")
    if (!is.null(ophy) && ophy == "pruningwise") {
        NEWROOT <- phy$edge[N - 2L, 1L]
        EDGEROOT <- c(N, N - 1L)
    } else {

### ... otherwise, we remove one of the edges coming from
### the root, and eventually adding the branch length to
### the other one also coming from the root.
### In all cases, the node deleted is the 2nd one (numbered
### nb.tip+2 in 'edge'), so we simply need to renumber the
### nodes by adding 1, except the root (this remains the
### origin of the tree).

        EDGEROOT <- which(phy$edge[, 1L] == ROOT)
        NEWROOT <- ROOT + 1L
    }

    ## make sure EDGEROOT is ordered as described above:
    if (phy$edge[EDGEROOT[1L], 2L] != NEWROOT)
        EDGEROOT <- EDGEROOT[2:1]

    phy$edge <- phy$edge[-EDGEROOT[1L], ]

    s <- phy$edge == NEWROOT # renumber the new root
    phy$edge[s] <- ROOT

    s <- phy$edge > NEWROOT # renumber all nodes greater than the new root
    phy$edge[s] <- phy$edge[s] - 1L

    if (!is.null(phy$edge.length)) {
        phy$edge.length[EDGEROOT[2L]] <-
            phy$edge.length[EDGEROOT[2L]] + phy$edge.length[EDGEROOT[1L]]
        phy$edge.length <- phy$edge.length[-EDGEROOT[1L]]
    }

    phy$Nnode <- phy$Nnode - 1L

    if (!is.null(phy$node.label)) {
        if (NEWROOT == n + 2L)
            phy$node.label <- phy$node.label[-1]
        else {
            lbs <- phy$node.label
            tmp <- lbs[NEWROOT - n]
            lbs <- lbs[-c(1, NEWROOT)]
            phy$node.label <- c(tmp, lbs)
        }
    }
    phy
}

root <- function(phy, outgroup, node = NULL,
                 resolve.root = FALSE, interactive = FALSE)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    phy <- reorder(phy)
    n <- length(phy$tip.label)
    ROOT <- n + 1L
    if (interactive) {
        node <- identify(phy)$nodes
        cat("You have set resolve.root =", resolve.root, "\n")
    }
    if (!is.null(node)) {
        if (node <= n)
            stop("incorrect node#: should be greater than the number of taxa")
        outgroup <- NULL
        newroot <- node
    } else {
        if (is.numeric(outgroup)) {
            if (any(outgroup > n))
                stop("incorrect taxa#: should not be greater than the number of taxa")
            outgroup <- sort(outgroup) # used below
        }
        if (is.character(outgroup))
            outgroup <- which(phy$tip.label %in% outgroup)
        if (length(outgroup) == n) return(phy)

        ## First check that the outgroup is monophyletic--
        ## unless there's only one tip specified of course
        if (length(outgroup) > 1) {
            pp <- prop.part(phy)
            ingroup <- (1:n)[-outgroup]
            newroot <- 0L
            for (i in 2:phy$Nnode) {
                if (identical(pp[[i]], ingroup)) {
                    newroot <- i + n
                    break
                }
                if (identical(pp[[i]], outgroup)) {
                    newroot <- phy$edge[which(phy$edge[, 2] == i + n), 1]
                    break
                }
            }
            if (!newroot)
                stop("the specified outgroup is not monophyletic")
        } else newroot <- phy$edge[which(phy$edge[, 2] == outgroup), 1]
    }
    N <- Nedge(phy)
    oldNnode <- phy$Nnode
    if (newroot == ROOT) {
        if (resolve.root) {
            snw <- which(phy$edge[, 1] == newroot)
            if (length(snw) > 2) {
                a <- snw[1]:(snw[2] - 1)
                b <- snw[2]:N
                newnod <- oldNnode + n + 1
                phy$edge[snw[-1], 1] <- newnod
                phy$edge <- rbind(phy$edge[a, ], c(ROOT, newnod),
                                  phy$edge[b, ])
                if (!is.null(phy$edge.length))
                    phy$edge.length <-
                        c(phy$edge.length[a], 0, phy$edge.length[b])
                phy$Nnode <- phy$Nnode + 1L
                ## node renumbering (see comments below)
                newNb <- integer(n + oldNnode)
                newNb[newroot] <- n + 1L
                sndcol <- phy$edge[, 2] > n
                phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <-
                    (n + 2):(n + phy$Nnode)
                phy$edge[, 1] <- newNb[phy$edge[, 1]]
            }
        }
        return(phy)
    }

    phy$root.edge <- NULL # just in case...
    Nclade <- tabulate(phy$edge[, 1])[ROOT] # degree of the root node
    ## if only 2 edges connect to the root, we have to fuse them:
    fuseRoot <- Nclade == 2

    start <- which(phy$edge[, 1] == ROOT)
    end <- c(start[-1] - 1, N)
    o <- integer(N)
    INV <- logical(N)

    w <- which(phy$edge[, 2] == newroot)
    nod <- phy$edge[w, 1]
    i <- w
    NEXT <- 1L

    ## The loop below starts from the new root and goes up in the edge
    ## matrix reversing the edges that need to be as well as well
    ## inverting their order. The other edges must not be changed, so
    ## their indices are stored in `stack'.
    ## We then bind the other edges in a straightforward way.

    if (nod != ROOT) {
        ## it is important that the 3 next lines
        ## are inside this "if" statement
        o[NEXT] <- w
        NEXT <- NEXT + 1L
        INV[w] <- TRUE
        i <- w - 1L
        stack <- 0L
        repeat {
            if (phy$edge[i, 2] == nod) {
                if (stack) {
                    o[NEXT:(NEXT + stack - 1L)] <- (i + 1L):(i + stack)
                    NEXT <- NEXT + stack
                    stack <- 0L
                }
                if (phy$edge[i, 1] == ROOT) break
                o[NEXT] <- i
                NEXT <- NEXT + 1L
                INV[i] <- TRUE
                nod <- phy$edge[i, 1]
            } else stack <- stack + 1L
            i <- i - 1L
        }
    }

    ## we keep the edge leading to the old root if needed:
    if (!fuseRoot) {
        o[NEXT] <- i
        INV[i] <- TRUE
        NEXT <- NEXT + 1L
    }

    endOfOutgroup <- which(phy$edge[(w+1):N, 1] < newroot)[1] + w - 1
    if (is.na(endOfOutgroup)) endOfOutgroup <- N
    endOfClade <- end[end >= endOfOutgroup][1]

    ## bind the other clades...
    for (j in 1:Nclade) {
        if (end[j] == endOfClade) next
        ## do we have to fuse the two basal edges?
        if (fuseRoot) {
            phy$edge[start[j], 1] <- phy$edge[i, 2]
            if (!is.null(phy$edge.length))
                phy$edge.length[start[j]] <- phy$edge.length[start[j]] +
                    phy$edge.length[i]
        } #else {
          #  o[NEXT] <- i#start[j]
          #  NEXT <- NEXT + 1L
        #}
        s <- start[j]:end[j]
        ne <- length(s)
        o[NEXT:(NEXT + ne - 1L)] <- s
        NEXT <- NEXT + ne
    }

    ## possibly bind the edges below the outgroup till the end of the clade
    if (all(endOfOutgroup != end)) {
        j <- (endOfOutgroup + 1L):endOfClade
        ## we must take care that the branch inversions done above
        ## may have changed the hierarchy of clades here, so we
        ## travel from the bottom of this series of edges
        stack <- 0L
        inverted <- phy$edge[INV, 1] # <- fails if ', 2]' is used
        for (k in rev(j)) {
            if (any(phy$edge[k, 1] == inverted)) {
                o[NEXT] <- k
                NEXT <- NEXT + 1L
                if (stack){
                    o[NEXT:(NEXT + stack - 1L)] <- (k + 1L):(k + stack)
                    NEXT <- NEXT + stack
                    stack <- 0L
                }
            } else stack <- stack + 1L
        }
    }

    ## ... and the outgroup
    s <- (w + 1L):endOfOutgroup
    ne <- length(s)
    o[NEXT:(NEXT + ne - 1L)] <- s

    if (fuseRoot) {
        phy$Nnode <- oldNnode - 1L
        N <- N - 1L
    }
    phy$edge[INV, ] <- phy$edge[INV, 2:1]
    phy$edge <- phy$edge[o, ]
    if (!is.null(phy$edge.length))
        phy$edge.length <- phy$edge.length[o]

    if (resolve.root) {
        newnod <- oldNnode + n + 1
        if (length(outgroup) == 1L) {
            wh <- which(phy$edge[, 2] == outgroup)
            phy$edge[1] <- newnod
            phy$edge <-
                rbind(c(newroot, newnod), phy$edge[-wh, ], phy$edge[wh, ])
            snw <- which(phy$edge[, 1] == newroot)
            phy$edge[snw[length(snw) - 1], 1] <- newnod
            if (!is.null(phy$edge.length)) {
                phy$edge.length <-
                    c(0, phy$edge.length[-wh], phy$edge.length[wh])
            }
        } else {
            wh <- which(phy$edge[, 1] == newroot)
            phy$edge[wh[-1], 1] <- newnod
            s1 <- 1:(wh[2] - 1)
            s2 <- wh[2]:N
            phy$edge <-
                rbind(phy$edge[s1, ], c(newroot, newnod), phy$edge[s2, ])
            if (!is.null(phy$edge.length)) {
                tmp <- phy$edge.length[1]
                phy$edge.length[1] <- 0
                phy$edge.length <-
                    c(phy$edge.length[s1], tmp, phy$edge.length[s2])
            }
        }
        ## N <- N + 1L ... not needed
        phy$Nnode <- phy$Nnode + 1L
    }

    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format
    newNb <- integer(n + oldNnode)
    newNb[newroot] <- n + 1L
    sndcol <- phy$edge[, 2] > n
    ## executed from right to left, so newNb is modified before phy$edge:
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- n + 2:phy$Nnode
    phy$edge[, 1] <- newNb[phy$edge[, 1]]

    if (!is.null(phy$node.label)) {
        #browser()
        newNb <- newNb[-(1:n)]
        if (fuseRoot) {
            newNb <- newNb[-1]
            phy$node.label <- phy$node.label[-1]
        }
        phy$node.label <- phy$node.label[order(newNb)]
        if (resolve.root) {
            phy$node.label[is.na(phy$node.label)] <- phy$node.label[1]
            phy$node.label[1] <- NA
            ##phy$node.label <- c(phy$node.label[1], NA, phy$node.label[-1])
            ##phy$node.label <- c("NA", phy$node.label)
        }
    }
    phy
}

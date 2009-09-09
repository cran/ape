## drop.tip.R (2009-09-09)

##   Remove Tips in a Phylogenetic Tree

## Copyright 2003-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

extract.clade <- function(phy, node, root.edge = 0)
{
    Ntip <- length(phy$tip.label)
    ROOT <- Ntip + 1
    Nedge <- dim(phy$edge)[1]
    wbl <- !is.null(phy$edge.length)
    if (length(node) > 1) {
        node <- node[1]
        warning("only the first value of 'node' has been considered")
    }
    if (is.character(node)) {
        if (is.null(phy$node.label))
            stop("the tree has no node labels")
        node <- which(phy$node.label %in% node) + Ntip
    }
    if (node <= Ntip)
        stop("node number must be greater than the number of tips")
    if (node == ROOT) return(phy)
    phy <- reorder(phy) # insure it is in cladewise order
    root.node <- which(phy$edge[, 2] == node)
    start <- root.node + 1 # start of the clade looked for
    anc <- phy$edge[root.node, 1] # the ancestor of 'node'
    next.anc <- which(phy$edge[-(1:start), 1] <= anc) # find the next occurence of 'anc' or an 'older' node

    keep <- if (length(next.anc)) start + 0:(next.anc[1] - 1) else start:Nedge

    if (root.edge) {
        NewRootEdge <- phy$edge.length[root.node]
        root.edge <- root.edge - 1
        while (root.edge) {
            if (anc == ROOT) break
            i <- which(phy$edge[, 2] ==  anc)
            NewRootEdge <- NewRootEdge + phy$edge.length[i]
            root.edge <- root.edge - 1
            anc <- phy$edge[i, 1]
        }
        if (root.edge && !is.null(phy$root.edge))
            NewRootEdge <- NewRootEdge + phy$root.edge
        phy$root.edge <- NewRootEdge
    }

    phy$edge <- phy$edge[keep, ]
    if (wbl) phy$edge.length <- phy$edge.length[keep]
    TIPS <- phy$edge[, 2] <= Ntip
    tip <- phy$edge[TIPS, 2]
    phy$tip.label <- phy$tip.label[tip]
    ## keep the ordering so no need to reorder tip.label:
    phy$edge[TIPS, 2] <- order(tip)
    if (!is.null(phy$node.label))
        phy$node.label <- phy$node.label[sort(unique(phy$edge[, 1])) - Ntip]
    Ntip <- length(phy$tip.label)
    phy$Nnode <- dim(phy$edge)[1] - Ntip + 1L
    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format -- same as in root()
    newNb <- integer(Ntip + phy$Nnode)
    newNb[node] <- Ntip + 1L
    sndcol <- phy$edge[, 2] > Ntip
    ## executed from right to left, so newNb is modified before phy$edge:
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <-
        (Ntip + 2):(Ntip + phy$Nnode)
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    phy
}

drop.tip <-
    function(phy, tip, trim.internal = TRUE, subtree = FALSE,
             root.edge = 0, rooted = is.rooted(phy))
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')

    Ntip <- length(phy$tip.label)
    ## find the tips to drop:
    if (is.character(tip))
        tip <- which(phy$tip.label %in% tip)

    if (!rooted && subtree) {
        phy <- root(phy, (1:Ntip)[-tip][1])
        root.edge <- 0
    }

    phy <- reorder(phy)
    NEWROOT <- ROOT <- Ntip + 1
    Nnode <- phy$Nnode
    Nedge <- dim(phy$edge)[1]
    if (subtree) {
        trim.internal <- TRUE
        tr <- reorder(phy, "pruningwise")
        N <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
                as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]),
                as.integer(Nedge), double(Ntip + Nnode),
                DUP = FALSE, PACKAGE = "ape")[[6]]
    }
    wbl <- !is.null(phy$edge.length)
    edge1 <- phy$edge[, 1] # local copies
    edge2 <- phy$edge[, 2] #
    keep <- !logical(Nedge)

    ## find the tips to drop:
    if (is.character(tip))
        tip <- which(phy$tip.label %in% tip)

    if (!rooted && subtree) {
        phy <- root(phy, (1:Ntip)[-tip][1])
        root.edge <- 0
    }

    ## delete the terminal edges given by `tip':
    keep[match(tip, edge2)] <- FALSE

    if (trim.internal) {
        ints <- edge2 > Ntip
        ## delete the internal edges that do not have anymore
        ## descendants (ie, they are in the 2nd col of `edge' but
        ## not in the 1st one)
        repeat {
            sel <- !(edge2 %in% edge1[keep]) & ints & keep
            if (!sum(sel)) break
            keep[sel] <- FALSE
        }
        if (subtree) {
            ## keep the subtending edge(s):
            subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
            keep[subt] <- TRUE
        }
        if (root.edge && wbl) {
            degree <- tabulate(edge1[keep])
            if (degree[ROOT] == 1) {
                j <- integer(0) # will store the indices of the edges below the new root
                repeat {
                    i <- which(edge1 == NEWROOT & keep)
                    j <- c(i, j)
                    NEWROOT <- edge2[i]
                    degree <- tabulate(edge1[keep])
                    if (degree[NEWROOT] > 1) break
                }
                keep[j] <- FALSE
                if (length(j) > root.edge) j <- 1:root.edge
                NewRootEdge <- sum(phy$edge.length[j])
                if (length(j) < root.edge && !is.null(phy$root.edge))
                    NewRootEdge <- NewRootEdge + phy$root.edge
                phy$root.edge <- NewRootEdge
            }
        }
    }

    if (!root.edge) phy$root.edge <- NULL

    ## drop the edges
    phy$edge <- phy$edge[keep, ]
    if (wbl) phy$edge.length <- phy$edge.length[keep]

    ## find the new terminal edges (works whatever 'subtree' and 'trim.internal'):
    TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])

    ## get the old No. of the nodes and tips that become tips:
    oldNo.ofNewTips <- phy$edge[TERMS, 2]

    n <- length(oldNo.ofNewTips) # the new number of tips in the tree

    ## the tips may not be sorted in increasing order of their
    ## in the 2nd col of edge, so no need to reorder $tip.label
    phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])

    ## make new tip labels if necessary:
    if (subtree || !trim.internal) {
        ## get the logical indices of the tips that are kept within 'oldNo.ofNewTips':
        tips.kept <- oldNo.ofNewTips <= Ntip & !(oldNo.ofNewTips %in% tip)
        new.tip.label <- character(n)
        new.tip.label[tips.kept] <- phy$tip.label[-tip]
        ## get the numbers of the nodes that become tips:
        node2tip <- oldNo.ofNewTips[!tips.kept]
        new.tip.label[!tips.kept] <- if (subtree) {
            paste("[", N[node2tip], "_tips]", sep = "")
        } else {
            if (is.null(phy$node.label)) rep("NA", length(node2tip))
            else phy$node.label[node2tip - Ntip]
        }
        if (!is.null(phy$node.label))
            phy$node.label <- phy$node.label[-(node2tip - Ntip)]
        phy$tip.label <- new.tip.label
    } else phy$tip.label <- phy$tip.label[-tip]

    ## update node.label if needed:
    if (!is.null(phy$node.label))
        phy$node.label <- phy$node.label[sort(unique(phy$edge[, 1])) - Ntip]

    phy$Nnode <- dim(phy$edge)[1] - n + 1L # update phy$Nnode

    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format -- same as in root()
    newNb <- integer(n + phy$Nnode)
    newNb[NEWROOT] <- n + 1L
    sndcol <- phy$edge[, 2] > n
    ## executed from right to left, so newNb is modified before phy$edge:
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <-
        (n + 2):(n + phy$Nnode)
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    storage.mode(phy$edge) <- "integer"
    collapse.singles(phy)
}

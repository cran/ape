### bind.tree.R (2006-10-06)
###
###     Bind Trees
###
### Copyright 2003-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

bind.tree <- function(x, y, where = "root", position = 0)
{
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    ROOT <- nb.tip + 1
    if (where == 0 || where == "root")
      where <- ROOT
    if (position < 0) position <- 0
    if (where > nb.tip + nb.node) stop("node number out of range for tree 'x'")
    nb.edge <- dim(x$edge)[1]
    yHasNoRootEdge <- is.null(y$root.edge)
    xHasNoRootEdge <- is.null(x$root.edge)

    ## check whether both trees have branch lengths:
    wbl <- TRUE
    noblx <- is.null(x$edge.length)
    nobly <- is.null(y$edge.length)
    if (noblx && nobly) wbl <- FALSE
    if (xor(noblx, nobly)) {
        if (nobly) x$edge.length <- NULL
        else y$edge.length <- NULL
        wbl <- FALSE
        warning("one tree has no branch lengths, they will be ignored")
    }

    ## To avoid problems with tips or nodes with indentical
    ## labels we substitute the one where `y' is grafted:
    if (where <= nb.tip) {
        Tip.Label.where <- x$tip.label[where]
        x$tip.label[where] <- "TheTipWhereToGraftY"
    }
    if (where > ROOT) {
        xHasNoNodeLabel <- TRUE
        if (is.null(x$node.label)) {
            x$node.label <- paste("NODE", 1:nb.node, sep = "")
            x$node.label[where - nb.tip] <- "TheNodeWhereToGraftY"
        } else {
            Node.Label.where <- x$node.label[where - nb.tip]
            x$node.label[where - nb.tip] <- "TheNodeWhereToGraftY"
            xHasNoNodeLabel <- FALSE
        }
    }

    ## if we bind `y' under a node or tip of `y', we first
    ## adjust the edge lengths if needed
    if (position && wbl) {
        if (where == ROOT) {
            if (xHasNoRootEdge) stop("tree 'x' has no root edge")
            if (x$root.edge < position)
              stop("argument 'position' is larger than the root edge.")
            x$root.edge <- x$root.edge - position
        } else {
            i <- which(x$edge[, 2] == where)
            if (x$edge.length[i] < position)
              stop("argument 'position' is larger than the specified edge.")
            x$edge.length[i] <- x$edge.length[i] - position
        }
        if (yHasNoRootEdge ) y$root.edge <- position
        else y$root.edge <- y$root.edge + position
    }

    X <- write.tree(x, multi.line = FALSE)
    Y <- write.tree(y, multi.line = FALSE)
    Y <- substr(Y, 1, nchar(Y) - 1)

    if (where <= nb.tip) {
        if (position)
          X <- gsub("TheTipWhereToGraftY",
                    paste("(", "TheTipWhereToGraftY", ",", Y, ")",
                          sep = ""), X)
        else
          X <- gsub("TheTipWhereToGraftY", Y, X)
    }
    if (where == ROOT) {
        rmvx <- if (xHasNoRootEdge) "\\);$" else ";$"
        X <- gsub(rmvx, "", X)
        Y <- gsub("^\\(", "", Y)
        if (!xHasNoRootEdge) X <- paste("(", X, sep = "")
        X <- paste(X, ",", Y, ";", sep = "")
    }
    if (where > ROOT) {
        if (position) {
            ## find where is the node in `X':
            ## below 19 is: nchar("TheNodeWhereToGraftY") - 1
            for (i in 1:nchar(X)) {
                if ("TheNodeWhereToGraftY" == substr(X, i, i + 19))
                  break
                i <- i + 1
            }
            ## now go back to find the left matching parentheses
            n.paren <- 1
            i <- i - 2
            while (n.paren > 0) {
                if (substr(X, i, i) == ")") n.paren <- n.paren + 1
                if (substr(X, i, i) == "(") n.paren <- n.paren - 1
                i <- i - 1
            }
            ## insert the left parenthesis:
            ## here 21 is: nchar("TheNodeWhereToGraftY") + 1
            X <- paste(substr(X, 1, i - 1), "(",
                       substr(X, i, 21), sep = "")
            ## and insert `y':
            X <- gsub("TheNodeWhereToGraftY",
                      paste("TheNodeWhereToGraftY", ",", Y,
                            sep = ""), X)
        } else {
            xx <- paste(")", "TheNodeWhereToGraftY", sep = "")
            X <- gsub(xx, paste(",", Y, xx, sep = ""), X)
        }
    }
    phy <- read.tree(text = X)
    ## restore the labels:
    if (where <= nb.tip)
      phy$tip.label[which(phy$tip.label == "TheTipWhereToGraftY")] <-
        Tip.Label.where
    if (where > ROOT) {
        if (xHasNoNodeLabel) phy$node.label <- NULL
        else
          phy$node.label[which(phy$node.label == "TheNodeWhereToGraftY")] <-
            Node.Label.where
    }
    phy
}

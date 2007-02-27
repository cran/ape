### dist.topo.R (2007-02-27)
###
###      Topological Distances, Tree Bipartitions,
###   Consensus Trees, and Bootstrapping Phylogenies
###
### Copyright 2005-2007 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

dist.topo <- function(x, y, method = "PH85")
{
    if (method == "BHV01" && (is.null(x$edge.length) || is.null(y$edge.length)))
      stop("trees must have branch lengths for Billera et al.'s distance.")
    n <- length(x$tip.label)
    bp1 <- .Call("bipartition", x$edge[, 1], x$edge[, 2],
                  n, x$Nnode, PACKAGE = "ape")
    bp1 <- lapply(bp1, function(xx) sort(x$tip.label[xx]))
    bp2 <- .Call("bipartition", y$edge[, 1], y$edge[, 2],
                 n, y$Nnode, PACKAGE = "ape")
    bp2 <- lapply(bp2, function(xx) sort(y$tip.label[xx]))
    q1 <- length(bp1)
    q2 <- length(bp2)
    if (method == "PH85") {
        p <- 0
        for (i in 1:q1) {
            for (j in 1:q2) {
                if (identical(all.equal(bp1[[i]], bp2[[j]]), TRUE)) {
                    p <- p + 1
                    break
                }
            }
        }
        dT <- if (q1 == q2) 2*(q1 - p) else 2*(min(q1, q2) - p) + abs(q1 - q2)
    }
    if (method == "BHV01") {
        dT <- 0
        found1 <- FALSE
        found2 <- logical(q2)
        found2[1] <- TRUE
        for (i in 2:q1) {
            for (j in 2:q2) {
                if (identical(bp1[[i]], bp2[[j]])) {
                    dT <- dT + abs(x$edge.length[which(x$edge[, 2] == n + i)] -
                                   y$edge.length[which(y$edge[, 2] == n + j)])
                    found1 <- found2[j] <- TRUE
                    break
                }
            }
            if (found1) found1 <- FALSE
            else dT <- dT + x$edge.length[which(x$edge[, 2] == n + i)]
        }
        if (!all(found2))
          dT <- dT + sum(y$edge.length[y$edge[, 2] %in% (n + which(!found2))])
    }
    dT
}

prop.part <- function(...)
{
    obj <- list(...)
    if (length(obj) == 1 && class(obj[[1]]) != "phylo")
      obj <- unlist(obj, recursive = FALSE)
    ntree <- length(obj)
    bp <- .Call("bipartition", obj[[1]]$edge[, 1], obj[[1]]$edge[, 2],
                length(obj[[1]]$tip.label), obj[[1]]$Nnode,
                PACKAGE = "ape")
    clades <- lapply(bp, function(xx) sort(obj[[1]]$tip.label[xx]))
    no <- rep(1, length(clades))

    if (ntree > 1) {
        for (k in 2:ntree) {
            bp <- .Call("bipartition", obj[[k]]$edge[, 1], obj[[k]]$edge[, 2],
                        length(obj[[k]]$tip.label), obj[[k]]$Nnode,
                        PACKAGE = "ape")
            bp <- lapply(bp, function(xx) sort(obj[[k]]$tip.label[xx]))
            for (i in 1:length(bp)) {
                done <- FALSE
                for (j in 1:length(clades)) {
                    if (identical(all.equal(bp[[i]], clades[[j]]), TRUE)) {
                        no[j] <- no[j] + 1
                        done <- TRUE
                        break
                    }
                }
                if (!done) {
                    clades <- c(clades, bp[i])
                    no <- c(no, 1)
                }
            }
        }
    }
    attr(clades, "number") <- no
    class(clades) <- "prop.part"
    clades
}

print.prop.part <- function(x, ...)
{
    for (i in 1:length(x)) {
        cat("==>", attr(x, "number")[i], "time(s):")
        print(x[[i]], quote = FALSE)
    }
}

summary.prop.part <- function(object, ...) attr(object, "number")

prop.clades <- function(phy, ..., part = NULL)
{
    if (is.null(part)) {
        obj <- list(...)
        if (length(obj) == 1 && class(obj[[1]]) != "phylo")
          obj <- unlist(obj, recursive = FALSE)
        part <- prop.part(obj)
    }
    bp <- .Call("bipartition", phy$edge[, 1], phy$edge[, 2],
                length(phy$tip.label), phy$Nnode,
                PACKAGE = "ape")
    bp <- lapply(bp, function(xx) sort(phy$tip.label[xx]))
    n <- numeric(length(bp))
    for (i in 1:length(bp)) {
        for (j in 1:length(part)) {
            if (identical(all.equal(bp[[i]], part[[j]]), TRUE)) {
                n[i] <- attr(part, "number")[j]
                done <-  TRUE
                break
            }
        }
    }
    n
}

boot.phylo <- function(phy, x, FUN, B = 100, block = 1)
{
    if (is.list(x)) {
        nm <- names(x)
        n <- length(x)
        x <- unlist(x)
        nL <- length(x)
        x <- matrix(x, n, nL/n, byrow = TRUE)
        rownames(x) <- nm
    }
    boot.tree <- list()
    length(boot.tree) <- B
    for (i in 1:B) {
        if (block > 1) {
            y <- seq(block, ncol(x), block)
            boot.i <- sample(y, replace = TRUE)
            boot.samp <- numeric(ncol(x))
            boot.samp[y] <- boot.i
            for (j in 1:(block - 1))
              boot.samp[y - j] <- boot.i - j
        } else boot.samp <- sample(ncol(x), replace = TRUE)
        boot.tree[[i]] <- FUN(x[, boot.samp])
    }
    prop.clades(phy, boot.tree)
}

consensus <- function(..., p = 1)
{
    obj <- list(...)
    if (length(obj) == 1 && class(obj[[1]]) != "phylo")
      obj <- unlist(obj, recursive = FALSE)
    ntree <- length(obj)
    ## Get all observed partitions and their frequencies:
    pp <- prop.part(obj)
    ## Drop the partitions whose frequency is less than 'p':
    pp <- pp[attr(pp, "number") >= p * ntree]
    ## Get the order of the remaining partitions by decreasing size:
    ind <- rev(sort(unlist(lapply(pp, length)),
                    index.return = TRUE)$ix)
    STRING <- paste(pp[[1]], collapse = ",")
    STRING <- paste("(", STRING, ");", sep = "")
    for (i in ind[-1]) {
        ## 1. Delete all tips in the focus partition:
        STRING <- unlist(strsplit(STRING, paste(pp[[i]], collapse = "|")))
        ## 2. Put the partition in any of the created gaps:
        STRING <- c(STRING[1],
                    paste("(", paste(pp[[i]], collapse = ","), ")", sep = ""),
                    STRING[-1])
        ## 3. Stick back the Newick string:
        STRING <- paste(STRING, collapse = "")
    }
    ## Remove the extra commas:
    STRING <- gsub(",{2,}", ",", STRING)
    STRING <- gsub("\\(,", "\\(", STRING)
    STRING <- gsub(",\\)", "\\)", STRING)
    read.tree(text = STRING)
}

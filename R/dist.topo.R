## dist.topo.R (2007-07-04)

##      Topological Distances, Tree Bipartitions,
##   Consensus Trees, and Bootstrapping Phylogenies

## Copyright 2005-2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

dist.topo <- function(x, y, method = "PH85")
{
    if (method == "BHV01" && (is.null(x$edge.length) || is.null(y$edge.length)))
      stop("trees must have branch lengths for Billera et al.'s distance.")
    n <- length(x$tip.label)
    bp1 <- .Call("bipartition", x$edge, n, x$Nnode, PACKAGE = "ape")
    bp1 <- lapply(bp1, function(xx) sort(x$tip.label[xx]))
    bp2 <- .Call("bipartition", y$edge, n, y$Nnode, PACKAGE = "ape")
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

prop.part <- function(..., check.labels = FALSE)
{
    obj <- list(...)
    if (length(obj) == 1 && class(obj[[1]]) != "phylo")
      obj <- unlist(obj, recursive = FALSE)
    ntree <- length(obj)
    if (!check.labels) {
        for (i in 1:ntree) storage.mode(obj[[i]]$Nnode) <- "integer"
        clades <- .Call("prop_part", obj, ntree, TRUE, PACKAGE = "ape")
        attr(clades, "number") <- attr(clades, "number")[1:length(clades)]
        attr(clades, "labels") <- obj[[1]]$tip.label
    } else {
        bp <- .Call("bipartition", obj[[1]]$edge, length(obj[[1]]$tip.label),
                    obj[[1]]$Nnode, PACKAGE = "ape")
        clades <- lapply(bp, function(xx) sort(obj[[1]]$tip.label[xx]))
        no <- rep(1, length(clades))

        if (ntree > 1) {
            for (k in 2:ntree) {
                bp <- .Call("bipartition", obj[[k]]$edge,
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
    }
    class(clades) <- "prop.part"
    clades
}

print.prop.part <- function(x, ...)
{
    if (is.null(attr(x, "labels"))) {
        for (i in 1:length(x)) {
            cat("==>", attr(x, "number")[i], "time(s):")
            print(x[[i]], quote = FALSE)
        }
    } else {
        for (i in 1:length(attr(x, "labels")))
          cat(i, ": ", attr(x, "labels")[i], "\n", sep = "")
        cat("\n")
        for (i in 1:length(x)) {
            cat("==>", attr(x, "number")[i], "time(s):")
            print(x[[i]], quote = FALSE)
        }
    }
}

summary.prop.part <- function(object, ...) attr(object, "number")

plot.prop.part <- function(x, barcol = "blue", leftmar = 4, ...)
{
    if (is.null(attr(x, "labels")))
      stop("cannot plot this partition object; see ?prop.part for details.")
    L <- length(x)
    n <- length(attr(x, "labels"))
    layout(matrix(1:2, 2, 1), heights = c(1, 3))
    par(mar = c(0.1, leftmar, 0.1, 0.1))
    plot(1:L, attr(x, "number"), type = "h", col = barcol, xlim = c(1, L),
         xlab = "", ylab = "Number", xaxt = "n", bty = "n")
    plot(0, type = "n", xlim = c(1, L), ylim = c(1, n),
         xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    for (i in 1:L) points(rep(i, length(x[[i]])), x[[i]], ...)
    mtext(attr(x, "labels"), side = 2, at = 1:n, las = 1)
}

prop.clades <- function(phy, ..., part = NULL)
{
    if (is.null(part)) {
        obj <- list(...)
        if (length(obj) == 1 && class(obj[[1]]) != "phylo")
          obj <- unlist(obj, recursive = FALSE)
        part <- prop.part(obj, check.labels = TRUE)
    }
    bp <- .Call("bipartition", phy$edge, length(phy$tip.label),
                phy$Nnode, PACKAGE = "ape")
    if (!is.null(attr(part, "labels")))
      for (i in 1:length(part))
        part[[i]] <- sort(attr(part, "labels")[part[[i]]])
    bp <- lapply(bp, function(xx) sort(phy$tip.label[xx]))
    n <- numeric(phy$Nnode)
    for (i in 1:phy$Nnode) {
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
    boot.tree <- vector("list", B)
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
    for (i in 1:B) storage.mode(boot.tree[[i]]$Nnode) <- "integer"
    storage.mode(phy$Nnode) <- "integer"
    attr(.Call("prop_part", c(list(phy), boot.tree), B + 1, FALSE,
               PACKAGE = "ape"), "number") - 1
}

consensus <- function(..., p = 1)
{
    obj <- list(...)
    if (length(obj) == 1 && class(obj[[1]]) != "phylo")
      obj <- unlist(obj, recursive = FALSE)
    ntree <- length(obj)
    ## Get all observed partitions and their frequencies:
    pp <- prop.part(obj, check.labels = TRUE)
    ## Drop the partitions whose frequency is less than 'p':
    pp <- pp[attr(pp, "number") >= p * ntree]
    ## Get the order of the remaining partitions by decreasing size:
    ind <- rev(sort(unlist(lapply(pp, length)),
                    index.return = TRUE)$ix)
    pp <- lapply(pp, function(xx) paste("IMPROBABLE_PREFIX", xx,
                                        "IMPROBABLE_SUFFIX", sep = "_"))
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
    STRING <- gsub("IMPROBABLE_PREFIX_", "", STRING)
    STRING <- gsub("_IMPROBABLE_SUFFIX", "", STRING)
    read.tree(text = STRING)
}

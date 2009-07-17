## dist.topo.R (2009-07-06)

##      Topological Distances, Tree Bipartitions,
##   Consensus Trees, and Bootstrapping Phylogenies

## Copyright 2005-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

dist.topo <- function(x, y, method = "PH85")
{
    if (method == "BHV01" && (is.null(x$edge.length) || is.null(y$edge.length)))
      stop("trees must have branch lengths for Billera et al.'s distance.")
    n <- length(x$tip.label)
    bp1 <- .Call("bipartition", x$edge, n, x$Nnode, PACKAGE = "ape")
    bp1 <- lapply(bp1, function(xx) sort(x$tip.label[xx]))
    ## fix by Tim Wallstrom:
    bp2.tmp <- .Call("bipartition", y$edge, n, y$Nnode, PACKAGE = "ape")
    bp2 <- lapply(bp2.tmp, function(xx) sort(y$tip.label[xx]))
    bp2.comp <- lapply(bp2.tmp, function(xx) setdiff(1:n, xx))
    bp2.comp <- lapply(bp2.comp, function(xx) sort(y$tip.label[xx]))
    ## End
    q1 <- length(bp1)
    q2 <- length(bp2)
    if (method == "PH85") {
        p <- 0
        for (i in 1:q1) {
            for (j in 1:q2) {
                if (identical(bp1[[i]], bp2[[j]]) |
                    identical(bp1[[i]], bp2.comp[[j]])) {
                    p <- p + 1
                    break
                }
            }
        }
        dT <- q1 + q2 - 2 * p # same than:
        ##dT <- if (q1 == q2) 2*(q1 - p) else 2*(min(q1, q2) - p) + abs(q1 - q2)
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

.compressTipLabel <- function(x)
{
    ## 'x' is a list of objects of class "phylo" possibly with no class
    if (!is.null(attr(x, "TipLabel"))) return(x)
    ref <- x[[1]]$tip.label
    if (any(table(ref) != 1))
        stop("some tip labels are duplicated in tree no. 1")
    n <- length(ref)
    for (i in 2:length(x)) {
        if (identical(x[[i]]$tip.label, ref)) next
        ilab <- match(x[[i]]$tip.label, ref)
        ## can use tabulate here because 'ilab' contains integers
        if (any(tabulate(ilab) > 1))
            stop(paste("some tip labels are duplicated in tree no.", i))
        if (any(is.na(ilab)))
            stop(paste("tree no.", i, "has different tip labels"))
        ie <- match(1:n, x[[i]]$edge[, 2])
        x[[i]]$edge[ie, 2] <- ilab
    }
    for (i in 1:length(x)) x[[i]]$tip.label <- NULL
    attr(x, "TipLabel") <- ref
    x
}

prop.part <- function(..., check.labels = TRUE)
{
    obj <- list(...)
    if (length(obj) == 1 && class(obj[[1]]) != "phylo")
        obj <- obj[[1]]
    ## <FIXME>
    ## class(obj) <- NULL # needed?
    ## </FIXME>
    ntree <- length(obj)
    if (ntree == 1) check.labels <- FALSE
    if (check.labels) obj <- .compressTipLabel(obj)
    for (i in 1:ntree) storage.mode(obj[[i]]$Nnode) <- "integer"
    ## <FIXME>
    ## The 1st must have tip labels
    ## Maybe simply pass the number of tips to the C code??
    if (!is.null(attr(obj, "TipLabel")))
        for (i in 1:ntree) obj[[i]]$tip.label <- attr(obj, "TipLabel")
    ## </FIXME>
    clades <- .Call("prop_part", obj, ntree, TRUE, PACKAGE = "ape")
    attr(clades, "number") <- attr(clades, "number")[1:length(clades)]
    attr(clades, "labels") <- obj[[1]]$tip.label
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
         xlab = "", ylab = "Frequency", xaxt = "n", bty = "n")
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

boot.phylo <- function(phy, x, FUN, B = 100, block = 1, trees = FALSE)
{
    if (is.list(x) && !is.data.frame(x)) {
        if (inherits(x, "DNAbin")) x <- as.matrix(x)
        else {
            nm <- names(x)
            n <- length(x)
            x <- unlist(x)
            nL <- length(x)
            x <- matrix(x, n, nL/n, byrow = TRUE)
            rownames(x) <- nm
        }
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
    ans <- attr(.Call("prop_part", c(list(phy), boot.tree),
                      B + 1, FALSE, PACKAGE = "ape"), "number") - 1
    if (trees) ans <- list(BP = ans, trees = boot.tree)
    ans
}

consensus <- function(..., p = 1, check.labels = TRUE)
{
    foo <- function(ic, node) {
        ## ic: index of 'pp'
        ## node: node number in the final tree
        pool <- pp[[ic]]
        if (ic < m) {
            for (j in (ic + 1):m) {
                wh <- match(pp[[j]], pool)
                if (!any(is.na(wh))) {
                    edge[pos, 1] <<- node
                    pool <- pool[-wh]
                    edge[pos, 2] <<- nextnode <<- nextnode + 1L
                    pos <<- pos + 1L
                    foo(j, nextnode)
                }
            }
        }
        size <- length(pool)
        if (size) {
            ind <- pos:(pos + size - 1)
            edge[ind, 1] <<- node
            edge[ind, 2] <<- pool
            pos <<- pos + size
        }
    }
    obj <- list(...)
    if (length(obj) == 1) {
        ## better than unlist(obj, recursive = FALSE)
        ## because "[[" keeps the class of 'obj':
        obj <- obj[[1]]
        if (class(obj) == "phylo") return(obj)
    }
    if (!is.null(attr(obj, "TipLabel")))
        labels <- attr(obj, "TipLabel")
    else {
        labels <- obj[[1]]$tip.label
        if (check.labels) obj <- .compressTipLabel(obj)
    }
    ntree <- length(obj)
    ## Get all observed partitions and their frequencies:
    pp <- prop.part(obj, check.labels = FALSE)
    ## Drop the partitions whose frequency is less than 'p':
    pp <- pp[attr(pp, "number") >= p * ntree]
    ## Get the order of the remaining partitions by decreasing size:
    ind <- sort(unlist(lapply(pp, length)), decreasing = TRUE,
                index.return = TRUE)$ix
    pp <- pp[ind]
    n <- length(labels)
    m <- length(pp)
    edge <- matrix(0L, n + m - 1, 2)
    if (m == 1) {
        edge[, 1] <- n + 1L
        edge[, 2] <- 1:n
    } else {
        nextnode <- n + 1L
        pos <- 1L
        foo(1, nextnode)
    }
    structure(list(edge = edge, tip.label = labels,
              Nnode = m), class = "phylo")
}

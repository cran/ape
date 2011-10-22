## ltt.plot.R (2011-09-24)

##    Lineages Through Time Plot

## Copyright 2002-2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

ltt.plot.coords <- function(phy, backward = TRUE, tol = 1e-6)
{
    if (is.ultrametric(phy, tol)) {
        if (is.binary.tree(phy)) {
            N <- numeric(phy$Nnode + 1)
            N[] <- 1
        } else {
            node.order <- tabulate(phy$edge[, 1])
            N <- node.order[-(1:length(phy$tip.label))] - 1
        }
        bt <- branching.times(phy)
        names(bt) <- NULL
        o <- order(bt, decreasing = TRUE)
        time <- c(-bt[o], 0)
        if (!is.binary.tree(phy)) N <- c(1, N[o])
    } else {
        if (!is.binary.tree(phy)) phy <- multi2di(phy)
        n <- Ntip(phy)
        m <- phy$Nnode
        ROOT <- n + 1L
        event <- time.event <- numeric(n + m)

        time.event[ROOT] <- 0
        phy <- reorder(phy)

        for (i in 1:nrow(phy$edge))
            time.event[phy$edge[i, 2]] <- time.event[phy$edge[i, 1]] + phy$edge.length[i]

        present <- max(time.event)
        event[1:n] <- -1
        event[ROOT:(n + m)] <- 1

        ## delete the events that are too close to present:
        past.event <- present - time.event > tol
        event <- event[past.event]
        time.event <- time.event[past.event]

        ## reorder wrt time:
        o <- order(time.event)
        time.event <- time.event[o]
        event <- event[o]

        time <- c(time.event - present, 0)
        N <- c(1, event)
    }
    N <- cumsum(N)
    if (!is.null(phy$root.edge)) {
        time <- c(time[1] - phy$root.edge, time)
        N <- c(1, N)
    }
    if (!backward) time <- time - time[1]
    cbind(time, N)
}

ltt.plot <- function(phy, xlab = "Time", ylab = "N",
                     backward = TRUE, tol = 1e-6, ...)
{
    if (!inherits(phy, "phylo"))
        stop("object \"phy\" is not of class \"phylo\"")

    xy <- ltt.plot.coords(phy, backward, tol)

    plot.default(xy, xlab = xlab, ylab = ylab, xaxs = "r",
                 yaxs = "r", type = "S", ...)
}

ltt.lines <- function(phy, backward = TRUE, tol = 1e-6, ...)
{
    xy <- ltt.plot.coords(phy, backward, tol)
    lines(xy, type = "S", ...)
}

mltt.plot <-
    function(phy, ..., dcol = TRUE, dlty = FALSE, legend = TRUE,
             xlab = "Time", ylab = "N", log = "", backward = TRUE, tol = 1e-6)
{
    if (inherits(phy, "phylo")) { # if a tree of class "phylo"
        TREES <- list(ltt.plot.coords(phy, backward, tol))
        names(TREES) <- deparse(substitute(phy))
    } else { # a list of trees
        TREES <- lapply(phy, ltt.plot.coords, backward = backward, tol = tol)
        names(TREES) <- names(phy)
        if (is.null(names(TREES)))
            names(TREES) <-
                paste(deparse(substitute(phy)), "-", 1:length(TREES))
    }
    dts <- list(...)
    n <- length(dts)
    if (n) {
        mc <- as.character(match.call())[-(1:2)]
        nms <- mc[1:n]
        for (i in 1:n) {
            if (class(dts[[i]]) == "phylo") {
                a <- list(ltt.plot.coords(dts[[i]], backward, tol))
                names(a) <- nms[i]
            } else { # a list of trees
                a <- lapply(dts[[i]], ltt.plot.coords, backward = backward, tol = tol)
                names(a) <- names(dts[[i]])
                if (is.null(names(a)))
                    names(a) <- paste(deparse(substitute(phy)), "-", 1:length(a))
            }
            TREES <- c(TREES, a)
        }
    }
    n <- length(TREES)
    xl <- c(min(unlist(lapply(TREES, function(x) min(x[, 1])))), 0)
    yl <- c(1, max(unlist(lapply(TREES, function(x) max(x[, 2])))))

    plot.default(NA, type = "n", xlim = xl, ylim = yl, xaxs = "r",
                 yaxs = "r", xlab = xlab, ylab = ylab, log = log)

    lty <- if (!dlty) rep(1, n) else 1:n
    col <- if (!dcol) rep(1, n) else topo.colors(n)

    for (i in 1:n)
        lines(TREES[[i]], col = col[i], lty = lty[i], type = "S")

    if (legend)
        legend(xl[1], yl[2], legend = names(TREES),
               lty = lty, col = col, bty = "n")
}

ltt.coplot <- function(phy, backward = TRUE, ...)
{
    layout(matrix(1:2, 2))
    par(mar = c(0, 3, 0.5, 0.5))
    o <- plot(phy, root.edge = TRUE, ...)
    par(mar = c(3, 3, 0, 0.5))
    ltt.plot(phy, xlim = o$x.lim, backward = FALSE, xaxt = "n")
    if (backward) axisPhylo() else axis(1)
}

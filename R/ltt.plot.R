## ltt.plot.R (2008-02-22)

##    Lineages Through Time Plot

## Copyright 2002-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

ltt.plot <- function(phy, xlab = "Time", ylab = "N", ...)
{
    if (class(phy) != "phylo") stop('object "phy" is not of class "phylo"')
    time <- sort(branching.times(phy), decreasing = TRUE)
    N <- 1:(length(time) + 1)
    plot.default(-c(time, 0), N, xlab = xlab, ylab = ylab,
                 xaxs = "r", yaxs = "r", type = "S", ...)
}

ltt.lines <- function(phy, ...)
{
    time <- sort(branching.times(phy), decreasing = TRUE)
    N <- 1:(length(time) + 1)
    lines(-c(time, 0), N, type = "S", ...)
}

mltt.plot <- function(phy, ..., dcol = TRUE, dlty = FALSE, legend = TRUE,
                      xlab = "Time", ylab = "N", log = "")
{
    ltt.xy <- function(phy) {
        x <- -c(sort(branching.times(phy), decreasing = TRUE), 0)
        names(x) <- NULL
        y <- 1:length(x)
        cbind(x, y)
    }
    if (class(phy) == "phylo") {
        TREES <- list(ltt.xy(phy))
        names(TREES) <- deparse(substitute(phy))
    } else { # a list of trees
        TREES <- lapply(phy, ltt.xy)
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
                a <- list(ltt.xy(dts[[i]]))
                names(a) <- nms[i]
            } else { # a list of trees
                a <- lapply(dts[[i]], ltt.xy)
                names(a) <- names(dts[[i]])
                if (is.null(names(a)))
                    names(a) <-
                        paste(deparse(substitute(phy)), "-", 1:length(a))
            }
            TREES <- c(TREES, a)
        }
    }
    n <- length(TREES)
    xl <- c(min(unlist(lapply(TREES, function(x) min(x[, 1])))), 0)
    yl <- c(1, max(unlist(lapply(TREES, function(x) max(x[, 2])))))

    plot.default(1, 1, type = "n", xlim = xl, ylim = yl, xaxs = "r",
                 yaxs = "r", xlab = xlab, ylab = ylab, log = log)

    lty <- if (!dlty) rep(1, n) else 1:n
    col <- if (!dcol) rep(1, n) else topo.colors(n)

    for (i in 1:n)
      lines(TREES[[i]], col = col[i], lty = lty[i], type = "S")

    if (legend)
      legend(xl[1], yl[2], legend = names(TREES),
             lty = lty, col = col, bty = "n")
}

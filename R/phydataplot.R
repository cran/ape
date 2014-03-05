## phydataplot.R (2014-01-20)

##   Annotate Phylogenies

## Copyright 2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

polar2rect <- function(r, angle)
    list(x = r * cos(angle), y = r * sin(angle))

rect2polar <- function(x, y)
    list(r = sqrt(x^2 + y^2), angle = atan2(y, x))

.matchDataPhylo <- function(x, phy)
{
    msg <- "'x' has no (row)names: data are assumed to be in the same order than the tips of the tree"
    labs <- phy$tip.label
    if (is.vector(x)) { # also for lists
        if (is.null(names(x))) warning(msg) else x <- x[labs]
    } else {
        if (is.null(rownames(x))) warning(msg) else x <- x[labs, ]
    }
    x
}

ring <- function(x, phy, style = "ring", offset = 1, ...)
{
    style <- match.arg(style, c("ring", "segments", "arrows"))
    x <- .matchDataPhylo(x, phy)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    n <- lastPP$Ntip
    one2n <- seq_len(n)
    tmp <- rect2polar(lastPP$xx[one2n], lastPP$yy[one2n])
    theta <- tmp$angle
    r0 <- max(tmp$r) + offset
    r1 <- r0 + x
    s0 <- polar2rect(rep.int(r0, 100L), seq(0, 2*pi, length.out = 100L))
    s1 <- polar2rect(r1, theta)

    switch(style, ring = {
        if (length(x) < n) x <- rep_len(x, n)
        dx <- dim(x)
        if (is.null(dx)) dim(x) <- dx <- c(n, 1L)
        nc <- dx[2]
        col <- list(...)$col
        if (is.null(col)) col <- "grey"
        if (nc == 1) {
            col <- rep_len(col, n)
        } else {
            colvar <- col
            col <- rep(col[1], n)
        }
        iangle <- min(diff(sort(theta)))
        iangle2 <- iangle / 2
        for (i in one2n) {
            R <- rep(r0, 100)
            THETA <- seq(theta[i] - iangle2, theta[i] + iangle2, length.out = 100)
            xy1 <- polar2rect(R, THETA)
            xy2 <- polar2rect(R + x[i, 1], THETA)
            polygon(c(xy1$x, rev(xy2$x)), c(xy1$y, rev(xy2$y)), col = col[i], border = NA)
            if (nc > 1) {
                for (j in 2:nc) {
                    xy1 <- xy2
                    xy2 <- polar2rect(R + sum(x[i, 1:j]), THETA)
                    polygon(c(xy1$x, rev(xy2$x)), c(xy1$y, rev(xy2$y)), col = colvar[j], border = NA)
                }
            }
        }
        ##polypath(c(s0$x, NA, s0$x), c(s0$y, NA, s1$y), rule = "evenodd",
        ##         border = 1, col = "transparent")
    }, segments = {
        s0 <- polar2rect(rep.int(r0, n), theta)
        segments(s0$x, s0$y, s1$x, s1$y, ...)
    },  arrows = {
        s0 <- polar2rect(rep.int(r0, n), theta)
        fancyarrows(s0$x, s0$y, s1$x, s1$y, ...)
    })
}

phydataplot <- function(x, phy, style = "bars", offset = 1, ...)
{
    style <- match.arg(style, c("bars", "segments", "image", "arrows"))
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    circular <- if (lastPP$type %in% c("radial", "fan")) TRUE else FALSE

    n <- length(phy$tip.label)
    one2n <- seq_len(n)
    x <- .matchDataPhylo(x, phy)

    if (!circular) {
        if (lastPP$direction != "rightwards")
            stop("for the moment, only rightwards trees are supported")
        x0 <- max(lastPP$xx[one2n]) + offset
        if (style != "image") x1 <- x0 + x
        y1 <- lastPP$yy[one2n]
        if (style %in% c("bars", "image")) {
            o <- order(y1)
            x <- if (style == "image") x[o, o] else
            if (is.vector(x)) x[o] else x[o, ]
        }
    }

    switch(style, bars = {
        if (circular)
            stop("style = \"bars\" not implemented with circular trees; see the function 'ring'")
        if (!is.null(dim(x))) x <- t(x)
        barplot(x, width = 1, add = TRUE, horiz = TRUE, offset = x0,
                axes = FALSE, axisnames = FALSE, space = c(0.5, rep(0, n - 1)), ...)
        px <- pretty(c(0, x))
        axis(1, px + x0, labels = px, line = 1)
    }, segments = {
        if (circular) ring(x, phy, style, offset, ...)
        else segments(x0, y1, x1, y1, ...)
    }, image = {
        if (circular)
            stop("style = \"image\" not implemented with circular trees")
        if (inherits(x, "DNAbin"))
            stop("object of class \"DNAbin\" not yet supported")
        ##image(x, show.labels = FALSE, add = TRUE, ...)
        x1 <- seq(x0, lastPP$x.lim[2], length.out = n)
        image(x1, y1[o], x, add = TRUE, ...)
        mtext(phy$tip.label[o], 1, 1, at = x1, font = lastPP$font,
              cex = lastPP$cex, col = lastPP$tip.color)
    }, arrows = {
        if (circular) ring(x, phy, style, offset, ...)
        else fancyarrows(x0, y1, x1, y1, ...)
    })
}

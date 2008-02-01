## nodelabels.R (2007-03-05)

##   Labelling Trees

## Copyright 2004-2007 Emmanuel Paradis, 2006 Ben Bolker, and 2006 Jim Lemon

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

## from JL:
## floating.pie() from plotrix with two changes:
## (1) aspect ratio fixed, so pies will appear circular
##     (`radius' is the radius in user coordinates along the x axis);
## (2) zero values allowed (but not negative).

floating.pie.asp <- function(xpos, ypos, x, edges = 200, radius = 1,
                             col = NULL, startpos = 0, ...)
{
    u <- par("usr")
    user.asp <- diff(u[3:4])/diff(u[1:2])
    p <- par("pin")
    inches.asp <- p[2]/p[1]
    asp <- user.asp/inches.asp
    if (!is.numeric(x) || any(is.na(x) | x < 0)) {
      ## browser()
      stop("floating.pie: x values must be non-negative")
    }
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    if (is.null(col)) col <- rainbow(nx)
    else if (length(col) < nx) col <- rep(col, nx)
    bc <- 2 * pi * (x[1:nx] + dx/2) + startpos
    for (i in 1:nx) {
        n <- max(2, floor(edges * dx[i]))
        t2p <- 2 * pi * seq(x[i], x[i + 1], length = n) + startpos
        xc <- c(cos(t2p) * radius + xpos, xpos)
        yc <- c(sin(t2p) * radius*asp + ypos, ypos)
        polygon(xc, yc, col = col[i], ...)
        ## t2p <- 2 * pi * mean(x[i + 0:1]) + startpos
        ## xc <- cos(t2p) * radius
        ## yc <- sin(t2p) * radius*asp
        ## lines(c(1, 1.05) * xc, c(1, 1.05) * yc)
    }
    ## return(bc)
}

BOTHlabels <- function(text, sel, XX, YY, adj, frame, pch, thermo,
                       pie, piecol, col, bg, ...)
{
    if (missing(text)) text <- NULL
    if (length(adj) == 1) adj <- c(adj, 0.5)
    if (is.null(text) && is.null(pch) && is.null(thermo) && is.null(pie))
      text <- as.character(sel)
    frame <- match.arg(frame, c("rect", "circle", "none"))
    args <- list(...)
    CEX <- if ("cex" %in% names(args)) args$cex else par("cex")
    if (frame != "none" && !is.null(text)) {
        if (frame == "rect") {
            width <- strwidth(text, units = "inches", cex = CEX)
            height <- strheight(text, units = "inches", cex = CEX)
            if ("srt" %in% names(args)) {
                args$srt <- args$srt %% 360 # just in case srt >= 360
                if (args$srt == 90 || args$srt == 270) {
                    tmp <- width
                    width <- height
                    height <- tmp
                } else if (args$srt != 0)
                  warning("only right angle rotation of frame is supported;\n         try  `frame = \"n\"' instead.\n")
            }
            width <- xinch(width)
            height <- yinch(height)
            xl <- XX - width*adj[1] - xinch(0.03)
            xr <- xl + width + xinch(0.03)
            yb <- YY - height*adj[2] - yinch(0.02)
            yt <- yb + height + yinch(0.05)
            rect(xl, yb, xr, yt, col = bg)
        }
        if (frame == "circle") {
            radii <- 0.8*apply(cbind(strheight(text, units = "inches", cex = CEX),
                                     strwidth(text, units = "inches", cex = CEX)), 1, max)
            symbols(XX, YY, circles = radii, inches = max(radii), add = TRUE, bg = bg)
        }
    }
    if (!is.null(thermo)) {
        parusr <- par("usr")
        width <- CEX * (parusr[2] - parusr[1]) / 40
        height <- CEX * (parusr[4] - parusr[3]) / 15
        if (is.vector(thermo)) thermo <- cbind(thermo, 1 - thermo)
        thermo <- height * thermo
        xl <- XX - width/2
        xr <- xl + width
        yb <- YY - height/2
        if (is.null(piecol)) piecol <- rainbow(ncol(thermo))
        ## draw the first rectangle:
        rect(xl, yb, xr, yb + thermo[, 1], border = NA, col = piecol[1])
        for (i in 2:ncol(thermo))
          rect(xl, yb + rowSums(thermo[, 1:(i - 1), drop = FALSE]),
               xr, yb + rowSums(thermo[, 1:i]),
               border = NA, col = piecol[i])
        rect(xl, yb, xr, yb + height, border = "black")
        segments(xl, YY, xl - width/5, YY)
        segments(xr, YY, xr + width/5, YY)
    }
    ## from BB:
    if (!is.null(pie)) {
        if (is.vector(pie)) pie <- cbind(pie, 1 - pie)
        xrad <- CEX * diff(par("usr")[1:2]) / 50
        for (i in 1:length(sel))
          floating.pie.asp(XX[i], YY[i], pie[i, ],
                           radius = xrad, col = piecol)
    }
    if (!is.null(text)) text(XX, YY, text, adj = adj, col = col, ...)
    if (!is.null(pch)) points(XX + adj[1] - 0.5, YY + adj[2] - 0.5,
                              pch = pch, col = col, bg = bg, ...)
}

nodelabels <- function(text, node, adj = c(0.5, 0.5), frame = "rect",
                       pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
                       col = "black", bg = "lightblue", ...)
{
    if (missing(node))
      node <- (.last_plot.phylo$Ntip + 1):length(.last_plot.phylo$xx)
    XX <- .last_plot.phylo$xx[node]
    YY <- .last_plot.phylo$yy[node]
    BOTHlabels(text, node, XX, YY, adj, frame, pch, thermo,
               pie, piecol, col, bg, ...)
}

tiplabels <- function(text, tip, adj = c(0.5, 0.5), frame = "rect",
                      pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
                      col = "black", bg = "yellow", ...)
{
    if (missing(tip)) tip <- 1:.last_plot.phylo$Ntip
    XX <- .last_plot.phylo$xx[tip]
    YY <- .last_plot.phylo$yy[tip]
    BOTHlabels(text, tip, XX, YY, adj, frame, pch, thermo,
               pie, piecol, col, bg, ...)
}

edgelabels <- function(text, edge, adj = c(0.5, 0.5), frame = "rect",
                      pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
                      col = "black", bg = "lightgreen", ...)
{
    if (missing(edge)) {
        sel <- 1:dim(.last_plot.phylo$edge)[1]
        subedge <- .last_plot.phylo$edge
    } else {
        sel <- edge
        subedge <- .last_plot.phylo$edge[sel, , drop = FALSE]
    }
    if (.last_plot.phylo$type == "phylogram") {
        if(.last_plot.phylo$direction %in% c("rightwards", "leftwards")) {
            XX <- (.last_plot.phylo$xx[subedge[, 1]] +
                   .last_plot.phylo$xx[subedge[, 2]]) / 2
            YY <- .last_plot.phylo$yy[subedge[, 2]]
        } else {
            XX <- .last_plot.phylo$xx[subedge[, 2]]
            YY <- (.last_plot.phylo$yy[subedge[, 1]] +
                   .last_plot.phylo$yy[subedge[, 2]]) / 2
        }
    } else {
        XX <- (.last_plot.phylo$xx[subedge[, 1]] +
               .last_plot.phylo$xx[subedge[, 2]]) / 2
        YY <- (.last_plot.phylo$yy[subedge[, 1]] +
               .last_plot.phylo$yy[subedge[, 2]]) / 2
    }
    BOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo,
               pie, piecol, col, bg, ...)
}

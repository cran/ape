## plot.phylo.R (2011-12-03)

##   Plot Phylogenies

## Copyright 2002-2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plot.phylo <-
    function(x, type = "phylogram", use.edge.length = TRUE,
             node.pos = NULL, show.tip.label = TRUE,
             show.node.label = FALSE, edge.color = "black",
             edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"),
             adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE,
             label.offset = 0, underscore = FALSE, x.lim = NULL,
             y.lim = NULL, direction = "rightwards", lab4ut = "horizontal",
             tip.color = "black", plot = TRUE, rotate.tree = 0, ...)
{
    Ntip <- length(x$tip.label)
    if (Ntip == 1) {
        warning("found only one tip in the tree")
        return(NULL)
    }
    if (any(tabulate(x$edge[, 1]) == 1))
      stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")

    .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy)
        .C("node_height", as.integer(Ntip), as.integer(Nnode),
           as.integer(edge[, 1]), as.integer(edge[, 2]),
           as.integer(Nedge), as.double(yy),
           DUP = FALSE, PACKAGE = "ape")[[6]]

    .nodeDepth <- function(Ntip, Nnode, edge, Nedge)
        .C("node_depth", as.integer(Ntip), as.integer(Nnode),
           as.integer(edge[, 1]), as.integer(edge[, 2]),
           as.integer(Nedge), double(Ntip + Nnode),
           DUP = FALSE, PACKAGE = "ape")[[6]]

    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, edge.length)
        .C("node_depth_edgelength", as.integer(Ntip),
           as.integer(Nnode), as.integer(edge[, 1]),
           as.integer(edge[, 2]), as.integer(Nedge),
           as.double(edge.length), double(Ntip + Nnode),
           DUP = FALSE, PACKAGE = "ape")[[7]]

    Nedge <- dim(x$edge)[1]
    Nnode <- x$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan",
                              "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards",
                                        "upwards", "downwards"))
    if (is.null(x$edge.length)) use.edge.length <- FALSE

    ## the order of the last two conditions is important:
    if (type %in% c("unrooted", "radial") || !use.edge.length ||
        is.null(x$root.edge) || !x$root.edge) root.edge <- FALSE
    if (type == "fan" && root.edge) {
        warning("drawing root edge with type = 'fan' is not yet supported")
        root.edge <- FALSE
    }

    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- x$edge # to save
    if (phyloORclado) {
        ## we first compute the y-coordinates of the tips.
        phyOrder <- attr(x, "order")
        ## make sure the tree is in cladewise order:
        if (is.null(phyOrder) || phyOrder != "cladewise") {
            x <- reorder(x) # fix from Klaus Schliep (2007-06-16)
            if (!identical(x$edge, xe)) {
                ## modified from Li-San Wang's fix (2007-01-23):
                ereorder <- match(x$edge[, 2], xe[, 2])
                if (length(edge.color) > 1) {
                    edge.color <- rep(edge.color, length.out = Nedge)
                    edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                    edge.width <- rep(edge.width, length.out = Nedge)
                    edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                    edge.lty <- rep(edge.lty, length.out = Nedge)
                    edge.lty <- edge.lty[ereorder]
                }
            }
        }
### By contrats to ape (< 2.4), the arguments edge.color, etc., are
### not elongated before being passed to segments(), except if needed
### to be reordered
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    ## 'z' is the tree in pruningwise order used in calls to .C
    z <- reorder(x, order = "pruningwise")

    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) node.pos <- 2
        }
        if (node.pos == 1)
            yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
        else {
          ## node_height_clado requires the number of descendants
          ## for each node, so we compute `xx' at the same time
          ans <- .C("node_height_clado", as.integer(Ntip),
                    as.integer(Nnode), as.integer(z$edge[, 1]),
                    as.integer(z$edge[, 2]), as.integer(Nedge),
                    double(Ntip + Nnode), as.double(yy),
                    DUP = FALSE, PACKAGE = "ape")
          xx <- ans[[6]] - 1
          yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 1
            xx <- max(xx) - xx
        } else  {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, z$edge.length)
        }
    } else {
    rotate.tree <- 2 * pi * rotate.tree/360
    switch(type, "fan" = {
        ## if the tips are not in the same order in tip.label
        ## and in edge[, 2], we must reorder the angles: we
        ## use `xx' to store temporarily the angles
        TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
        xx <- seq(0, 2*pi*(1 - 1/Ntip), 2*pi/Ntip)
        theta <- double(Ntip)
        theta[TIPS] <- xx
        theta <- c(theta, numeric(Nnode))
        theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, theta)
        if (use.edge.length) {
            r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, z$edge.length)
        } else {
            r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            r <- 1/r
        }
        theta <- theta + rotate.tree
        xx <- r*cos(theta)
        yy <- r*sin(theta)
    }, "unrooted" = {
        nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
        XY <- if (use.edge.length)
            unrooted.xy(Ntip, Nnode, z$edge, z$edge.length, nb.sp, rotate.tree)
        else
            unrooted.xy(Ntip, Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
        ## rescale so that we have only positive values
        xx <- XY$M[, 1] - min(XY$M[, 1])
        yy <- XY$M[, 2] - min(XY$M[, 2])
    }, "radial" = {
        X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
        X[X == 1] <- 0
        ## radius:
        X <- 1 - X/Ntip
        ## angle (1st compute the angles for the tips):
        yy <- c((1:Ntip)*2*pi/Ntip, rep(0, Nnode))
        Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
        xx <- X * cos(Y + rotate.tree)
        yy <- X * sin(Y + rotate.tree)
    })}
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == "rightwards") xx <- xx + x$root.edge
            if (direction == "upwards") yy <- yy + x$root.edge
        }
    }
    if (no.margin) par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                pin1 <- par("pin")[1] # width of the device in inches
                strWi <- strwidth(x$tip.label, "inches") # id. for the tip labels
                ## 1.04 comes from that we are using a regular axis system
                ## with 4% on both sides of the range of x:
                xx.tips <- xx[1:Ntip] * 1.04
                ## 'alp' is the conversion coefficient from
                ## user coordinates to inches:
                alp <- try(uniroot(function(a) max(a*xx.tips + strWi) - pin1,
                                   c(0, 1e6))$root, silent = TRUE)
                ## if the above fails, give 1/3 of the device for the tip labels:
                if (is.character(alp)) tmp <- max(xx.tips)*1.5 else {
                    tmp <- if (show.tip.label) max(xx.tips + strWi/alp) else max(xx.tips)
                }
                x.lim[2] <- tmp
            } else x.lim <- c(1, Ntip)
        } else switch(type, "fan" = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
                x.lim <- c(min(xx) - offset, max(xx) + offset)
            } else x.lim <- c(min(xx), max(xx))
        }, "unrooted" = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
                x.lim <- c(0 - offset, max(xx) + offset)
            } else x.lim <- c(0, max(xx))
        }, "radial" = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset)
            } else x.lim <- c(-1, 1)
        })
    } else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label)
          x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
        if (type == "radial")
          x.lim[1] <-
            if (show.tip.label) -1 - max(nchar(x$tip.label) * 0.03 * cex)
            else -1
    }
    ## mirror the xx:
    if (phyloORclado && direction == "leftwards") xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) y.lim <- c(1, Ntip) else {
                y.lim <- c(0, NA)
                pin2 <- par("pin")[2] # height of the device in inches
                strWi <- strwidth(x$tip.label, "inches")
                ## 1.04 comes from that we are using a regular axis system
                ## with 4% on both sides of the range of x:
                yy.tips <- yy[1:Ntip] * 1.04
                ## 'alp' is the conversion coefficient from
                ## user coordinates to inches:
                alp <- try(uniroot(function(a) max(a*yy.tips + strWi) - pin2,
                                   c(0, 1e6))$root, silent = TRUE)
                ## if the above fails, give 1/3 of the device for the tip labels:
                if (is.character(alp)) tmp <- max(yy.tips)*1.5 else {
                    tmp <- if (show.tip.label) max(yy.tips + strWi/alp) else max(yy.tips)
                }
                y.lim[2] <- tmp
            }
        } else switch(type, "fan" = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
                y.lim <- c(min(yy) - offset, max(yy) + offset)
            } else y.lim <- c(min(yy), max(yy))
        }, "unrooted" = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
                y.lim <- c(0 - offset, max(yy) + offset)
            } else y.lim <- c(0, max(yy))
        }, "radial" = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            } else y.lim <- c(-1, 1)
        })
    } else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label)
          y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
        if (type == "radial")
          y.lim[1] <- if (show.tip.label) -1 - max(nchar(x$tip.label) * 0.018 * max(yy) * cex) else -1
    }
    ## mirror the yy:
    if (phyloORclado && direction == "downwards") yy <- y.lim[2] - yy
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") y.lim[2] <- y.lim[2] + x$root.edge
    }
    asp <- if (type %in% c("fan", "radial", "unrooted")) 1 else NA # fixes by Klaus Schliep (2008-03-28 and 2010-08-12)
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, ann = FALSE, axes = FALSE, asp = asp, ...)

if (plot) {
    if (is.null(adj))
        adj <- if (phyloORclado && direction == "leftwards") 1 else 0
    if (phyloORclado && show.tip.label) {
        MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
        loy <- 0
        if (direction == "rightwards") {
            lox <- label.offset + MAXSTRING * 1.05 * adj
        }
        if (direction == "leftwards") {
            lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
            ##xx <- xx + MAXSTRING
        }
        if (!horizontal) {
            psr <- par("usr")
            MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - psr[1])
            loy <- label.offset + MAXSTRING * 1.05 * adj
            lox <- 0
            srt <- 90 + srt
            if (direction == "downwards") {
                loy <- -loy
                ##yy <- yy + MAXSTRING
                srt <- 180 + srt
            }
        }
    }
    if (type == "phylogram") {
        phylogram.plot(x$edge, Ntip, Nnode, xx, yy,
                       horizontal, edge.color, edge.width, edge.lty)
    } else {
        if (type == "fan") {
            ereorder <- match(z$edge[, 2], x$edge[, 2])
            if (length(edge.color) > 1) {
                edge.color <- rep(edge.color, length.out = Nedge)
                edge.color <- edge.color[ereorder]
            }
            if (length(edge.width) > 1) {
                edge.width <- rep(edge.width, length.out = Nedge)
                edge.width <- edge.width[ereorder]
            }
            if (length(edge.lty) > 1) {
                edge.lty <- rep(edge.lty, length.out = Nedge)
                edge.lty <- edge.lty[ereorder]
            }
            circular.plot(z$edge, Ntip, Nnode, xx, yy, theta,
                          r, edge.color, edge.width, edge.lty)
        } else
        cladogram.plot(x$edge, xx, yy, edge.color, edge.width, edge.lty)
    }
    if (root.edge)
      switch(direction,
             "rightwards" = segments(0, yy[ROOT], x$root.edge, yy[ROOT]),
             "leftwards" = segments(xx[ROOT], yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT]),
             "upwards" = segments(xx[ROOT], 0, xx[ROOT], x$root.edge),
             "downwards" = segments(xx[ROOT], yy[ROOT], xx[ROOT], yy[ROOT] + x$root.edge))
    if (show.tip.label) {
        if (is.expression(x$tip.label)) underscore <- TRUE
        if (!underscore) x$tip.label <- gsub("_", " ", x$tip.label)

        if (phyloORclado)
            text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, adj = adj,
                 font = font, srt = srt, cex = cex, col = tip.color)

        if (type == "unrooted") {
            if (lab4ut == "horizontal") {
                y.adj <- x.adj <- numeric(Ntip)
                sel <- abs(XY$axe) > 0.75 * pi
                x.adj[sel] <- -strwidth(x$tip.label)[sel] * 1.05
                sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * pi
                x.adj[sel] <- -strwidth(x$tip.label)[sel] * (2 * abs(XY$axe)[sel] / pi - 0.5)
                sel <- XY$axe > pi / 4 & XY$axe < 0.75 * pi
                y.adj[sel] <- strheight(x$tip.label)[sel] / 2
                sel <- XY$axe < -pi / 4 & XY$axe > -0.75 * pi
                y.adj[sel] <- -strheight(x$tip.label)[sel] * 0.75
                text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + y.adj * cex,
                     x$tip.label, adj = c(adj, 0), font = font,
                     srt = srt, cex = cex, col = tip.color)
            } else { # if lab4ut == "axial"
                adj <- abs(XY$axe) > pi/2
                srt <- 180 * XY$axe / pi
                srt[adj] <- srt[adj] - 180
                adj <- as.numeric(adj)
                xx.tips <- xx[1:Ntip]
                yy.tips <- yy[1:Ntip]
                if (label.offset) {
                    xx.tips <- xx.tips + label.offset * cos(XY$axe)
                    yy.tips <- yy.tips + label.offset * sin(XY$axe)
                }
                ## `srt' takes only a single value, so can't vectorize this:
                ## (and need to 'elongate' these vectors:)
                font <- rep(font, length.out = Ntip)
                tip.color <- rep(tip.color, length.out = Ntip)
                cex <- rep(cex, length.out = Ntip)
                for (i in 1:Ntip)
                    text(xx.tips[i], yy.tips[i], cex = cex[i],
                         x$tip.label[i], adj = adj[i], font = font[i],
                         srt = srt[i], col = tip.color[i])
            }
        }
        if (type %in% c("fan", "radial")) {
            xx.tips <- xx[1:Ntip]
            yy.tips <- yy[1:Ntip]
            ## using atan2 considerably facilitates things compared to acos...
            angle <- atan2(yy.tips, xx.tips) # in radians
            if (label.offset) {
                xx.tips <- xx.tips + label.offset * cos(angle)
                yy.tips <- yy.tips + label.offset * sin(angle)
            }
            s <- xx.tips < 0
            angle <- angle * 180/pi # switch to degrees
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
            ## `srt' takes only a single value, so can't vectorize this:
            ## (and need to 'elongate' these vectors:)
            font <- rep(font, length.out = Ntip)
            tip.color <- rep(tip.color, length.out = Ntip)
            cex <- rep(cex, length.out = Ntip)
            for (i in 1:Ntip)
                text(xx.tips[i], yy.tips[i], x$tip.label[i], font = font[i],
                     cex = cex[i], srt = angle[i], adj = adj[i],
                     col = tip.color[i])
        }
    }
    if (show.node.label)
        text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)],
             x$node.label, adj = adj, font = font, srt = srt, cex = cex)
}
    L <- list(type = type, use.edge.length = use.edge.length,
              node.pos = node.pos, show.tip.label = show.tip.label,
              show.node.label = show.node.label, font = font,
              cex = cex, adj = adj, srt = srt, no.margin = no.margin,
              label.offset = label.offset, x.lim = x.lim, y.lim = y.lim,
              direction = direction, tip.color = tip.color,
              Ntip = Ntip, Nnode = Nnode)
    assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)),
           envir = .PlotPhyloEnv)
    invisible(L)
}

phylogram.plot <- function(edge, Ntip, Nnode, xx, yy, horizontal,
                           edge.color, edge.width, edge.lty)
{
    nodes <- (Ntip + 1):(Ntip + Nnode)
    if (!horizontal) {
        tmp <- yy
        yy <- xx
        xx <- tmp
    }
    ## un trait vertical à chaque noeud...
    x0v <- xx[nodes]
    y0v <- y1v <- numeric(Nnode)
    ## store the index of each node in the 1st column of edge:
    NodeInEdge1 <- vector("list", Nnode)
    for (i in nodes) {
        ii <- i - Ntip
        j <- NodeInEdge1[[ii]] <- which(edge[, 1] == i)
        tmp <- range(yy[edge[j, 2]])
        y0v[ii] <- tmp[1]
        y1v[ii] <- tmp[2]
    }
    ## ... et un trait horizontal partant de chaque tip et chaque noeud
    ##  vers la racine
    x0h <- xx[edge[, 1]]
    x1h <- xx[edge[, 2]]
    y0h <- yy[edge[, 2]]

    nc <- length(edge.color)
    nw <- length(edge.width)
    nl <- length(edge.lty)

    if (nc + nw + nl == 3) {
        color.v <- edge.color
        width.v <- edge.width
        lty.v <- edge.lty
    } else {
        Nedge <- dim(edge)[1]
        edge.color <- rep(edge.color, length.out = Nedge)
        edge.width <- rep(edge.width, length.out = Nedge)
        edge.lty <- rep(edge.lty, length.out = Nedge)
        DF <- data.frame(edge.color, edge.width, edge.lty, stringsAsFactors = FALSE)
        color.v <- rep("black", Nnode)
        width.v <- rep(1, Nnode)
        lty.v <- rep(1, Nnode)
        for (i in 1:Nnode) {
            br <- NodeInEdge1[[i]]
            if (length(br) > 2) {
                x <- unique(DF[br, 1])
                if (length(x) == 1) color.v[i] <- x
                x <- unique(DF[br, 2])
                if (length(x) == 1) width.v[i] <- x
                x <- unique(DF[br, 3])
                if (length(x) == 1) lty.v[i] <- x
            } else {
                A <- br[1]
                B <- br[2]
                if (any(DF[A, ] != DF[B, ])) {
                    color.v[i] <- edge.color[B]
                    width.v[i] <- edge.width[B]
                    lty.v[i] <- edge.lty[B]
                    ## add a new line:
                    y0v <- c(y0v, y0v[i])
                    y1v <- c(y1v, yy[i + Ntip])
                    x0v <- c(x0v, x0v[i])
                    color.v <- c(color.v, edge.color[A])
                    width.v <- c(width.v, edge.width[A])
                    lty.v <- c(lty.v, edge.lty[A])
                    ## shorten the line:
                    y0v[i] <- yy[i + Ntip]
                } else {
                    color.v[i] <- edge.color[A]
                    width.v[i] <- edge.width[A]
                    lty.v[i] <- edge.lty[A]
                }
            }
        }
    }

    if (horizontal) {
        segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width, lty = edge.lty) # draws horizontal lines
        segments(x0v, y0v, x0v, y1v, col = color.v, lwd = width.v, lty = lty.v) # draws vertical lines
    } else {
        segments(y0h, x0h, y0h, x1h, col = edge.color, lwd = edge.width, lty = edge.lty) # draws vertical lines
        segments(y0v, x0v, y1v, x0v, col = color.v, lwd = width.v, lty = lty.v) # draws horizontal lines
    }
}

cladogram.plot <- function(edge, xx, yy, edge.color, edge.width, edge.lty)
    segments(xx[edge[, 1]], yy[edge[, 1]], xx[edge[, 2]], yy[edge[, 2]],
             col = edge.color, lwd = edge.width, lty = edge.lty)

circular.plot <- function(edge, Ntip, Nnode, xx, yy, theta,
                          r, edge.color, edge.width, edge.lty)
### 'edge' must be in pruningwise order
{
    r0 <- r[edge[, 1]]
    r1 <- r[edge[, 2]]
    theta0 <- theta[edge[, 2]]
    costheta0 <- cos(theta0)
    sintheta0 <- sin(theta0)

    x0 <- r0 * costheta0
    y0 <- r0 * sintheta0
    x1 <- r1 * costheta0
    y1 <- r1 * sintheta0

    segments(x0, y0, x1, y1, col = edge.color, lwd = edge.width, lty = edge.lty)

    tmp <- which(diff(edge[, 1]) != 0)
    start <- c(1, tmp + 1)
    Nedge <- dim(edge)[1]
    end <- c(tmp, Nedge)

    ## function dispatching the features to the arcs
    foo <- function(edge.feat, default) {
        if (length(edge.feat) == 1) return(rep(edge.feat, Nnode))
        else {
            edge.feat <- rep(edge.feat, length.out = Nedge)
            feat.arc <- rep(default, Nnode)
            for (k in 1:Nnode) {
                tmp <- edge.feat[start[k]]
                if (tmp == edge.feat[end[k]]) feat.arc[k] <- tmp
            }
        }
        feat.arc
    }
    co <- foo(edge.color, "black")
    lw <- foo(edge.width, 1)
    ly <- foo(edge.lty, 1)

    for (k in 1:Nnode) {
        i <- start[k]
        j <- end[k]
        X <- rep(r[edge[i, 1]], 100)
        Y <- seq(theta[edge[i, 2]], theta[edge[j, 2]], length.out = 100)
        lines(X*cos(Y), X*sin(Y), col = co[k], lwd = lw[k], lty = ly[k])
    }
}

unrooted.xy <- function(Ntip, Nnode, edge, edge.length, nb.sp, rotate.tree)
{
    foo <- function(node, ANGLE, AXIS) {
        ind <- which(edge[, 1] == node)
        sons <- edge[ind, 2]
        start <- AXIS - ANGLE/2
        for (i in 1:length(sons)) {
            h <- edge.length[ind[i]]
            angle[sons[i]] <<- alpha <- ANGLE*nb.sp[sons[i]]/nb.sp[node]
            axis[sons[i]] <<- beta <- start + alpha/2
            start <- start + alpha
            xx[sons[i]] <<- h*cos(beta) + xx[node]
            yy[sons[i]] <<- h*sin(beta) + yy[node]
        }
        for (i in sons)
            if (i > Ntip) foo(i, angle[i], axis[i])
    }
    Nedge <- dim(edge)[1]
    yy <- xx <- numeric(Ntip + Nnode)
    ## `angle': the angle allocated to each node wrt their nb of tips
    ## `axis': the axis of each branch
    axis <- angle <- numeric(Ntip + Nnode)
    ## start with the root...
    foo(Ntip + 1L, 2*pi, 0 + rotate.tree)

    M <- cbind(xx, yy)
    axe <- axis[1:Ntip] # the axis of the terminal branches (for export)
    axeGTpi <- axe > pi
    ## insures that returned angles are in [-PI, +PI]:
    axe[axeGTpi] <- axe[axeGTpi] - 2*pi
    list(M = M, axe = axe)
}

node.depth <- function(phy)
{
    n <- length(phy$tip.label)
    m <- phy$Nnode
    N <- dim(phy$edge)[1]
    phy <- reorder(phy, order = "pruningwise")
    .C("node_depth", as.integer(n), as.integer(m),
       as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
       as.integer(N), double(n + m), DUP = FALSE, PACKAGE = "ape")[[6]]
}

plot.multiPhylo <- function(x, layout = 1, ...)
{
    if (layout > 1)
      layout(matrix(1:layout, ceiling(sqrt(layout)), byrow = TRUE))
    else layout(matrix(1))
    if (!par("ask")) {
        par(ask = TRUE)
        on.exit(par(ask = FALSE))
    }
    for (i in 1:length(x)) plot(x[[i]], ...)
}

trex <- function(phy, title = TRUE, subbg = "lightyellow3",
                 return.tree = FALSE, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    devmain <- dev.cur() # where the main tree is plotted

    restore <- function() {
        dev.set(devmain)
        assign("last_plot.phylo", lastPP, envir = .PlotPhyloEnv)
    }

    on.exit(restore())
    NEW <- TRUE
    cat("Click close to a node. Right-click to exit.\n")
    repeat {
        x <- identify.phylo(phy, quiet = TRUE)
        if (is.null(x)) return(invisible(NULL)) else {
            x <- x$nodes
            if (is.null(x)) cat("Try again!\n") else {
                if (NEW) {
                    dev.new()
                    par(bg = subbg)
                    devsub <- dev.cur()
                    NEW <- FALSE
                } else dev.set(devsub)

                tr <- extract.clade(phy, x)
                plot(tr, ...)
                if (is.character(title)) title(title)
                else if (title) {
                     tl <-
                         if (is.null(phy$node.label))
                         paste("From node #", x, sep = "")
                         else paste("From", phy$node.label[x - Ntip(phy)])
                     title(tl)
                }
                if (return.tree) return(tr)
                restore()
            }
        }
    }
}

kronoviz <- function(x, layout = length(x), horiz = TRUE, ...)
{
    par(mar = rep(0.5, 4), oma = rep(2, 4))
    rts <- sapply(x, function(x) branching.times(x)[1])
    maxrts <- max(rts)
    lim <- cbind(rts - maxrts, rts)
    Ntree <- length(x)
    Ntips <- sapply(x, Ntip)
    if (horiz) {
        nrow <- layout
        w <- 1
        h <- Ntips
    } else {
        nrow <- 1
        w <- Ntips
        h <- 1
    }
    layout(matrix(1:layout, nrow), widths = w, heights = h)
    if (layout > Ntree && !par("ask")) {
        par(ask = TRUE)
        on.exit(par(ask = FALSE))
    }
    if (horiz) {
        for (i in 1:Ntree)
            plot(x[[i]], x.lim = lim[i, ], ...)
    } else {
        for (i in 1:Ntree)
            plot(x[[i]], y.lim = lim[i, ], direction = "u", ...)
    }
    axisPhylo(if (horiz) 1 else 4) # better if the deepest tree is last ;)
}

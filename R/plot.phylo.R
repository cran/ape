### plot.phylo.R (2006-11-13)
###
###          Plot Phylogenies
###
### Copyright 2002-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

plot.phylo <- function(x, type = "phylogram", use.edge.length = TRUE,
                       node.pos = NULL, show.tip.label = TRUE,
                       show.node.label = FALSE, edge.color = "black",
                       edge.width = 1, font = 3, cex = par("cex"),
                       adj = NULL, srt = 0, no.margin = FALSE,
                       root.edge = FALSE, label.offset = 0, underscore = FALSE,
                       x.lim = NULL, y.lim = NULL, direction = "rightwards",
                       lab4ut = "horizontal", tip.color = "black", ...)
{
    nb.tip <- length(x$tip.label)
    if (nb.tip == 1) stop("found only one tip in the tree!")
    nb.edge <- dim(x$edge)[1]
    if (any(tabulate(x$edge[, 1]) == 1))
      stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles().")
    nb.node <- x$Nnode
    root <- nb.tip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards",
                                        "upwards", "downwards"))
    if (is.null(x$edge.length)) use.edge.length <- FALSE

    ## <FIXME> Maybe it's possible to simplify this:
    if (!show.tip.label) x$tip.label <- rep("", nb.tip)
    ## </FIXME>

    if (type == "unrooted" || !use.edge.length) root.edge <- FALSE
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    if (phyloORclado) {
        ## we first compute the y-coordinates of the tips, this
        ## doesn't require to have the tips ordered (2006-10-13)
        yy <- numeric(nb.tip + nb.node)
        TIPS <- x$edge[x$edge[, 2] <= nb.tip, 2]
        yy[TIPS] <- 1:nb.tip
    }
    x <- reorder(x, order = "pruningwise")
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) node.pos <- 2
        }
        if (node.pos == 1)
          yy <- .C("node_height", as.integer(nb.tip), as.integer(nb.node),
                   as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
                   as.integer(nb.edge), as.double(yy),
                   DUP = FALSE, PACKAGE = "ape")[[6]]
        else {
          ## node_height_clado requires the number of descendants
          ## for each node, so we compute `xx' at the same time
          ans <- .C("node_height_clado", as.integer(nb.tip),
                    as.integer(nb.node), as.integer(x$edge[, 1]),
                    as.integer(x$edge[, 2]), as.integer(nb.edge),
                    double(nb.tip + nb.node), as.double(yy),
                    DUP = FALSE, PACKAGE = "ape")
          xx <- ans[[6]] - 1
          yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if(node.pos != 2)
              xx <- .C("node_depth", as.integer(nb.tip), as.integer(nb.node),
                       as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
                       as.integer(nb.edge), double(nb.tip + nb.node),
                       DUP = FALSE, PACKAGE = "ape")[[6]] - 1
            xx <- max(xx) - xx
          } else  {
            xx <- .C("node_depth_edgelength", as.integer(nb.tip),
                     as.integer(nb.node), as.integer(x$edge[, 1]),
                     as.integer(x$edge[, 2]), as.integer(nb.edge),
                     as.double(x$edge.length), double(nb.tip + nb.node),
                     DUP = FALSE, PACKAGE = "ape")[[7]]
        }
    }
    if (type == "unrooted") {
        XY <- if (use.edge.length) unrooted.xy(nb.tip, nb.node, x$edge, x$edge.length) else unrooted.xy(nb.tip, nb.node, x$edge, rep(1, nb.edge))
        ## rescale so that we have only positive values
        xx <- XY$M[, 1] - min(XY$M[, 1])
        yy <- XY$M[, 2] - min(XY$M[, 2])
    }
    if (type == "radial") {
        X <- .C("node_depth", as.integer(nb.tip), as.integer(nb.node),
                as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
                as.integer(nb.edge), double(nb.tip + nb.node),
                DUP = FALSE, PACKAGE = "ape")[[6]]
        X[X == 1] <- 0
        ## radius:
        X <- 1 - X/nb.tip
        ## angle (1st compute the angles for the tips):
        yy <- c((1:nb.tip)*2*pi/nb.tip, rep(0, nb.node))
        Y <- .C("node_height", as.integer(nb.tip), as.integer(nb.node),
                as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
                as.integer(nb.edge), as.double(yy),
                DUP = FALSE, PACKAGE = "ape")[[6]]
        xx <- X * cos(Y)
        yy <- X * sin(Y)
    }
    if (phyloORclado && direction != "rightwards") {
        if (direction == "leftwards") {
            xx <- -xx
            xx <- xx - min(xx)
        }
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
            if (direction == "downwards") {
                yy <- -yy
                yy <- yy - min(yy)
            }
        }
    }
    if (phyloORclado && root.edge) {
        if (direction == "rightwards") xx <- xx + x$root.edge
        if (direction == "upwards") yy <- yy + x$root.edge
    }
    if (no.margin) par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                tmp <- nchar(x$tip.label) * 0.018 * max(xx) * cex
                x.lim[2] <-
                  if (direction == "leftwards") max(xx[root] + tmp) else max(xx[1:nb.tip] + tmp)
            } else x.lim <- c(1, nb.tip)
        }
        if (type == "unrooted") {
            offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
            x.lim <- c(0 - offset, max(xx) + offset)
        }
        if (type == "radial") {
            offset <- max(nchar(x$tip.label) * 0.03 * cex)
            x.lim <- c(-1 - offset, 1 + offset)
        }
    } else  {
        if (length(x.lim) == 1) {
            if (phyloORclado)
              x.lim <- if (horizontal) c(0, x.lim) else c(1, x.lim)
            if (type == "unrooted") {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
                x.lim <- c(0 - offset, x.lim)
            }
            if (type == "radial") {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, x.lim)
            }
        }
    }
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (!horizontal) {
                y.lim <- c(0, NA)
                y.lim[2] <- if (direction == "downwards")
                  max(yy[root] + nchar(x$tip.label) * 0.018 * max(yy) * cex)
                else
                  max(yy[1:nb.tip] + nchar(x$tip.label) * 0.018 * max(yy) * cex)
            } else y.lim <- c(1, nb.tip)
        }
        if (type == "unrooted") {
            offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
            y.lim <- c(0 - offset, max(yy) + offset)
        }
        if (type == "radial") {
            offset <- max(nchar(x$tip.label) * 0.03 * cex)
            y.lim <- c(-1 - offset, 1 + offset)
        }
    } else {
        if (length(y.lim) == 1) {
            if (phyloORclado)
              y.lim <- if (horizontal)
                c(1, y.lim) else c(0, y.lim)
            if (type == "unrooted") {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
                y.lim <- c(0 - offset, y.lim)
            }
            if (type == "radial") {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * cex)
                y.lim <- c(-1 - offset, y.lim)
            }
        }
    }
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") y.lim[2] <- y.lim[2] + x$root.edge
    }

    edge.color <- rep(edge.color, length.out = nb.edge)
    edge.width <- rep(edge.width, length.out = nb.edge)
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "",
         ylab = "", xaxt = "n", yaxt = "n", bty = "n", ...)
    if (is.null(adj))
      adj <- if (phyloORclado && direction == "leftwards") 1 else 0
    if (phyloORclado) {
        MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
        if (direction == "rightwards") {
            lox <- label.offset + MAXSTRING * 1.05 * adj
            loy <- 0
        }
        if (direction == "leftwards") {
            lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
            loy <- 0
            xx <- xx + MAXSTRING
        }
        if (!horizontal) {
            psr <- par("usr")
            MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3]) / (psr[2] - psr[1])
            loy <- label.offset + MAXSTRING * 1.05 * adj
            lox <- 0
            srt <- 90 + srt
            if (direction == "downwards") {
                loy <- -loy
                yy <- yy + MAXSTRING
                srt <- 180 + srt
            }
        }
    }
    if (type == "phylogram")
      phylogram.plot(x$edge, nb.tip, nb.node, xx, yy,
                     horizontal, edge.color, edge.width)
    else cladogram.plot(x$edge, xx, yy, edge.color, edge.width)
    if (root.edge)
      switch(direction,
             "rightwards" = segments(0, yy[root], x$root.edge, yy[root]),
             "leftwards" = segments(xx[root], yy[root], xx[root] + x$root.edge, yy[root]),
             "upwards" = segments(xx[root], 0, xx[root], x$root.edge),
             "downwards" = segments(xx[root], yy[root], xx[root], yy[root] + x$root.edge))
    if (!underscore) x$tip.label <- gsub("_", " ", x$tip.label)
    if (phyloORclado) {
        text(xx[1:nb.tip] + lox, yy[1:nb.tip] + loy, x$tip.label, adj = adj,
             font = font, srt = srt, cex = cex, col = tip.color)
    }
    if (type == "unrooted") {
        if (lab4ut == "horizontal") {
            y.adj <- x.adj <- numeric(nb.tip)
            sel <- abs(XY$axe) > 0.75 * pi
            x.adj[sel] <- -strwidth(x$tip.label)[sel] * 1.05
            sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * pi
            x.adj[sel] <- -strwidth(x$tip.label)[sel] * (2 * abs(XY$axe)[sel] / pi - 0.5)
            sel <- XY$axe > pi / 4 & XY$axe < 0.75 * pi
            y.adj[sel] <- strheight(x$tip.label)[sel] / 2
            sel <- XY$axe < -pi / 4 & XY$axe > -0.75 * pi
            y.adj[sel] <- -strheight(x$tip.label)[sel] * 0.75
            text(xx[1:nb.tip] + x.adj, yy[1:nb.tip] + y.adj, x$tip.label,
                 adj = c(adj, 0), font = font, srt = srt, cex = cex, col = tip.color)
        } else { # if lab4ut == "axial"
            adj <- as.numeric(abs(XY$axe) > pi/2)
            srt <- 180*XY$axe/pi
            srt[as.logical(adj)] <- srt[as.logical(adj)] - 180
            ## <FIXME> temporary check of the values of `srt':
            ## set to 0 if "-0.000001 < srt < 0"
            sel <- srt > -1e-6 & srt < 0
            if (any(sel)) srt[sel] <- 0
            ## </FIXME>
            ## `srt' takes only a single value, so we cannot vectorize this:
            for (i in 1:nb.tip)
              text(xx[i], yy[i], cex = cex, x$tip.label[i], adj = adj[i],
                   font = font, srt = srt[i], col = tip.color[i])
        }
    }
    if (type == "radial") {
        angle <- acos(xx[1:nb.tip]) * 180 / pi
        s1 <- angle > 90 & yy[1:nb.tip] > 0
        s2 <- angle < 90 & yy[1:nb.tip] < 0
        s3 <- angle > 90 & yy[1:nb.tip] < 0
        angle[s1] <- angle[s1] + 180
        angle[s2] <- -angle[s2]
        angle[s3] <- 180 - angle[s3]
        adj <- numeric(nb.tip)
        adj[xx[1:nb.tip] < 0] <- 1
        ## `srt' takes only a single value, so we cannot vectorize this:
        for (i in 1:nb.tip)
          text(xx[i], yy[i], x$tip.label[i], font = font, cex = cex,
               srt = angle[i], adj = adj[i], col = tip.color[i])
    }
    if (show.node.label)
      text(xx[root:length(xx)] + label.offset, yy[root:length(yy)],
           x$node.label, adj = adj, font = font, srt = srt, cex = cex)
    ## <FIXME>
    ## Maybe not all elements in `L' are useful...
    ## </FIXME>
    L <- list(type = type, use.edge.length = use.edge.length,
              node.pos = node.pos, show.tip.label = show.tip.label,
              show.node.label = show.node.label,
              edge.color = edge.color, edge.width = edge.width,
              font = font, cex = cex, adj = adj, srt = srt,
              no.margin = no.margin, label.offset = label.offset,
              x.lim = x.lim, y.lim = y.lim, direction = direction,
              tip.color = tip.color, Ntip = nb.tip, Nnode = nb.node)
    .last_plot.phylo <<- c(L, list(xx = xx, yy = yy))
    invisible(L)
}

phylogram.plot <- function(edge, nb.tip, nb.node, xx, yy,
                           horizontal, edge.color, edge.width)
{
    nodes <- (nb.tip + 1):(nb.tip + nb.node)
    if (!horizontal) {
        tmp <- yy
        yy <- xx
        xx <- tmp
    }
    ## un trait vertical à chaque noeud...
    x0v <- xx[nodes]
    y0v <- y1v <- numeric(nb.node)
    for (i in nodes) {
        j <- edge[which(edge[, 1] == i), 2]
        y0v[i - nb.tip] <- min(yy[j])
        y1v[i - nb.tip] <- max(yy[j])
    }
    ## ... et un trait horizontal partant de chaque tip et chaque noeud
    ##  vers la racine
    sq <- if (nb.node == 1) 1:nb.tip else c(1:nb.tip, nodes[-1])
    y0h <- yy[sq]
    x1h <- xx[sq]
    ## match() is very useful here becoz each element in edge[, 2] is
    ## unique (not sure this is so useful in edge[, 1]; needs to be checked)
    ## `pos' gives for each element in `sq' its index in edge[, 2]
    pos <- match(sq, edge[, 2])
    x0h <- xx[edge[pos, 1]]
    ## we need to reorder `edge.color' and `edge.width':
    ## <FIXME> Need to be fixed wrt to reordering pruningwise:
    edge.width <- edge.width[pos]
    edge.color <- edge.color[pos]
    ## </FIXME>
    e.w <- unique(edge.width)
    if (length(e.w) == 1) width.v <- rep(e.w, nb.node)
    else {
        width.v <- rep(1, nb.node)
        for (i in 1:nb.node) {
            br <- edge[which(edge[, 1] == i + nb.tip), 2]
            width <- unique(edge.width[br])
            if (length(width) == 1) width.v[i] <- width
        }
    }
    e.c <- unique(edge.color)
    if (length(e.c) == 1) color.v <- rep(e.c, nb.node)
    else {
        color.v <- rep("black", nb.node)
        for (i in 1:nb.node) {
            br <- edge[which(edge[, 1] == i + nb.tip), 2]
            color <- unique(edge.color[br])
            if (length(color) == 1) color.v[i] <- color
        }
    }
    if (horizontal) {
        segments(x0v, y0v, x0v, y1v, col = color.v, lwd = width.v) # draws vertical lines
        segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width) # draws horizontal lines
    } else {
        segments(y0v, x0v, y1v, x0v, col = color.v, lwd = width.v) # draws horizontal lines
        segments(y0h, x0h, y0h, x1h, col = edge.color, lwd = edge.width) # draws vertical lines
    }
}

cladogram.plot <- function(edge, xx, yy, edge.color, edge.width)
  segments(xx[edge[, 1]], yy[edge[, 1]], xx[edge[, 2]], yy[edge[, 2]],
           col = edge.color, lwd = edge.width)

unrooted.xy <- function(nb.tip, nb.node, edge, edge.length)
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
          if (i > nb.tip) foo(i, angle[i], axis[i])
    }
    root <- nb.tip + 1
    nb.edge <- dim(edge)[1]
    yy <- xx <- numeric(nb.tip + nb.node)
    nb.sp <- .C("node_depth", as.integer(nb.tip), as.integer(nb.node),
                as.integer(edge[, 1]), as.integer(edge[, 2]),
                as.integer(nb.edge), double(nb.tip + nb.node),
                DUP = FALSE, PACKAGE = "ape")[[6]]
    ## `angle': the angle allocated to each node wrt their nb of tips
    ## `axis': the axis of each branch
    axis <- angle <- numeric(nb.tip + nb.node)
    ## start with the root...
    ## xx[root] <- yy[root] <- 0 # already set!
    foo(root, 2*pi, 0)

    M <- cbind(xx, yy)
    axe <- axis[1:nb.tip] # the axis of the terminal branches (for export)
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

plot.multi.tree <- function(x, layout = 1, ...)
{
    if (layout > 1)
      layout(matrix(1:layout, ceiling(sqrt(layout)), byrow = TRUE))
    if (!par("ask")) {
        par(ask = TRUE)
        on.exit(par(ask = FALSE))
    }
    for (i in x) plot(i, ...)
}

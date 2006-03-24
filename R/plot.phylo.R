### plot.phylo.R (2006-01-06)
###
###          Plot Phylogenies
###
### Copyright 2002-2006 Emmanuel Paradis
###
### This file is part of the `ape' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

plot.phylo <- function(x, type = "phylogram", use.edge.length = TRUE,
                       node.pos = NULL, show.tip.label = TRUE,
                       show.node.label = FALSE, edge.color = "black",
                       edge.width = 1, font = 3,cex = par("cex"),
                       adj = NULL, srt = 0, no.margin = FALSE,
                       root.edge = FALSE, label.offset = 0, underscore = FALSE,
                       x.lim = NULL, y.lim = NULL, direction = "rightwards",
                       lab4ut = "horizontal", ...)
{
    type <- match.arg(type, c("phylogram", "cladogram", "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards",
                                        "upwards", "downwards"))
    if (is.null(x$edge.length)) use.edge.length <- FALSE
    tmp <- as.numeric(x$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    if (!show.tip.label) x$tip.label <- rep("", nb.tip)
    if (type == "unrooted" || !use.edge.length) root.edge <- FALSE
    if (type %in% c("phylogram", "cladogram")) {
        if (is.null(node.pos)) {
            if (type == "phylogram") node.pos <- 1
            if (type == "cladogram") node.pos <- if (!use.edge.length) 2 else 1
        }
        yy <- if (node.pos == 1) node.height(x$edge) else node.height.clado(x$edge)
        if (!use.edge.length) {
            xx <- node.depth(x$edge) - 1
            xx <- max(xx) - xx
        } else  {
            nms <- c(-(1:nb.node), 1:nb.tip)
            xx <- .C("node_depth_edgelength", as.integer(nb.tip),
                     as.integer(nb.node), as.integer(x$edge[, 1]),
                     as.integer(x$edge[, 2]), as.integer(nms),
                     as.double(x$edge.length), as.double(numeric(nb.tip + nb.node)),
                     DUP = FALSE, PACKAGE = "ape")[[7]]
            names(xx) <- as.character(nms)
        }
    }
    if (type == "unrooted") {
        XY <- if (use.edge.length) unrooted.xy(nb.tip, nb.node, x$edge, x$edge.length) else unrooted.xy(nb.tip, nb.node, x$edge, rep(1, dim(x$edge)[1]))
        ## rescale so that we have only positive values
        xx <- XY$M[, 1] - min(XY$M[, 1])
        yy <- XY$M[, 2] - min(XY$M[, 2])
    }
    if (type == "radial") {
        X <- node.depth(x$edge)
        X[X == 1] <- 0
        ## radius:
        X <- 1 - X / nb.tip
        ## angle:
        Y <- node.height(x$edge) * 2 * pi / nb.tip
        xx <- X * cos(Y)
        yy <- X * sin(Y)
    }
    if (type %in% c("phylogram", "cladogram") && direction != "rightwards") {
        if (direction == "leftwards") {
            xx <- -xx
            xx <- xx - min(xx)
        }
        if (direction %in% c("upwards", "downwards")) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
            if (direction == "downwards") {
                yy <- -yy
                yy <- yy - min(yy)
            }
        }
    }
    if (type %in% c("phylogram", "cladogram") && root.edge) {
        if (direction == "rightwards") xx <- xx + x$root.edge
        if (direction == "upwards") yy <- yy + x$root.edge
    }
    if (no.margin) par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (type %in% c("phylogram", "cladogram")) {
            if (direction %in% c("rightwards", "leftwards")) {
                x.lim <- c(0, NA)
                x.lim[2] <- if (direction == "leftwards")
                  max(xx["-1"] + nchar(x$tip.label) * 0.018 * max(xx) * cex)
                else  max(xx[as.character(1:nb.tip)] +
                          nchar(x$tip.label) * 0.018 * max(xx) * cex)
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
            if (type %in% c("phylogram", "cladogram"))
              x.lim <- if (direction %in% c("rightwards", "leftwards"))
                c(0, x.lim) else c(1, x.lim)
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
        if (type %in% c("phylogram", "cladogram")) {
            if (direction %in% c("upwards", "downwards")) {
                y.lim <- c(0, NA)
                y.lim[2] <- if (direction == "downwards")
                  max(yy["-1"] + nchar(x$tip.label) * 0.018 * max(yy) * cex)
                else max(yy[as.character(1:nb.tip)] +
                         nchar(x$tip.label) * 0.018 * max(yy) * cex)
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
            if(type %in% c("phylogram", "cladogram"))
              y.lim <- if (direction %in% c("rightwards", "leftwards"))
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
    if (type %in% c("phylogram", "cladogram") && root.edge) {
        if (direction == "leftwards") x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") y.lim[2] <- y.lim[2] + x$root.edge
    }

    edge.color <- rep(edge.color, length.out = dim(x$edge)[1])
    edge.width <- rep(edge.width, length.out = dim(x$edge)[1])
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "",
         ylab = "", xaxt = "n", yaxt = "n", bty = "n", ...)
    if (is.null(adj))
      adj <-
        if (type %in% c("phylogram", "cladogram") && direction == "leftwards") 1
        else 0
    if (type %in% c("phylogram", "cladogram")) {
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
        if (direction %in% c("upwards", "downwards")) {
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
                     direction, edge.color, edge.width)
    else cladogram.plot(x$edge, xx, yy, edge.color, edge.width)
    if (root.edge)
      switch(direction, "rightwards" = segments(0, yy["-1"], x$root.edge, yy["-1"]),
             "leftwards" = segments(xx["-1"], yy["-1"], xx["-1"] + x$root.edge, yy["-1"]),
             "upwards" = segments(xx["-1"], 0, xx["-1"], x$root.edge),
             "downwards" = segments(xx["-1"], yy["-1"], xx["-1"], yy["-1"] + x$root.edge))
    if (!underscore) x$tip.label <- gsub("_", " ", x$tip.label)
    if (type %in% c("phylogram", "cladogram")) {
        text(xx[as.character(1:nb.tip)] + lox, yy[as.character(1:nb.tip)] + loy,
             x$tip.label, adj = adj, font = font, srt = srt, cex = cex)
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
            text(xx[as.character(1:nb.tip)] + x.adj, yy[as.character(1:nb.tip)] + y.adj,
                 x$tip.label, adj = c(adj, 0), font = font, srt = srt, cex = cex)
        } else { # if lab4ut == "axial"
            adj <- as.numeric(abs(XY$axe) > pi / 2)
            srt <- 180 * XY$axe / pi
            srt[as.logical(adj)] <- srt[as.logical(adj)] - 180
            for (i in 1:nb.tip)
              text(xx[as.character(i)], yy[as.character(i)], cex = cex,
                   x$tip.label[i], adj = adj[i], font = font, srt = srt[i])
        }
    }
    if (type == "radial") {
        angle <- acos(xx[as.character(1:nb.tip)]) * 180 / pi
        s1 <- angle > 90 & yy[as.character(1:nb.tip)] > 0
        s2 <- angle < 90 & yy[as.character(1:nb.tip)] < 0
        s3 <- angle > 90 & yy[as.character(1:nb.tip)] < 0
        angle[s1] <- angle[s1] + 180
        angle[s2] <- -angle[s2]
        angle[s3] <- 180 - angle[s3]
        adj <- numeric(nb.tip)
        adj[xx[as.character(1:nb.tip)] < 0] <- 1
        for (i in 1:nb.tip)
          text(xx[as.character(i)], yy[as.character(i)], x$tip.label[i],
               font = font, cex = cex, srt = angle[i], adj = adj[i])
    }
    if (show.node.label) text(xx[as.character(-(1:nb.node))] + label.offset,
                              yy[as.character(-(1:nb.node))], x$node.label,
                              adj = adj, font = font, srt = srt, cex = cex)
    L <- list(type = type, use.edge.length = use.edge.length,
              node.pos = node.pos, show.tip.label = show.tip.label,
              show.node.label = show.node.label,
              edge.color = edge.color, edge.width = edge.width,
              font = font, cex = cex, adj = adj, srt = srt,
              no.margin = no.margin, label.offset = label.offset,
              x.lim = x.lim, y.lim = y.lim, direction = direction)
    .last_plot.phylo <<- c(L, list(xx = xx, yy = yy))
    invisible(L)
}

phylogram.plot <- function(edge, nb.tip, nb.node, xx, yy,
                           direction, edge.color, edge.width)
{
    if (direction %in% c("upwards", "downwards")) {
        tmp <- yy
        yy <- xx
        xx <- tmp
    }
    ## un trait vertical à chaque noeud...
    x0v <- xx[1:nb.node]
    y0v <- y1v <- numeric(nb.node)
    for (i in names(x0v)) {
        j <- edge[which(edge[, 1] == i), 2]
        y0v[-as.numeric(i)] <- min(yy[j])
        y1v[-as.numeric(i)] <- max(yy[j])
    }
    ## ... et un trait horizontal partant de chaque tip et chaque noeud
    ##  vers la racine
    sq <- if (nb.node == 1) 1:nb.tip else c(1:nb.tip, -(nb.node:2))
    sq <- as.character(sq)
    y0h <- yy[as.character(sq)]
    x1h <- xx[as.character(sq)]
    x0h <- numeric(nb.tip + nb.node - 1)
    j <- 1
    for (i in sq) {
        x0h[j] <- xx[edge[which(edge[, 2] == i), 1]]
        j <- j + 1
    }
    ## we need to reorder "edge.color" and "edge.width":
    names(edge.width) <- names(edge.color) <- edge[, 2]
    edge.width <- edge.width[sq]
    edge.color <- edge.color[sq]
    color.v <- rep("black", length(x0v))
    width.v <- rep(1, length(x0v))
    for (i in 1:length(x0v)) {
        br <- edge[which(edge[, 1] == names(x0v)[i]), 2]
        color <- unique(edge.color[br])
        width <- unique(edge.width[br])
        if (length(color) == 1) color.v[i] <- color
        if (length(width) == 1) width.v[i] <- width
    }
    names(x0v) <- NULL
    if (direction %in% c("rightwards", "leftwards")) {
        segments(x0v, y0v, x0v, y1v, col = color.v, lwd = width.v) # draws vertical lines
        segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width) # draws horizontal lines
    } else {
        segments(y0v, x0v, y1v, x0v, col = color.v, lwd = width.v) # draws horizontal lines
        segments(y0h, x0h, y0h, x1h, col = edge.color, lwd = edge.width) # draws vertical lines
    }
}

cladogram.plot <- function(edge, xx, yy, edge.color, edge.width)
{
    segments(xx[edge[, 1]], yy[edge[, 1]], xx[edge[, 2]], yy[edge[, 2]],
             col = edge.color, lwd = edge.width)
}

node.depth <- function(x)
### This function returns the depth of nodes and tips given by
### the number of descendants (1 is returned for tips).
### The input is the matrix `edge' of an object of class "phylo".
{
    tmp <- as.numeric(x)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    xx[(nb.node + 1):(nb.tip + nb.node)] <- 1
    ## `unused' says if the node or tip has NOT been used to compute
    ## the `xx' value of its ancestor
    unused <- rep(TRUE, nb.tip + nb.node)
    names(unused) <- names(xx)

    while(sum(unused) > 1) {
        term <- names(xx[!is.na(xx) & unused])
        ind <- x[, 2] %in% term
        term.br <- matrix(x[ind], length(term), 2)
        ## extract the nodes with at least 2 branches above
        basal <- names(which(table(term.br[, 1]) >= 2))
        for (nod in basal) {
            pair.ind <- which(x[, 1] == nod)
            pairs <- x[pair.ind, 2]
            ## Here we need to check that all the branches found in the next
            ## few lines just above are `available' for `clustering'; this may
            ## not be the case if other sister-branches have daughter-branches
            ## which are not yet defined, for instance if there is a multichotomy.
            if (all(pairs %in% term)) {
                xx[nod] <- sum(xx[pairs])
                unused[pairs] <- FALSE
            }
        }
    }
    xx
}

node.height <- function(x)
### This function returns the `height' of nodes and tips given by
### the average of the heights of the direct descendants. The tips
### have the heights 1 to the number of tips. This is not intended
### to be normally used directly by the user.
### The input is the matrix `edge' of an object of class "phylo".
{
    tmp <- as.numeric(x)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    ## yy: vecteur donnant l'ordonnée des lignes horizontales partant des tips
    ##     et des noeuds et allant vers la racine
    yy <- as.numeric(rep(NA, nb.tip + nb.node))
    names(yy) <- as.character(c(-(1:nb.node), 1:nb.tip))
    yy[(nb.node + 1):(nb.tip + nb.node)] <- 1:nb.tip
    ## `unused' says if the node or tip has NOT been used to compute
    ## the `yy' value of its ancestor
    unused <- rep(TRUE, nb.tip + nb.node)
    names(unused) <- names(yy)

    while(sum(unused) > 1) {
        term <- names(yy[!is.na(yy) & unused])
        ind <- x[, 2] %in% term
        term.br <- matrix(x[ind], length(term), 2)

        ## extract the nodes with at least 2 branches above
        basal <- names(which(table(term.br[, 1]) >= 2))
        for (nod in basal) {
            pair.ind <- which(x[, 1] == nod)
            pairs <- x[pair.ind, 2]
            ## Here we need to check that all the branches found in the next
            ## few lines just above are `available' for `clustering'; this may
            ## not be the case if other sister-branches have daughter-branches
            ## which are not yet defined, for instance if there is a multichotomy.
            if (all(pairs %in% term)) {
                yy[nod] <- sum(yy[pairs])/length(yy[pairs])
                unused[pairs] <- FALSE
            }
        }
    }
    yy
}

node.height.clado <- function(x)
### This function returns the `height' of nodes and tips given by
### the average of the heights of the tips. This is
### most suitable for cladogram with no edge lengths. The tips
### have the heights 1 to the number of tips. This is not intended
### to be normally used directly by the user.
### The input is the matrix `edge' of an object of class "phylo".
{
    tmp <- as.numeric(x)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    ## yy: vecteur donnant l'ordonnée des lignes horizontales partant des tips
    ##     et des noeuds et allant vers la racine
    yy <- as.numeric(rep(NA, nb.tip + nb.node))
    names(yy) <- as.character(c(-(1:nb.node), 1:nb.tip))
    N <- yy ## to weigh with the number of tips
    yy[(nb.node + 1):(nb.tip + nb.node)] <- 1:nb.tip
    N[(nb.node + 1):(nb.tip + nb.node)] <- 1
    ## `unused' says if the node or tip has NOT been used to compute
    ## the `yy' value of its ancestor
    unused <- rep(TRUE, nb.tip + nb.node)
    names(unused) <- names(yy)

    while(sum(unused) > 1) {
        term <- names(yy[!is.na(yy) & unused])
        ind <- x[, 2] %in% term
        term.br <- matrix(x[ind], length(term), 2)

        ## extract the nodes with at least 2 branches above
        basal <- names(which(table(term.br[, 1]) >= 2))
        for (nod in basal) {
            pair.ind <- which(x[, 1] == nod)
            pairs <- x[pair.ind, 2]
            ## Here we need to check that all the branches found in the next
            ## few lines just above are `available' for `clustering'; this may
            ## not be the case if other sister-branches have daughter-branches
            ## which are not yet defined, for instance if there is a multichotomy.
            if (all(pairs %in% term)) {
                N[nod] <- sum(N[pairs])
                yy[nod] <- sum(yy[pairs] * N[pairs]) / N[nod]
                unused[pairs] <- FALSE
            }
        }
    }
    yy
}

unrooted.xy <- function(nb.tip, nb.node, edge, edge.length)
{
    ## `xx' and `yy' are the coordinates of the nodes and tips
    xx <- as.numeric(rep(NA, (nb.tip + nb.node)))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    yy <- xx
    nb.sp <- node.depth(edge)[1:nb.node]
    angle <- as.numeric(rep(NA, nb.node)) # the angle allocated to each node wrt their nb of tips
    names(angle) <- names(nb.sp)
    axis <- angle # the axis of each branch
    axe <- numeric(nb.tip) # the axis of the terminal branches (for export)
    names(axe) <- as.character(1:nb.tip)
    ## start with the root...
    xx["-1"] <- yy["-1"] <- 0
    ind <- which(edge[, 1] == "-1")
    sons <- edge[ind, 2]
    alpha <- 2 * pi * nb.sp[sons] / nb.tip
    alpha[is.na(alpha)] <- 2 * pi / nb.tip # for tips
    start <- -alpha[1] / 2
    for (i in 1:length(sons)) {
        h <- edge.length[ind[i]]
        beta <- start + alpha[i] / 2
        xx[sons[i]] <- h * cos(beta)
        yy[sons[i]] <- h * sin(beta)
        if (as.numeric(sons[i]) < 0) {
            axis[sons[i]] <- beta
            angle[sons[i]] <- alpha[i]
        } else axe[sons[i]] <- beta
        start <- start + alpha[i]
    }
    next.node <- sons[as.numeric(sons) < 0]
    ## ... and go on with the other nodes and the tips
    while (length(next.node)) {
        future.next.node <- NULL
        for (ancestor in next.node) {
            ind <- which(edge[, 1] == ancestor)
            sons <- edge[ind, 2]
            ## angle where to start to draw the lines:
            start <- axis[ancestor] - angle[ancestor] / 2
            for (i in 1:length(sons)) {
                h <- edge.length[ind[i]]
                if (as.numeric(sons[i]) < 0) {
                    angle[sons[i]] <- angle[ancestor] * nb.sp[sons[i]] / nb.sp[ancestor]
                    beta <- axis[sons[i]] <- start + angle[sons[i]] / 2
                    start <- start + angle[sons[i]]
                } else {
                    beta <- start + angle[ancestor] / nb.sp[ancestor] / 2
                    start <- start + angle[ancestor] / nb.sp[ancestor]
                    axe[sons[i]] <- beta
                }
                xx[sons[i]] <- h * cos(beta) + xx[ancestor]
                yy[sons[i]] <- h * sin(beta) + yy[ancestor]
            }
            ## In the following line we avoid tips but a vector of mode character
            ## is returned, hence the evaluation in `while(...)'.
            future.next.node <- c(future.next.node, sons[as.numeric(sons) < 0])
        }
        next.node <- future.next.node
    }
    M <- matrix(c(xx, yy), ncol = 2)
    rownames(M) <- names(xx)
    axe[axe > pi] <- axe[axe > pi] - 2 * pi # insures that returned angles are in [-PI, +PI]
    list(M = M, axe = axe)
}

plot.multi.tree <- function(x, layout = 1, ...)
{
    if (layout > 1)
      layout(matrix(1:layout, ceiling(sqrt(layout)), byrow = TRUE))
    if (!par("ask")) {
        changedpar <- TRUE
        par(ask = TRUE)
    } else changedpar <- FALSE
    for (i in x) plot(i, ...)
    if (changedpar) par(ask = FALSE)
}

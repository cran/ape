### plot.phylo.R  (2004-02-04)
###
###     Plot Phylogenies
###
### Copyright 2003 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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
                       node.pos = NULL, show.node.label = FALSE,
                       edge.color = NULL, edge.width = NULL, font = 3,
                       adj = 0, srt = 0, no.margin = FALSE, root.edge = FALSE,
                       label.offset = NULL, underscore = FALSE, x.lim = NULL, ...)
{
    phy <- x
    rm(x)
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    type <- match.arg(type, c("phylogram", "cladogram", "unrooted"))
    if (is.null(phy$edge.length)) use.edge.length <- FALSE
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)

    if (type == "unrooted" | !use.edge.length) root.edge <- FALSE

    if (type %in% c("phylogram", "cladogram")) {
        if (is.null(node.pos)) {
            if (type == "phylogram") node.pos <- 1
            if (type == "cladogram") {
                if (!use.edge.length) node.pos <- 2 else node.pos <- 1
            }
        }
        if (node.pos == 1) yy <- node.height(phy$edge)
        if (node.pos == 2) yy <- node.height.clado(phy$edge)
        if (!use.edge.length) {
            xx <- node.depth(phy$edge) - 1
            xx <- max(xx) - xx
        } else {
            xx <- node.depth.edgelength(phy$edge, phy$edge.length)
        }
        if (root.edge) xx <- xx + phy$root.edge
    } else { # if type == "unrooted"
        if (use.edge.length) XY <- unrooted.xy(nb.tip, nb.node, phy$edge, phy$edge.length) else
        XY <- unrooted.xy(nb.tip, nb.node, phy$edge, rep(1, dim(phy$edge)[1]))
        ## rescale so that we have only positive values
        xx <- XY[, 1] - min(XY[, 1])
        yy <- XY[, 2] - min(XY[, 2])
    }

    if (no.margin) op <- par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (type %in% c("phylogram", "cladogram")) {
            x.lim <- max(xx[as.character(1:nb.tip)] +
                         nchar(phy$tip.label) * 0.018 * max(xx) * par("cex"))
        } else { # if type == "unrooted"
            offset <- max(nchar(phy$tip.label) * 0.018 * max(xx) * par("cex"))
            x.lim <- max(xx)
            y.lim <- max(yy)
        }
    }
    if (is.null(label.offset)) label.offset <- 0.1
    if (is.null(edge.color)) edge.color <- rep("black", dim(phy$edge)[1]) else {
        names(edge.color) <- phy$edge[, 2]
        edge.color <- edge.color[as.character(c(1:nb.tip, -(nb.node:2)))]
    }
    if (is.null(edge.width)) edge.width <- rep(1, dim(phy$edge)[1]) else {
        names(edge.width) <- phy$edge[, 2]
        edge.width <- edge.width[as.character(c(1:nb.tip, -(nb.node:2)))]
    }
    if (type %in% c("phylogram", "cladogram")) {
        plot(0, type="n", xlim=c(0, x.lim), ylim=c(1, nb.tip),
             xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
    } else { # if type == "unrooted"
        plot(0, type="n", xlim=c(0 - offset, x.lim + offset), ylim=c(0 - offset, y.lim + offset),
             xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
    }

    switch(type, "phylogram" = phylogram.plot(phy$edge, nb.tip, nb.node, xx, yy, edge.color, edge.width),
                 "cladogram" = cladogram.plot(phy$edge, xx, yy, edge.color, edge.width),
                 "unrooted" = unrooted.plot(phy$edge, xx, yy, edge.color, edge.width))

    if (root.edge) segments(0, yy["-1"], phy$root.edge, yy["-1"])

    if (!underscore) phy$tip.label <- gsub("_", " ", phy$tip.label)
    if (type %in% c("phylogram", "cladogram")) {
        text(xx[as.character(1:nb.tip)] + label.offset, 1:nb.tip, phy$tip.label,
             adj = adj, font = font, srt = srt)
    } else { # if type == "unrooted"
        text(xx[as.character(1:nb.tip)], yy[as.character(1:nb.tip)],
             phy$tip.label, adj = adj, font = font, srt = srt)
    }
    if (show.node.label) text(xx[as.character(-(1:nb.node))] + label.offset,
                              yy[as.character(-(1:nb.node))], phy$node.label,
                              adj = adj, font = font, srt = srt)
    invisible(list(type = type, use.edge.length = use.edge.length, node.pos = node.pos,
                   show.node.label = show.node.label, edge.color = edge.color,
                   edge.width = edge.width, font = font, adj = adj, srt = srt,
                   no.margin = no.margin, label.offset = label.offset, x.lim = x.lim))
}

phylogram.plot <- function(edge, nb.tip, nb.node, xx, yy, edge.color, edge.width)
{
    if (nb.node == 1) {
        x0v <- 0
        y0v <- 1
        y1v <- nb.tip
        x0h <- rep(0, nb.tip)
        x1h <- xx[as.character(1:nb.tip)]
        y0h <- 1:nb.tip
    } else {
        ## un trait vertical à chaque noeud...
        x0v <- xx[1:nb.node]
        y0v <- y1v <- numeric(nb.node)
        for (i in as.numeric(names(x0v))) {
            pair <- as.character(edge[which(edge[, 1] == i), 2])
            y0v[-i] <- min(yy[pair])
            y1v[-i] <- max(yy[pair])
        }
        names(x0v) <- NULL
        ## ... et un trait horizontal partant de chaque tip et chaque noeud
        ##  vers la racine
        x0h <- x1h <- numeric(nb.tip + nb.node - 1)
        y0h <- numeric(nb.tip + nb.node - 1)
        j <- 1
        for (i in c(1:nb.tip, -(nb.node:2))) {
            y0h[j] <- yy[as.character(i)]
            x0h[j] <- xx[as.character(edge[which(edge[, 2] == i), 1])]
            x1h[j] <- xx[as.character(i)]
            j <- j + 1
        }
    }
    segments(x0v, y0v, x0v, y1v) # draws vertical lines
    segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width) # draws horizontal lines
}

cladogram.plot <- function(edge, xx, yy, edge.color, edge.width)
{
    segments(xx[edge[, 1]], yy[edge[, 1]], xx[edge[, 2]], yy[edge[, 2]],
             col = edge.color, lwd = edge.width)
}
### `cladogram.plot' and `unrooted.plot' are currently identical (2003-11-23)
unrooted.plot <- function(edge, xx, yy, edge.color, edge.width)
{
    segments(xx[edge[, 1]], yy[edge[, 1]], xx[edge[, 2]], yy[edge[, 2]],
             col = edge.color, lwd = edge.width)
}

node.depth.edgelength <- function(x, el)
### This function returns the distance from the root to
### each node and tip with respect to the edge lengths.
### The input is the matrix `edge' of an object of class
### "phylo", and the corresponding vector edge.length.
{
    tmp <- as.numeric(x)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    ## xx: vecteur donnant la distance d'un noeud ou tip à partir de la racine
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(x[, 2] == nod)
        base <- x[ind, 1]
        xx[i] <- xx[base] + el[ind]
    }
    xx
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
        ind <- as.logical(match(x[, 2], term))
        ind[is.na(ind)] <- FALSE
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
        ind <- as.logical(match(x[, 2], term))
        ind[is.na(ind)] <- FALSE
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
        ind <- as.logical(match(x[, 2], term))
        ind[is.na(ind)] <- FALSE
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
    angle <- as.numeric(rep(NA, (nb.node))) # the angle allocated to each node wrt their nb of tips
    names(angle) <- names(nb.sp)
    axis <- angle # the axis of each branch
    ## start with the root...
    xx["-1"] <- yy["-1"] <- 0
    ind <- which(edge[, 1] == "-1")
    sons <- edge[ind, 2]
    alpha <- 2 * pi / length(sons)
    for (i in 1:length(sons)) {
        h <- edge.length[ind[i]]
        xx[sons[i]] <- h * cos(i * alpha)
        yy[sons[i]] <- h * sin(i * alpha)
        if (as.numeric(sons[i]) < 0) {
            axis[sons[i]] <- i * alpha
            angle[sons[i]] <- 2 * pi * nb.sp[sons[i]] / nb.tip
        }
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
    M
}

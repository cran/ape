### plot.phylo.R  (2003-02-04)
###
###     Plot Phylogenies
###
### Copyright 2002 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

plot.phylo <- function(x, show.node.label = FALSE, edge.color = NULL,
                       edge.width = NULL, font = 3, adj = 0, srt = 0,
                       no.margin = FALSE, label.offset = NULL,
                       x.lim = NULL, ...)
{
    phy <- x
    rm(x)
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)

    if (nb.node == 1) {
        x0v <- 0
        y0v <- 1
        y1v <- nb.tip
        x0h <- rep(0, nb.tip)
        x1h <- phy$edge.length
        y0h <- 1:nb.tip
        xx <- c(0, x1h)
        names(xx) <- as.character(c(-1, 1:nb.tip))
    } else {
        ## yy: vecteur donnant l'ordonnée des lignes horizontales partant des tips
        ##     et des noeuds et allant vers la racine
        yy <- as.numeric(rep(NA, nb.tip + nb.node))
        names(yy) <- as.character(c(-(1:nb.node), 1:nb.tip))
        xx <- yy # will be used later...
        yy[as.character(1:nb.tip)] <- 1:nb.tip
        ## `unused' says if the node or tip has NOT been used to compute
        ## the `yy' value of its ancestor
        unused <- rep(TRUE, nb.tip + nb.node)
        names(unused) <- names(yy)

        while(sum(unused) > 1) {
            term <- names(yy[!is.na(yy) & unused])
            ind <- as.logical(match(phy$edge[, 2], term))
            ind[is.na(ind)] <- FALSE
            term.br <- matrix(phy$edge[ind], length(term), 2)
   
            ## extract the nodes with at least 2 branches above
            basal <- names(which(table(term.br[, 1]) >= 2))
            for (nod in basal) {
                pair.ind <- which(phy$edge[, 1] == nod)
                pairs <- phy$edge[pair.ind, 2]
                ## Here we need to check that all the branches found in the next
                ## few lines just above are `available' for `clustering'; this may
                ## not be the case if other sister-branches have daughter-branches
                ## which are not yet defined, for instance if there is a multichotomy.
                ## This fixes a bug where trees with multichotomies at intermediate
                ## levels (i.e. neither at the root, nor terminal) where not plotted
                ## EP (26-09-2002)
                if (all(pairs %in% term)) {
                    yy[nod] <- sum(yy[pairs])/length(yy[pairs])
                    unused[pairs] <- FALSE
                }
            }
        }

        ## xx: vecteur donnant la distance d'un noeud ou tip à partir de la racine
        xx["-1"] <- 0
        for (i in 2:length(xx)) {
            nod <- names(xx[i])
            ind <- which(phy$edge[, 2] == nod)
            base <- phy$edge[ind, 1]
            xx[i] <- xx[base] + phy$edge.length[ind]
        }

        ## un trait vertical à chaque noeud...
        x0v <- xx[1:nb.node]
        y0v <- y1v <- numeric(nb.node)
        for (i in as.numeric(names(x0v))) {
            pair <- as.character(phy$edge[which(phy$edge[, 1] == i), 2])
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
            x0h[j] <- xx[as.character(phy$edge[which(phy$edge[, 2] == i), 1])]
            x1h[j] <- xx[as.character(i)]
            j <- j + 1
        }
    }
    
    if (no.margin) op <- par(mai = rep(0, 4))
    if (is.null(x.lim)) x.lim <- max(xx[as.character(1:nb.tip)] +
                                     nchar(phy$tip.label) * 0.018 * max(xx) * par("cex"))
    if (is.null(label.offset)) label.offset <- 0.1
    if (is.null(edge.color)) edge.color <- rep("black", length(phy$edge.length))
    else {
        names(edge.color) <- phy$edge[, 2]
        edge.color <- edge.color[as.character(c(1:nb.tip, -(nb.node:2)))]
    }
    if (is.null(edge.width)) edge.width <- rep(1, length(phy$edge.length))
    else {
        names(edge.width) <- phy$edge[, 2]
        edge.width <- edge.width[as.character(c(1:nb.tip, -(nb.node:2)))]
    }
    plot(0, type="n", xlim=c(0, x.lim), ylim=c(1, nb.tip),
         xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
    segments(x0v, y0v, x0v, y1v) # draws vertical lines
    segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width) # draws horizontal lines
    text(xx[as.character(1:nb.tip)] + label.offset, 1:nb.tip, phy$tip.label,
         adj = adj, font = font, srt = srt)
    if (show.node.label) text(xx[as.character(-(1:nb.node))] + label.offset,
                              yy[as.character(-(1:nb.node))], phy$node.label,
                              adj = adj, font = font, srt = srt)
    invisible(list(show.node.label = show.node.label, edge.color = edge.color,
                   edge.width = edge.width, font = font, adj = adj, srt = srt,
                   no.margin = no.margin, label.offset = label.offset, x.lim = x.lim))
}

### drop.tip.R  (2005-04-14)
###
###     Remove Tips in a Phylogenetic Tree
###
### Copyright 2005 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

drop.tip <- function(phy, tip, trim.internal = TRUE, subtree = FALSE,
                     root.edge = 0)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    if (subtree) {
        trim.internal <- TRUE
        edge.bak <- phy$edge
    }
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    nobr <- is.null(phy$edge.length)
    if (is.numeric(tip)) tip <- phy$tip.label[tip]
    ## find the tips to drop...:
    del <- phy$tip.label %in% tip
    ## ... and the corresponding terminal branches:
    ind <- which(phy$edge[, 2] %in% as.character(which(del)))
    ## drop them...:
    phy$edge <- phy$edge[-ind, ]
    ## ... and the lengths if applies:
    if (!nobr) phy$edge.length <- phy$edge.length[-ind]
    ## drop the tip labels:
    phy$tip.label <- phy$tip.label[!del]
    if (trim.internal) {
        if (root.edge) {
            ## find the MRCA of the remaining tips:
            seq.nod <- list()
            ## This is modified since some tips were deleted!!
            for (i in phy$edge[, 2][as.numeric(phy$edge[, 2]) > 0]) {
                vec <- i
                j <- i
                while (j != "-1") {
                    ind <- which(phy$edge[, 2] == j)
                    j <- phy$edge[ind, 1]
                    vec <- c(vec, j)
                }
                seq.nod[[i]] <- vec
            }
            sn <- lapply(seq.nod, rev)
            i <- 1
            x <- unlist(lapply(sn, function(x) x[i]))
            while (length(unique(x)) == 1) {
                x <- unlist(lapply(sn, function(x) x[i]))
                i <-  i + 1
            }
            MRCA <- sn[[1]][i - 2]
            newrootedge <- if (is.null(phy$root.edge)) 0 else phy$root.edge
            for (i in 1:root.edge) {
                ind <- which(phy$edge[, 2] == MRCA)
                newrootedge <- newrootedge + phy$edge.length[ind]
                MRCA <- phy$edge[ind, 1]
                if (MRCA == "-1" & i < root.edge) {
                    newrootedge <- newrootedge
                    break
                }
            }
            phy$root.edge <- newrootedge
        } else {
            if (!is.null(phy$root.edge)) phy$root.edge <- NULL
        }
        while (!all(phy$edge[, 2][as.numeric(phy$edge[, 2]) < 0] %in% phy$edge[, 1])) {
            temp <- phy$edge[, 2][as.numeric(phy$edge[, 2]) < 0]
            k <- temp %in% phy$edge[, 1]
            ind <- phy$edge[, 2] %in% temp[!k]
            phy$edge <- phy$edge[!ind, ]
            if (!nobr) phy$edge.length <- phy$edge.length[!ind]
        }
    } else {
        temp <- phy$edge[, 2][as.numeric(phy$edge[, 2]) < 0]
        k <- temp %in% phy$edge[, 1]
        ind <- phy$edge[, 2] %in% temp[!k]
        phy$edge[which(ind), 2] <- as.character(nb.tip + (1:sum(ind)))
        if (is.null(phy$node.label)) new.tip.label <- rep("NA", sum(ind)) else {
            new.tip.label <- phy$node.label[!k]
            phy$node.label <- phy$node.label[k]
        }
        phy$tip.label <- c(phy$tip.label, new.tip.label)
    }
    useless.nodes <- names(which(table(phy$edge[, 1]) == 1))
    if (subtree) {
        if (!nobr) mnbr <- mean(phy$edge.length)
        if (length(useless.nodes) == 1) n <- length(tip) else {
            seq.nod <- list()
            wh <- numeric(0)
            for (i in as.character(which(del))) { # it is not needed to loop through all tips!
                vec <- i
                j <- i
                while (!(j %in% useless.nodes)) {
                    ind <- which(edge.bak[, 2] == j)
                    wh <- c(wh, ind)
                    j <- edge.bak[ind, 1]
                    vec <- c(vec, j)
                }
                seq.nod[[i]] <- vec
            }
            n <- table(unlist(lapply(seq.nod, function(x) rev(x)[1])))
        }
        new.lab <- paste("[", n, "_tips]", sep = "")
        for (i in 1:length(useless.nodes)) {
            wh <- which(phy$edge[, 1] == useless.nodes[i])
            phy$tip.label <- c(phy$tip.label, new.lab[i])
            if (wh == dim(phy$edge)[1]) {
                phy$edge <- rbind(phy$edge, c(useless.nodes[i], as.character(nb.tip + i)))
                if (!nobr) phy$edge.length <- c(phy$edge.length, mnbr)
            } else {
                phy$edge <- rbind(phy$edge[1:wh, ],
                                  c(useless.nodes[i], as.character(nb.tip + i)),
                                  phy$edge[(wh + 1):dim(phy$edge)[1], ])
                if (!nobr) phy$edge.length <- c(phy$edge.length[1:wh], mnbr,
                                                phy$edge.length[(wh + 1):(dim(phy$edge)[1] - 1)])
            }
        }
    } else {
        for (i in useless.nodes) {
            ind1 <- which(phy$edge[, 1] == i)
            ind2 <- which(phy$edge[, 2] == i)
            phy$edge[ind2, 2] <- phy$edge[ind1, 2]
            phy$edge <- phy$edge[-ind1, ]
            if (!nobr) {
                phy$edge.length[ind2] <- phy$edge.length[ind2] + phy$edge.length[ind1]
                phy$edge.length <- phy$edge.length[-ind1]
            }
        }
    }
    tmp <- as.numeric(phy$edge)
    n <- length(tmp)
    nodes <- tmp < 0
    ind.nodes <- (1:n)[nodes]
    ind.tips <- (1:n)[!nodes]
    new.nodes <- -as.numeric(factor(-tmp[nodes]))
    new.tips <- as.numeric(factor(tmp[!nodes]))
    tmp[ind.nodes] <- new.nodes
    tmp[ind.tips] <- new.tips
    dim(tmp) <- c(n / 2, 2)
    mode(tmp) <- "character"
    phy$edge <- tmp
    if (!trim.internal | subtree) {
        S <- write.tree(phy, multi.line = FALSE)
        phy <- if (nobr) clado.build(S) else tree.build(S)
    }
    phy
}

### drop.tip.R  (2003-01-20)
###
###     Remove Tips in a Phylogenetic Tree
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

drop.tip <- function(phy, tip, trim.internal = TRUE)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    nb.tip <- max(as.numeric(phy$edge))
    del <- phy$tip.label %in% tip
    ind <- which(phy$edge[, 2] %in% as.character(which(del)))
    phy$edge <- phy$edge[-ind, ]
    phy$edge.length <- phy$edge.length[-ind]
    phy$tip.label <- phy$tip.label[!del]
    if (trim.internal) {
        while (!all(phy$edge[, 2][as.numeric(phy$edge[, 2]) < 0] %in% phy$edge[, 1])) {
            temp <- phy$edge[, 2][as.numeric(phy$edge[, 2]) < 0]
            k <- temp %in% phy$edge[, 1]
            ind <- phy$edge[, 2] %in% temp[!k]
            phy$edge <- phy$edge[!ind, ]
            phy$edge.length <- phy$edge.length[!ind]
        }
    }
    else {
        temp <- phy$edge[, 2][as.numeric(phy$edge[, 2]) < 0]
        k <- temp %in% phy$edge[, 1]
        ind <- phy$edge[, 2] %in% temp[!k]
        phy$edge[which(ind), 2] <- as.character(nb.tip + (1:sum(ind)))
        if (is.null(phy$node.label)) new.tip.label <- rep("NA", sum(ind))
        else {
            new.tip.label <- phy$node.label[!k]
            phy$node.label <- phy$node.label[k]
        }
        phy$tip.label <- c(phy$tip.label, new.tip.label)
    }
    useless.nodes <- names(which(table(phy$edge[, 1]) == 1))
    for (i in useless.nodes) {
        ind1 <- which(phy$edge[, 1] == i)
        ind2 <- which(phy$edge[, 2] == i)
        phy$edge[ind2, 2] <- phy$edge[ind1, 2]
        phy$edge <- phy$edge[-ind1, ]
        phy$edge.length[ind2] <- phy$edge.length[ind2] + phy$edge.length[ind1]
        phy$edge.length <- phy$edge.length[-ind1]
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
    if (!trim.internal) phy <- tree.build(write.tree(phy))
    return(phy)
}

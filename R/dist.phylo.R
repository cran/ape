### dist.phylo.R  (2002-08-28)
###
###     Pairwise Distances from a Phylogenetic Tree
###
### Copyright 2002 Julien Claude <claude@isem.univ-montp2.fr>
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

dist.phylo <- function(phy)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)

    ## xx: vecteur donnant la distance d'un noeud ou tip à partir de la racine
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(phy$edge[, 2] == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    ## seq.nod (liste de vecteurs): séquence des noeuds allant de chaque tip vers la racine
    seq.nod <- list()
    for (i in as.character(1:nb.tip)) {
        vec <- i
        j <- i
        while (j != "-1") {
            ind <- which(phy$edge[, 2] == j)
            j <- phy$edge[ind, 1]
            vec <- c(vec, j)
        }
        seq.nod[[i]] <- vec
    }
    dist <- diag(rep(0, nb.tip))
    for (i in as.character(1:(nb.tip - 1))) {
        for (j in as.character(2:nb.tip)) {
            ind <- min(match(seq.nod[[i]], seq.nod[[j]]), na.rm = TRUE)
            k <- as.numeric(i)
            l <- as.numeric(j)
            dist[k, l] <- dist[l, k] <- xx[i] + xx[j] - 2 * xx[seq.nod[[j]][ind]]
        }
    }
    rownames(dist) <- colnames(dist) <- phy$tip.label
    return(dist)
}

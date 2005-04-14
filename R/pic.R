### pic.R  (2004-10-12)
###
###     Phylogenetically Independent Contrasts
###
### Copyright 2004 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

pic <- function(x, phy, scaled = TRUE, var.contrasts = FALSE)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    if (nb.node != nb.tip - 1) stop("\"phy\" is not fully dichotomous")
    if (length(x) != nb.tip) stop("length of phenotypic and of phylogenetic data do not match")

    phenotype <- as.numeric(rep(NA, nb.tip + nb.node))
    names(phenotype) <- as.character(c(1:nb.tip, -(1:nb.node)))
    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x
    } else {
        if(!any(is.na(match(names(x), phy$tip.label))))
          for (i in 1:nb.tip) phenotype[i] <- x[phy$tip.label[i]]
        else {
            phenotype[1:nb.tip] <- x
            warning("the names of argument \"x\" and the names of the tip labels
did not match: the former were ignored in the analysis.")
        }
    }

    bl <- phy$edge.length   # copy the branch lengths to rescale them subsequently
    ## `unused' says if the phenotype has NOT been used to compute a contrast
    unused <- rep(TRUE, nb.tip + nb.node)
    names(unused) <- names(phenotype)
    contr <- numeric(nb.node)
    names(contr) <- as.character(-(1:nb.node))
    if (var.contrasts) var.con <- contr

    while(sum(unused) > 1) {
        term <- names(phenotype[!is.na(phenotype) & unused])
        ind <- as.logical(match(phy$edge[, 2], term))
        ind[is.na(ind)] <- FALSE
        ## extract the nodes with 2 branches above
        basal <- names(which(table(phy$edge[ind, 1]) == 2))
        for (nod in basal) {
            pair.ind <- which(phy$edge[, 1] == nod)
            i <- pair.ind[1]
            j <- pair.ind[2]
            pair <- phy$edge[pair.ind, 2]
            a <- pair[1]
            b <- pair[2]
            contr[nod] <- if (scaled) (phenotype[a] - phenotype[b])/sqrt(bl[i] + bl[j]) else phenotype[a] - phenotype[b]
            if (var.contrasts) var.con[nod] <- bl[i] + bl[j]
            unused[pair] <- FALSE
            phenotype[nod] <- (phenotype[a] * bl[j] + phenotype[b] * bl[i]) / (bl[i] + bl[j])
            k <- which(phy$edge[, 2] == nod)
            bl[k] <- bl[k] + bl[i] * bl[j] / (bl[i] + bl[j])
        }
    }
    if (var.contrasts) return(cbind(contr, var.con)) else return(contr)
}

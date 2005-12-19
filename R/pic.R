### pic.R (2005-12-16)
###
###     Phylogenetically Independent Contrasts
###
### Copyright 2002-2005 Emmanuel Paradis
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
    if (class(phy) != "phylo")
      stop("object \"phy\" is not of class \"phylo\"")
    if (is.null(phy$edge.length))
      stop("your tree has no branch lengths: you may consider setting them equal to one, or using the function `compute.brlen'.")
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    if (nb.node != nb.tip - 1)
      stop("\"phy\" is not fully dichotomous")
    if (length(x) != nb.tip)
      stop("length of phenotypic and of phylogenetic data do not match")
    if (any(is.na(x)))
      stop("the present method cannot (yet) be used directly with missing data: you may consider removing the species with missing data from your tree with the function `drop.tip'.")

    phenotype <- as.numeric(rep(NA, nb.tip + nb.node))
    names(phenotype) <- as.character(c(1:nb.tip, -(1:nb.node)))
    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x
    } else {
        if(!any(is.na(match(names(x), phy$tip.label))))
          phenotype[1:nb.tip] <- x[phy$tip.label]
        else {
            phenotype[1:nb.tip] <- x
            warning('the names of argument "x" and the names of the tip labels
did not match: the former were ignored in the analysis.')
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
        ind <- phy$edge[, 2] %in% term
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

### rotate.R  (2004-09-22)
###
###     Rotate an Internal Branch of a Tree
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

rotate <- function(phy, group)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    if (length(group) == 1) if (group == "all") {
        ind <- which(as.numeric(phy$edge[, 2]) > 0)
        phy$edge[ind, 2] <- phy$edge[rev(ind), 2]
        phy$tip.label <- rev(phy$tip.label)
        return(phy)
    }
    if (is.character(group)) tip <- as.character(which(phy$tip.label == group))
    if (is.numeric(group)) tip <- as.character(group)
    ## Check that the group is monophyletic
    msg <- "the specified group is not monophyletic!"
    if (!all(diff(as.numeric(tip)) == 1)) stop(msg)
    ## Find the MRCA of the tips given as outgroup
    ## The following loop is `borrowed' from vcv.phylo()
    seq.nod <- list()
    nb.tip <- max(as.numeric(phy$edge))
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
    sn <- lapply(seq.nod[tip], rev)
    i <- 1
    x <- unlist(lapply(sn, function(x) x[i]))
    while (length(unique(x)) == 1) {
        x <- unlist(lapply(sn, function(x) x[i]))
        i <-  i + 1
    }
    MRCA <- sn[[1]][i - 2]
    ## Then check that all descendants of this node
    ## are included in the outgroup
    desc <- names(unlist(lapply(seq.nod, function(x) which(x == MRCA))))
    if (length(tip) != length(desc)) stop(msg)
    if (!all(sort(tip) == sort(desc))) stop(msg)
    ind <- which(phy$edge[, 2] %in% tip)
    phy$edge[ind, 2] <- phy$edge[rev(ind), 2]
    ind <- as.numeric(tip)
    phy$tip.label[ind] <- phy$tip.label[rev(ind)]
    phy
}

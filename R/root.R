### root.R  (2004-09-20)
###
###                (re)Roots Phylogenetic Trees
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

root <- function(phy, outgroup)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    if (is.character(outgroup)) tip <- as.character(which(phy$tip.label == outgroup))
    if (is.numeric(outgroup)) tip <- as.character(outgroup)
    ## First check that the outgroup is monophyletic--
    ## unless there's only one tip specified of course
    if (length(tip) > 1) {
        msg <- "the specified outgroup is not monophyletic!"
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
        newroot <- sn[[1]][i - 2]
        ## Then check that all descendants of this node
        ## are included in the outgroup
        desc <- names(unlist(lapply(seq.nod, function(x) which(x == newroot))))
        if (length(tip) != length(desc)) stop(msg)
        if (!all(sort(tip) == sort(desc))) stop(msg)

    } else {
        newroot <- phy$edge[which(phy$edge[, 2] == tip), 1]
    }
    if (newroot == "-1") return(phy)

    ## Invert all branches from the new root to the old one
    i <- which(phy$edge[, 2] == newroot)
    nod <- phy$edge[i, 1]
    while (nod != "-1") {
        j <- which(phy$edge[, 2] == nod)
        phy$edge[i, 1] <- phy$edge[i, 2]
        phy$edge[i, 2] <- nod
        i <- j
        nod <- phy$edge[i, 1]
    }

    i.oroot <- which(phy$edge[, 1] == "-1")    
    ## Unroot the tree if there's a basal dichotomy...
    if (length(i.oroot) == 2) {
        j <- i.oroot[which(i.oroot != i)]
        phy$edge[j, 1] <- phy$edge[i, 2]
        phy$edge <- phy$edge[-i, ]
        phy$edge.length[j] <- phy$edge.length[j] + phy$edge.length[i]
        phy$edge.length <- phy$edge.length[-i]
        phy$edge[which(phy$edge == newroot)] <- "-1"
    } else {
        ## ... otherwise just invert the root with the newroot
        phy$edge[which(phy$edge == newroot)] <- "-1"
        phy$edge[i.oroot] <- newroot
    }
    read.tree(text = write.tree(phy, multi.line = FALSE))
}

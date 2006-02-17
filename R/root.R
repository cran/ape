### root.R (2006-01-11)
###
###            Root of Phylogenetic Trees
###
### Copyright 2004-2006 Emmanuel Paradis
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

is.rooted <- function(phy)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    if (!is.null(phy$root.edge)) return(TRUE)
    else if (table(phy$edge[, 1])["-1"] > 2) return(FALSE) else return(TRUE)
}

unroot <- function(phy)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    if (dim(phy$edge)[1] == 2)
      stop("cannot unroot a tree with two edges.")
    ## delete FIRST the root.edge (in case this is sufficient to
    ## unroot the tree, i.e. there is multichotomy at the root)
    if (!is.null(phy$root.edge)) phy$root.edge <- NULL
    if (!is.rooted(phy)) return(phy)
    i.oroot <- which(phy$edge[, 1] == "-1")
    i <- i.oroot[1]
    if (as.numeric(phy$edge[i, 2]) > 0) i <- i.oroot[2]
    newroot <- phy$edge[i, 2]
    j <- i.oroot[which(i.oroot != i)]
    phy$edge[j, 1] <- newroot
    phy$edge <- phy$edge[-i, ]
    if (!is.null(phy$edge.length)) {
        phy$edge.length[j] <- phy$edge.length[j] + phy$edge.length[i]
        phy$edge.length <- phy$edge.length[-i]
    }
    phy$edge[which(phy$edge == newroot)] <- "-1"
    read.tree(text = write.tree(phy, multi.line = FALSE))
}

root <- function(phy, outgroup)
{
    if (class(phy) != "phylo") stop('object "phy" is not of class "phylo"')
    if (is.character(outgroup))
      tip <- as.character(which(phy$tip.label %in% outgroup))
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

    } else newroot <- phy$edge[which(phy$edge[, 2] == tip), 1]

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
        if (!is.null(phy$edge.length)) {
            phy$edge.length[j] <- phy$edge.length[j] + phy$edge.length[i]
            phy$edge.length <- phy$edge.length[-i]
        }
        phy$edge[which(phy$edge == newroot)] <- "-1"
    } else {
        ## ... otherwise just invert the root with the newroot
        phy$edge[which(phy$edge == newroot)] <- "-1"
        phy$edge[i.oroot] <- newroot
        ## ... and invert finally! (fixed 2005-11-07)
        phy$edge[i, ] <- rev(phy$edge[i, ])
    }
    if (!is.null(phy$node.label)) {
        tmp <- phy$node.label[1]
        phy$node.label[1] <- phy$node.label[-as.numeric(newroot)]
        phy$node.label[-as.numeric(newroot)] <- tmp
    }
    read.tree(text = write.tree(phy, multi.line = FALSE))
}

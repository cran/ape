### which.edge.R (2004-12-10)
###
###        Identifies Edges of a Tree
###
### Copyright 2004 Emmanuel Paradis
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

which.edge <- function(phy, group)
{
    if (class(phy) != "phylo")
      stop("object \"phy\" is not of class \"phylo\"")
    if (is.character(group))
      tip <- as.character(which(phy$tip.label == group))
    if (is.numeric(group))
      tip <- as.character(group)
    ## Find the MRCA of the tips given as group
    ## The following loop is `borrowed' from vcv.phylo()
    seq.nod <- list()
    wh <- numeric(0)
    for (i in tip) { # it is not needed to loop through all tips!
        vec <- i
        j <- i
        while (j != "-1") {
            ind <- which(phy$edge[, 2] == j)
            wh <- c(wh, ind)
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
    wh <- sort(unique(wh))
    if (MRCA != "-1") {
        root <- "-1"
        while (root != MRCA) {
            i <- which(phy$edge[wh, 1] == root)
            root <- phy$edge[wh[i], 2]
            wh <- wh[-i]
        }
    }
    wh
}

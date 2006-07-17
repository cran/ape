### collapse.singles.R  (2006-07-15)
###
###     Collapse "Single" Nodes
###
### Copyright 2006 Ben Bolker
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

collapse.singles <- function(tree)
{
    elen <- tree$edge.length
    xmat <- matrix(as.numeric(tree$edge), nrow = nrow(tree$edge),
                   ncol = ncol(tree$edge))
    singles <- NA
    while (length(singles) > 0) {
        tx <- table(xmat[xmat < 0])
        singles <- as.numeric(names(tx)[tx < 3])
        if (length(singles) > 0) {
            i <- singles[1]
            prev.node <- which(xmat[, 2] == i)
            next.node <- which(xmat[,1] == i)
            xmat[prev.node, 2] <- xmat[next.node, 2]
            xmat <- xmat[xmat[, 1] != i, ] ## drop
            xmat[xmat < i] <- xmat[xmat < i] + 1 ## adjust indices
            elen[prev.node] <- elen[prev.node] + elen[next.node]
            elen <- elen[-next.node]
        }
    }
    tree$edge <- matrix(as.character(xmat),
                        nrow = nrow(xmat), ncol = ncol(xmat))
    tree$edge.length <- elen
    tree
}

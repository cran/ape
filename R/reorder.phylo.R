### reorder.phylo.R  (2006-08-25)
###
###     Internal Reordering of Trees
###
### Copyright 2006 Emmanuel Paradis
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

### See below for a different, but less efficient version.

reorder.phylo <- function(x, order = "cladewise", ...)
{
    order <- match.arg(order, c("cladewise", "pruningwise"))
    if (order == "cladewise") {
        x <- read.tree(text = write.tree(x, multi.line = FALSE))
        attr(x, "order") <- "cladewise"
    }
    else {
        mode(x$edge) <- "numeric" # temporary
        nb.tip <- max(x$edge)
        nb.node <- -min(x$edge)
        nb.edge <- dim(x$edge)[1]
        neworder <- numeric(nb.edge)
        j <- 1
        ## `ready' indicates whether an edge is ready to be
        ## collected; only the terminal edges are initially ready
        ready <- x$edge[, 2] > 0

        while(j < nb.edge) {
            ## We look at the edges ready:
            sel <- which(ready)
            k <- 1 # index along `sel'
            nodes <- unique(x$edge[sel, 1])
            nn <- tabulate(-x$edge[sel, 1])
            for (i in nodes) {
                nbr <- nn[-i]
                ## If the node has only one descendant, no need
                ## to go further; we also check that all edges
                ## originating from these nodes are included:
                if (nbr == 1 || i %in% x$edge[-sel, 1]) {
                    k <- k + nbr
                    next
                }
                neworder[j:(j + nbr - 1)] <- tmp <- sel[k:(k + nbr - 1)]
                j <- j + nbr
                ready[tmp] <- FALSE
                ## The next one will work if the tree is in clade-wise
                ## order. Compared to:
                ##   ready[which(x$edge[, 2] == i)] <- TRUE
                ## the execution time is halved for 1000 tips.
                ready[sel[k] - 1] <- TRUE
                k <- k + nbr
            }
        }

        x$edge <- x$edge[neworder, ]
        if (!is.null(x$edge.length))
          x$edge.length <- x$edge.length[neworder]
        attr(x, "order") <- "pruningwise"
    }
    mode(x$edge) <- "character" # temporary
    x
}

### For memory, I include below the first version I wrote for
### this function. It is substantially slower (ca. 1.74 s with
### 1000 tips against ca. 0.35 s)

### reorder.phylo <- function(x, order = "cladewise", ...)
### {
###     order <- match.arg(order, c("cladewise", "pruningwise"))
###     if (order == "cladewise")
###       x <- read.tree(text = write.tree(x, multi.line = FALSE))
###     else {
###         mode(x$edge) <- "numeric" # temporary
###         nb.tip <- max(x$edge)
###         nb.node <- -min(x$edge)
###         nb.edge <- dim(x$edge)[1]
###         neworder <- numeric(nb.edge)
###         j <- 1
###         ## `ready' indicates whether an edge is ready to be
###         ## collected; only the terminal edges are initially ready
###         ready <- x$edge[, 2] > 0
###
###         while(j < nb.edge) {
###             ## We look at the edges ready:
###             sel <- which(ready)
###             for (i in sel) {
###                 ## Here we check that all edges originating
###                 ## from these nodes are included
###                 if (x$edge[i, 1] %in% x$edge[-sel, 1]) next
###                 neworder[j] <- i
###                 j <- j + 1
###                 ready[i] <- FALSE
###                 ready[which(x$edge[, 2] == x$edge[i, 1])] <- TRUE
###             }
###         }
###
###         x$edge <- x$edge[neworder, ]
###         if (!is.null(x$edge.length))
###           x$edge.length <- x$edge.length[neworder]
###     }
###     mode(x$edge) <- "character" # temporary
###     x
### }

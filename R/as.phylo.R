### as.phylo.R  (2004-11-25)
###
###           Conversion Among Tree Objects
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

as.phylo <- function (x, ...) UseMethod("as.phylo")

as.phylo.hclust <- function(x, ...)
{
    N <- dim(x$merge)[1]
    edge <- matrix(NA, 2 * N, 2)
    edge.length <- numeric(2 * N)
    ## `node' gives the number of the node for the i-th row of x$merge
    node <- numeric(N)
    node[N] <- -1
    cur.nod <- -2
    j <- 1
    for (i in N:1) {
        edge[j:(j + 1), 1] <- node[i]
        for (l in 1:2) {
            k <- j + l - 1
            if (x$merge[i, l] > 0) {
                edge[k, 2] <- node[x$merge[i, l]] <- cur.nod
                cur.nod <- cur.nod - 1
                edge.length[k] <- x$height[i] - x$height[x$merge[i, l]]
            } else {
                edge[k, 2] <- -x$merge[i, l]
                edge.length[k] <- x$height[i]
            }
        }
        j <- j + 2
    }
    mode(edge) <- "character"
    obj <- list(edge = edge, edge.length = edge.length, tip.label = x$labels)
    class(obj) <- "phylo"
    read.tree(text = write.tree(obj, multi.line = FALSE))
}

as.phylo.phylog <- function(x, ...)
{
    tr <- read.tree(text = x$tre)
    edge.length <- numeric(dim(tr$edge)[1])
    term  <- which(as.numeric(tr$edge[, 2]) > 0)
    inte  <- which(as.numeric(tr$edge[, 2]) < 0)
    edge.length[term] <- x$leaves[tr$tip.label]
    edge.length[inte] <- x$nodes[tr$node.label][-1]
    tr$edge.length <- edge.length
    if (x$nodes["Root"] != 0) {
        edge.root <- x$nodes["Root"]
        names(edge.root) <- NULL
        tr$edge.root <- edge.root
    }
    tr
}

as.hclust.phylo <- function(x, ...)
{
    if (!is.ultrametric(x)) stop("the tree is not ultrametric")
    if (!is.binary.tree(x)) stop("the tree is not binary")
    bt <- rev(branching.times(x))
    N <- length(bt)
    tmp <- x$edge
    mode(tmp) <- "numeric"
    nm <- as.numeric(names(bt))
    merge <- matrix(NA, N, 2)
    for (i in 1:N) {
        ind <- which(tmp[, 1] == nm[i])
        for (k in 1:2)
          merge[i, k] <- if (tmp[ind[k], 2] > 0) -tmp[ind[k], 2]
          else which(nm == tmp[ind[k], 2])
    }
    names(bt) <- NULL
    obj <- list(merge = merge, height = bt, order = 1:(N + 1),
                labels = x$tip.label, call = match.call(),
                method = "unknown")
    class(obj) <- "hclust"
    obj
}

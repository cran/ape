### bind.tree.R  (2004-08-31)
###
###     Bind Trees
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

bind.tree <- function(x, y, node = -1, branch = NULL, position = NULL)
{
    if (node >= 0) stop("node number must be a negative integer.")
    tmp <- as.numeric(x$edge)
    nb.node <- -min(tmp)
    nb.tip <- max(tmp)
    if(!is.null(branch)) {
        if (!is.character(branch)) branch <- as.character(branch)
        new.node <- as.character(-nb.node - 1)
        if (branch == "-1") {
            if (is.null(x$root.edge)) stop("there is no root edge.")
            if (position > x$root.edge) stop ("argument \"position\" is larger than the root edge.")
            x$edge[which(x$edge[, 1] == "-1"), 1] <- new.node
            x$edge <- rbind(c("-1", new.node), x$edge)
            x$root.edge <- x$root.edge - position
            x$edge.length <- c(position, x$edge.length)
            node <- "-1"
        }
        else {
            ind <- which(x$edge[, 2] == branch)
            if (position > x$edge.length[ind])
              stop ("argument \"position\" is larger than the length of the specified branch.")
            x$edge[ind, 2] <- new.node
            x$edge <- rbind(x$edge, c(new.node, branch))
            x$edge.length[ind] <- x$edge.length[ind] - position
            x$edge.length <- c(x$edge.length, position)
            node <- new.node
        }
        nb.node <- nb.node + 1
    }
    if (!is.character(node)) node <- as.character(node)
    tmp <- as.numeric(y$edge)
    tmp <- ifelse(tmp > 0, tmp + nb.tip, tmp - nb.node)    
    if (is.null(y$root.edge)) tmp[which(tmp == -(nb.node + 1))] <- node # implicit mode conversion
    else mode(tmp) <- "character"
    tmp <- matrix(tmp, ncol = 2)
    if (!is.null(y$root.edge)) {
        tmp <- rbind(c(node, as.character(-nb.node - 1)), tmp)
        y$edge.length <- c(y$root.edge, y$edge.length)
    }
    obj <- list(edge = rbind(x$edge, tmp),
                edge.length = c(x$edge.length, y$edge.length),
                tip.label = c(x$tip.label, y$tip.label))
    if (!is.null(x$node.label)) {
        obj$node.label <- if (!is.null(y$node.label)) c(x$node.label, y$node.label) else c(x$node.label, rep(NA, -min(as.numeric(y$edge))))
    }
    else if (!is.null(y$node.label)) obj$node.label <- c(rep(NA, -min(as.numeric(x$edge))), y$node.label)
    if (!is.null(x$root.edge)) obj$root.edge <- x$root.edge
    class(obj) <- "phylo"
    tree.build(write.tree(obj))
}

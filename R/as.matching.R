### as.phylo.R (2005-04-13)
###
###           Conversion Between Phylo and Matching Objects
###
### Copyright 2005 Emmanuel Paradis
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

as.matching <- function(x, ...) UseMethod("as.matching")

as.matching.phylo <- function(x, labels = TRUE, ...)
{
    if (!is.binary.tree(x)) stop("the tree must be dichotomous.")
    mode(x$edge) <- "numeric"
    nb.tip <- max(x$edge)
    nb.node <- nb.tip - 1
    newlabel <- c(1:nb.tip, rep(NA, nb.node))
    done <- rep(FALSE, nb.tip + nb.node)
    lab <- c(1:nb.tip, -(1:nb.node))
    names(newlabel) <- names(done) <- as.character(lab)
    nextnode <- nb.tip + 1
    obj <- matrix(NA, 0, 3)
    while(sum(done) < nb.tip + nb.node - 1) {
        ## find the nodes and tips already labelled and not yet done:
        sel <- x$edge[, 2] %in% lab[which(!is.na(newlabel) & !done)]
        ## find those that are sibling pairs:
        nod.sel <- names(which(table(x$edge[sel, 1]) == 2))
        ## select the appropriate lines in the edge matrix:
        tmp <- x$edge[, 1] %in% nod.sel
        ## To insure that the nodes are in the same order than
        ## their descendants (apparently, table() does not respect
        ## this order):
        nod.sel <- unique(x$edge[tmp, 1])
        ## Since the siblings are already matched by pairs,
        ## we can put them in a matrix like this:
        m <- matrix(as.numeric(x$edge[tmp, 2]),
                    ncol = 2, byrow = TRUE)
        ## we set 'done' before changing the labels
        done[as.character(m)] <- TRUE
        if (any(m < 0)) m[m < 0] <- newlabel[as.character(m[m < 0])]
        ## we insure that both elements in each row are sorted
        ## (do not forget to transpose the returned matrix!):
        m <- t(apply(m, 1, sort))
        ## finally we sort each row according to the first col:
        o <- order(m[, 1])
        m <- m[o, , drop = FALSE]
        ## we also reorder the ancestors
        ## (before changing for the new numbers):
        nod.sel <- nod.sel[o]
        the.new.nodes <- nextnode:(nextnode + nrow(m) - 1)
        m <- cbind(m, the.new.nodes)
        newlabel[as.character(nod.sel)] <- the.new.nodes
        nextnode <- nextnode + nrow(m)
        obj <- rbind(obj, m)
    }
    dimnames(obj) <- NULL
    obj <- list(matching = obj)
    if (!is.null(x$edge.length)) {
        edge.length <-
          x$edge.length[order(newlabel[as.character(x$edge[, 2])])]
        names(edge.length) <- NULL
        obj$edge.length <- edge.length
    }
    class(obj) <- "matching"
    if (labels) {
        obj$tip.label <- x$tip.label
        if (!is.null(x$node.label)) obj$node.label <- x$node.label
    }
    obj
}

as.phylo.matching <- function(x, ...)
{
    N <- 2 * dim(x$matching)[1]
    edge <- matrix(NA, N, 2)
    if (!is.null(x$edge.length)) {
        edge.length <- numeric(N)
        br <- TRUE
    }
    new.nodes <- numeric(N + 1)
    new.nodes[N + 1] <- -1
    nb.tip <- (N + 2)/2
    nb.node <- nb.tip - 1
    nextnode <- -2
    j <- 1
    for (i in nb.node:1) {
        edge[j:(j + 1), 1] <- new.nodes[x$matching[i, 3]]
        for (k in 1:2) {
            if (x$matching[i, k] > nb.tip) {
                edge[j + k - 1, 2] <- new.nodes[x$matching[i, k]] <- nextnode
                nextnode <- nextnode - 1
            } else edge[j + k - 1, 2] <- x$matching[i, k]
            if (br)
              edge.length[j + k - 1] <- x$edge.length[x$matching[i, k]]
        }
        j <- j + 2
    }
    mode(edge) <- "character"
    obj <- list(edge = edge)
    if (!is.null(x$tip.label)) obj$tip.label <- x$tip.label
    else obj$tip.label <- as.character(1:nb.tip)
    if (br) obj$edge.length <- edge.length
    class(obj) <- "phylo"
    read.tree(text = write.tree(obj, multi.line = FALSE))
}
